# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
"""
Job management for the PeriLab HTTP wrapper.

Each submitted simulation gets its own working directory under JOBS_DIR,
containing the uploaded input files. A background thread runs the PeriLab
binary with that directory as its cwd. PeriLab writes its own log file
(named after the input deck) inside that directory; the monitor threads
tail it both to report progress and to serve the /log endpoint. Job state
is kept in memory and mirrored to `job.json` in the job directory so a
restart can at least recover job history (running jobs are marked
"interrupted" since the process itself is gone).
"""

from __future__ import annotations

import io
import json
import os
import re
import shutil
import subprocess
import threading
import time
import uuid
from dataclasses import dataclass, field, fields
from enum import Enum
from pathlib import Path
from typing import Optional


class JobStatus(str, Enum):
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"  # was running when the server restarted


@dataclass
class Job:
    id: str
    dir: Path
    command: list
    status: JobStatus = JobStatus.QUEUED
    exit_code: Optional[int] = None
    error: Optional[str] = None
    created_at: float = field(default_factory=time.time)
    started_at: Optional[float] = None
    finished_at: Optional[float] = None
    user_id: Optional[str] = None

    # Live progress parsed from the PeriLab log while the job runs.
    # "[Info] Step: 1 / 2 [2.500e-04 s]" -> step=1, total_steps=2, sim_time=2.5e-04
    step: Optional[int] = None
    total_steps: Optional[int] = None
    sim_time: Optional[float] = None

    # not persisted / not serialized
    _process: Optional[subprocess.Popen] = field(default=None, repr=False, compare=False)

    @property
    def progress(self) -> Optional[float]:
        """Fraction of completion in [0, 1], or None if unknown/unbounded."""
        if self.total_steps and self.step is not None:
            return min(self.step / self.total_steps, 1.0)
        return None

    def to_dict(self) -> dict:
        d = {}
        for f in fields(self):
            if f.name == "_process":
                continue  # Popen; not serializable (holds a thread lock)
            d[f.name] = getattr(self, f.name)
        d["dir"] = str(self.dir)
        d["status"] = self.status.value
        if self.progress is not None:
            d["progress"] = self.progress
        return d

    def save(self) -> None:
        meta_path = self.dir / "job.json"
        try:
            meta_path.write_text(json.dumps(self.to_dict(), indent=2))
        except OSError:
            pass  # best-effort persistence only


def _safe_relative_path(name: str) -> Path:
    """Sanitize an uploaded filename: strip leading slashes and reject '..'."""
    p = Path(name.replace("\\", "/"))
    parts = [part for part in p.parts if part not in ("", ".", "/")]
    if any(part == ".." for part in parts):
        raise ValueError(f"invalid path in uploaded filename: {name!r}")
    if not parts:
        raise ValueError(f"empty filename")
    return Path(*parts)


class JobManager:
    def __init__(self, base_dir: Path, binary_path: str, max_concurrent: int = 2,
                 mpi_launcher: Optional[str] = None, retention_days: float = 7.0):
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)
        self.binary_path = binary_path
        self.mpi_launcher = mpi_launcher  # e.g. "mpiexecjl", or None if unavailable
        self.retention_days = retention_days
        self.jobs: dict[str, Job] = {}
        self._sem = threading.Semaphore(max_concurrent)
        self._registry_lock = threading.Lock()
        # Prunes working dirs of finished jobs once they age out of the
        # retention window; daemon so it never blocks shutdown.
        threading.Thread(target=self._sweep_finished, daemon=True).start()

    def _sweep_finished(self) -> None:
        """Delete working dirs of finished jobs past the retention window.

        Terminal statuses set ``finished_at`` before the process terminates,
        so a live job is never eligible. Reuses ``delete()`` for cleanup.
        """
        while True:
            now = time.time()
            # Snapshot: delete() mutates the registry while iterating.
            for job_id in list(self.jobs):
                job = self.jobs[job_id]
                if job.status in (JobStatus.COMPLETED, JobStatus.FAILED,
                                  JobStatus.CANCELLED, JobStatus.INTERRUPTED):
                    finished = job.finished_at
                    if finished is not None and now - finished >= self.retention_days * 86400:
                        self.delete(job_id)
            time.sleep(60)

    # ------------------------------------------------------------------ #
    # Startup recovery
    # ------------------------------------------------------------------ #
    def load_existing(self) -> None:
        """Rebuild the in-memory registry from job.json files on disk."""
        if not self.base_dir.exists():
            return
        for job_dir in sorted(self.base_dir.iterdir()):
            meta_path = job_dir / "job.json"
            if not job_dir.is_dir() or not meta_path.exists():
                continue
            try:
                data = json.loads(meta_path.read_text())
            except (OSError, json.JSONDecodeError):
                continue
            status = JobStatus(data.get("status", JobStatus.INTERRUPTED))
            if status in (JobStatus.QUEUED, JobStatus.RUNNING):
                status = JobStatus.INTERRUPTED
            job = Job(
                id=data["id"],
                dir=job_dir,
                command=data.get("command", []),
                status=status,
                exit_code=data.get("exit_code"),
                error=data.get("error") or ("server restarted while job was active"
                                             if status == JobStatus.INTERRUPTED else None),
                created_at=data.get("created_at", time.time()),
                started_at=data.get("started_at"),
                finished_at=data.get("finished_at"),
            )
            self.jobs[job.id] = job
            job.save()

    # ------------------------------------------------------------------ #
    # Job creation / execution
    # ------------------------------------------------------------------ #
    def create_job(self, input_filename: str, files: dict, extra_args: list,
                   num_procs: int = 1, user_id: Optional[str] = None) -> Job:
        job_id = uuid.uuid4().hex[:12]
        job_dir = self.base_dir / job_id
        job_dir.mkdir(parents=True)

        main_rel = _safe_relative_path(input_filename)
        for name, content in files.items():
            rel = _safe_relative_path(name)
            dest = job_dir / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(content)

        command = [self.binary_path, str(main_rel), *extra_args]
        if num_procs and num_procs > 1:
            if not self.mpi_launcher:
                raise RuntimeError(
                    "num_procs > 1 requested but no MPI launcher is configured "
                    "in this image (set PERILAB_MPI_LAUNCHER)"
                )
            command = [self.mpi_launcher, "-n", str(num_procs), self.binary_path,
                       str(main_rel), *extra_args]

        job = Job(id=job_id, dir=job_dir, command=command, user_id=user_id)
        job.save()
        with self._registry_lock:
            self.jobs[job_id] = job
        threading.Thread(target=self._run, args=(job,), daemon=True).start()
        return job

    @staticmethod
    def log_path(job: Job) -> Optional[Path]:
        """PeriLab's own log: <input_stem>_<timestamp>.log in the job dir.

        e.g. running ``impact_glass.yaml`` yields
        ``<job_dir>/impact_glass_<timestamp>.log``. Returns ``None`` if
        PeriLab hasn't created it yet (it appears partway through the run).
        """
        main_rel = job.command[1] if len(job.command) > 1 else None
        if not main_rel:
            return None
        stem = Path(main_rel).stem
        matches = sorted(job.dir.glob(f"{stem}_*.log"))
        return matches[-1] if matches else None

    def _run(self, job: Job) -> None:
        with self._sem:
            job.status = JobStatus.RUNNING
            job.started_at = time.time()
            job.save()
            # The Step: progress lines live in PeriLab's own log, which
            # appears partway through the run; the monitor tails it live.
            try:
                proc = subprocess.Popen(
                    job.command,
                    cwd=job.dir,
                )
                job._process = proc
                monitor = threading.Thread(
                    target=self._monitor_log,
                    args=(job,),
                    daemon=True,
                )
                monitor.start()
                proc.wait()
                monitor.join()
                job.exit_code = proc.returncode
                job.status = JobStatus.COMPLETED if proc.returncode == 0 else JobStatus.FAILED
            except FileNotFoundError as e:
                job.error = f"could not launch PeriLab binary: {e}"
                job.status = JobStatus.FAILED
            except Exception as e:  # noqa: BLE001 - report any launch failure to the caller
                job.error = str(e)
                job.status = JobStatus.FAILED
            finally:
                job.finished_at = time.time()
                job._process = None
                job.save()

    def _monitor_log(self, job: Job) -> None:
        """Tail PeriLab's log and update step/progress fields live.

        PeriLab writes its own log (named after the input deck + a
        timestamp); the ``Step:`` progress lines live there. It appears
        partway through the run, so we re-resolve it each loop.

        Parses lines like "[Info] Step: 1 / 2 [2.500e-04 s]".
        """
        step_re = re.compile(
            r"Step:\s*(\d+)\s*/\s*(\d+)\s*\[([0-9eE.+-]+)\s*s\]"
        )
        while True:
            log_path = self.log_path(job)
            if log_path is not None:
                try:
                    with open(log_path, "r") as f:
                        f.seek(0, io.SEEK_END)
                        f.seek(0)
                        lines = f.readlines()
                except OSError:
                    lines = []
                for line in lines:
                    m = step_re.search(line)
                    if m:
                        job.step = int(m.group(1))
                        job.total_steps = int(m.group(2))
                        try:
                            job.sim_time = float(m.group(3))
                        except ValueError:
                            pass
                        job.save()
            if job.status not in (JobStatus.QUEUED, JobStatus.RUNNING):
                break
            time.sleep(0.2)

    # ------------------------------------------------------------------ #
    # Control
    # ------------------------------------------------------------------ #
    def cancel(self, job_id: str) -> bool:
        job = self.jobs.get(job_id)
        if not job or not job._process:
            return False
        job._process.terminate()
        try:
            job._process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            job._process.kill()
        job.status = JobStatus.CANCELLED
        job.finished_at = time.time()
        job.save()
        return True

    def delete(self, job_id: str) -> None:
        job = self.jobs.get(job_id)
        if not job:
            return
        if job.status == JobStatus.RUNNING:
            self.cancel(job_id)
        shutil.rmtree(job.dir, ignore_errors=True)
        with self._registry_lock:
            self.jobs.pop(job_id, None)
