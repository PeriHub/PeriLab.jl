# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
"""
HTTP wrapper around the PeriLab binary baked into this image.

Replaces workflows like:
    docker exec <container> PeriLab /app/simulations/model.yaml -v

with:
    POST /jobs                       (upload model.yaml + mesh files)
    GET  /jobs/{id}                  (poll status)
    GET  /jobs/{id}/log?tail=200     (see output)
    GET  /jobs/{id}/files            (list results)
    GET  /jobs/{id}/files/{path}     (download a result file)
"""

from __future__ import annotations

import asyncio
import os
import shutil
from pathlib import Path
from typing import List, Optional

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, PlainTextResponse, StreamingResponse

from .jobs import Job, JobManager, JobStatus

# --------------------------------------------------------------------------- #
# Configuration (override via environment variables / docker-compose)
# --------------------------------------------------------------------------- #
PERILAB_BIN = os.environ.get("PERILAB_BIN", "/app/PeriLab/bin/PeriLab")
PERILAB_VERSION_FILE = "/app/Project.toml"
JOBS_DIR = Path(os.environ.get("PERILAB_JOBS_DIR", "/app/simulations"))
MAX_CONCURRENT_JOBS = int(os.environ.get("MAX_CONCURRENT_JOBS", "2"))
RETENTION_DAYS = float(os.environ.get("PERILAB_JOB_RETENTION_DAYS", "7"))
# Set this to "mpiexecjl" (and install MPI in the image) to allow num_procs > 1.
MPI_LAUNCHER = os.environ.get("PERILAB_MPI_LAUNCHER") or (
    shutil.which("mpiexecjl") or shutil.which("mpiexec")
)

manager = JobManager(
    JOBS_DIR, PERILAB_BIN, MAX_CONCURRENT_JOBS, MPI_LAUNCHER, RETENTION_DAYS
)

app = FastAPI(
    title="PeriLab Simulation API",
    description="HTTP wrapper for running PeriLab peridynamic simulations.",
    version="1.0.0",
)


@app.on_event("startup")
def _startup() -> None:
    manager.load_existing()


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _get_or_404(job_id: str) -> Job:
    job = manager.jobs.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="job not found")
    return job


def _parse_perilab_version() -> Optional[str]:
    """Read the PeriLab version from its Project.toml."""
    try:
        for line in Path(PERILAB_VERSION_FILE).read_text().splitlines():
            line = line.strip()
            if line.startswith("version") and "=" in line:
                return line.split("=", 1)[1].strip().strip('"')
    except OSError:
        pass
    return None


def _resolve_within_job(job: Job, rel_path: str) -> Path:
    base = job.dir.resolve()
    target = (job.dir / rel_path).resolve()
    if base != target and base not in target.parents:
        raise HTTPException(status_code=400, detail="invalid path")
    return target

# --------------------------------------------------------------------------- #
# Health
# --------------------------------------------------------------------------- #
@app.get("/health", operation_id="health")
def health():
    """Report service health and the configured PeriLab runtime."""
    binary_found = Path(PERILAB_BIN).exists()
    return {
        "status": "ok" if binary_found else "degraded",
        "perilab_binary": PERILAB_BIN,
        "binary_found": binary_found,
        "max_concurrent_jobs": MAX_CONCURRENT_JOBS,
        "mpi_launcher": MPI_LAUNCHER,
        "active_jobs": sum(1 for j in manager.jobs.values() if j.status == JobStatus.RUNNING),
    }


@app.get("/version", operation_id="version")
def version():
    """Return the periLab runtime version this API serves."""
    return {"perilab_version": _parse_perilab_version()}


# --------------------------------------------------------------------------- #
# Job submission
# --------------------------------------------------------------------------- #
@app.post("/jobs", operation_id="submit_job", summary="Submit a job")
async def submit_job(
    input_file: UploadFile = File(..., description="Main YAML input deck for PeriLab"),
    extra_files: List[UploadFile] = File(
        default=[], description="Any mesh / auxiliary files referenced by the YAML deck"
    ),
    args: Optional[str] = Form(
        default="", description="Extra CLI flags for PeriLab, e.g. '-v --dryrun'"
    ),
    num_procs: int = Form(default=1, description="MPI ranks (requires MPI launcher in image)"),
    user_id: Optional[str] = Form(
        default=None, description="Owner id; jobs are tagged with it for filtering"
    ),
):
    """Submit a new PeriLab simulation job with upload files and CLI args."""
    if not input_file.filename:
        raise HTTPException(status_code=400, detail="input_file must have a filename")

    files = {input_file.filename: await input_file.read()}
    for f in extra_files:
        if f.filename:
            files[f.filename] = await f.read()

    extra_args = args.split() if args else []

    try:
        job = manager.create_job(input_file.filename, files, extra_args, num_procs, user_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e
    except RuntimeError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e

    return job.to_dict()


# --------------------------------------------------------------------------- #
# Job status / listing
# --------------------------------------------------------------------------- #
@app.get("/jobs", operation_id="list_jobs")
def list_jobs(user_id: Optional[str] = None):
    """List all jobs, optionally filtered by `user_id`."""
    jobs = sorted(manager.jobs.values(), key=lambda j: j.created_at, reverse=True)
    if user_id:
        jobs = [j for j in jobs if j.user_id == user_id]
    return [j.to_dict() for j in jobs]


@app.get("/jobs/{job_id}", operation_id="get_job")
def get_job(job_id: str):
    """Return a single job by id, or 404 if it does not exist."""
    return _get_or_404(job_id).to_dict()


@app.post("/jobs/{job_id}/cancel", operation_id="cancel_job")
def cancel_job(job_id: str):
    """Cancel a running job, or 400 if it is not currently running."""
    job = _get_or_404(job_id)
    if not manager.cancel(job_id):
        raise HTTPException(status_code=400, detail="job is not currently running")
    return job.to_dict()


@app.delete("/jobs/{job_id}", operation_id="delete_job")
def delete_job(job_id: str):
    """Delete a job and its working directory."""
    _get_or_404(job_id)
    manager.delete(job_id)
    return {"deleted": job_id}


# --------------------------------------------------------------------------- #
# Logs
# --------------------------------------------------------------------------- #
@app.get("/jobs/{job_id}/log", response_class=PlainTextResponse, operation_id="get_log")
def get_log(job_id: str, tail: Optional[int] = None):
    """Fetch the job's log text, optionally the last `tail` lines."""
    job = _get_or_404(job_id)
    log_path = manager.log_path(job)
    if not log_path:
        return ""
    text = log_path.read_text(errors="replace")
    if tail:
        text = "\n".join(text.splitlines()[-tail:])
    return text


@app.get("/jobs/{job_id}/log/stream", operation_id="stream_log")
async def stream_log(job_id: str):
    """Server-sent, append-only tail of the log until the job finishes."""
    job = _get_or_404(job_id)
    log_path = manager.log_path(job)

    async def generator():
        pos = 0
        # PeriLab's log appears partway through the run. It's append-only
        # once created, so resolve it once and only tail it.
        while True:
            if log_path is not None:
                with open(log_path, "r", errors="replace") as f:
                    f.seek(pos)
                    chunk = f.read()
                    pos = f.tell()
                if chunk:
                    yield chunk
            if job.status not in (JobStatus.QUEUED, JobStatus.RUNNING):
                break
            await asyncio.sleep(1)

    return StreamingResponse(generator(), media_type="text/plain")


# --------------------------------------------------------------------------- #
# Result files
# --------------------------------------------------------------------------- #
@app.get("/jobs/{job_id}/files", operation_id="list_files")
def list_files(job_id: str):
    job = _get_or_404(job_id)
    return sorted(str(p.relative_to(job.dir)) for p in job.dir.rglob("*") if p.is_file())


@app.get("/jobs/{job_id}/files/{file_path:path}", operation_id="download_file")
def download_file(job_id: str, file_path: str):
    job = _get_or_404(job_id)
    target = _resolve_within_job(job, file_path)
    if not target.is_file():
        raise HTTPException(status_code=404, detail="file not found")
    return FileResponse(target, filename=target.name)
