<!--
SPDX-FileCopyrightText: 2026 Christian Willberg

SPDX-License-Identifier: BSD-3-Clause
-->

# PeriLab benchmark suite

Solvers and models of PeriLab on the refined DCB (`create_dcb_mesh.jl`),
several MPI rank counts and, where available, the distributed matrix backends. Run it
once per checkout and compare:

```sh
# a second checkout, e.g. main
git worktree add --detach ../perilab_main main
julia --project=../perilab_main -e 'using Pkg; Pkg.instantiate()'

julia --project run_suite.jl project=../perilab_main label=main n_points=40000 ranks=1,8,28,56
julia --project run_suite.jl project=. label=branch n_points=40000 ranks=1,8,28,56 \
      gpu_project=<env with PeriLab + CUDA>
julia compare.jl results_main_<n>.csv results_branch_<n>.csv
```

| case | solver | models | ranks |
|---|---|---|---|
| `verlet_pd_solid` | Verlet | PD Solid Elastic | all |
| `verlet_correspondence` | Verlet | Correspondence Elastic | all |
| `verlet_pd_solid_damage` | Verlet | PD Solid Elastic + Critical Stretch | all |
| `verlet_pd_solid_thermal` | Verlet | PD Solid Elastic + Thermal Flow | all |
| `matrix_verlet_sparsearrays` | Verlet Matrix Based, SparseArrays | Correspondence | 1 |
| `linear_static_lu` | Linear Static Matrix Based, LU | Correspondence | 1 |
| `newmark_lu` | Newmark, LU | Correspondence | 1 |
| `matrix_verlet_cpu` | Verlet Matrix Based, Matrix Backend CPU | Correspondence | all |
| `matrix_verlet_cpu_damage` | same | + Critical Energy | all |
| `linear_static_bicgstab` | Linear Static, CPU + BiCGStab | Correspondence | all |
| `newmark_bicgstab` | Newmark, CPU + BiCGStab | Correspondence | all |
| `matrix_verlet_gpu` | Verlet Matrix Based, GPU | Correspondence | 1 |
| `matrix_verlet_mixed` | Verlet Matrix Based, Mixed | Correspondence | > 1 |

The last six cases need the licensed matrix elements and are skipped on checkouts without
the `Matrix Backend` option; GPU and Mixed need `gpu_project`. Every run is done twice in
one process; `t_init` (solver setup) and `t_step` (mean of the solver's per-step timer)
come from the second run. `u_norm` (norm of the displacements at the end) checks that
both checkouts and all rank counts compute the same; for `linear_static_*` the DCB has no
x constraint, so the rigid x translation (and with it `u_norm`) depends on round-off.

## Results: v2.2.5 vs v3.0.0 (main, d028818b5) vs feature/pd-matrix-integration, 40 018 points

2× Xeon Platinum 8280 (56 cores, 112 hardware threads), 1 GPU, `-bind-to core`;
t_step = seconds per time step (mean of the second run), t_init = `Solver_Manager.init`.
**The node was shared**: 8 CPU-bound processes of another user ran during the whole
benchmark, so some ranks shared a physical core and the 28/56-rank rows carry noise of
about ±30 % (see the identical explicit Verlet code paths).
v3.0.0 is the unreleased `main` (version 3.0.0 in `Project.toml`, no tag yet); v2.2.5 is
the latest release (945e3a44e), run with the same decks (v2 per-step timers mapped in
`step_timer`). Speedups are relative to v2.2.5; `u_norm` "=" means all checkouts agree to
1e-6.

```sh
julia compare.jl results_v2.2.5_40018.csv results_main_40018.csv results_branch_40018.csv
```

| case | ranks | t_step v2.2.5 [s] | t_step v3.0.0 [s] | t_step branch [s] | speedup v3.0.0 | speedup branch | u_norm |
|---|---:|---:|---:|---:|---:|---:|---|
| verlet_pd_solid | 1 | 0.14282 | 0.03386 | 0.03396 | 4.22 | 4.21 | = |
| verlet_pd_solid | 8 | 0.02791 | 0.00634 | 0.00648 | 4.41 | 4.31 | = |
| verlet_pd_solid | 28 | 0.02977 | 0.01174 | 0.01268 | 2.54 | 2.35 | = |
| verlet_pd_solid | 56 | 0.03490 | 0.02226 | 0.01658 | 1.57 | 2.11 | = |
| verlet_correspondence | 1 | 0.18563 | 0.04244 | 0.04324 | 4.37 | 4.29 | 1.765e+02 / 3.861e-03 / 3.861e-03 |
| verlet_correspondence | 8 | 0.08102 | 0.00707 | 0.00931 | 11.46 | 8.70 | 1.765e+02 / 3.861e-03 / 3.861e-03 |
| verlet_correspondence | 28 | 0.03452 | 0.00803 | 0.01131 | 4.30 | 3.05 | 1.765e+02 / 3.861e-03 / 3.861e-03 |
| verlet_correspondence | 56 | 0.03660 | 0.01361 | 0.01319 | 2.69 | 2.77 | 1.765e+02 / 3.861e-03 / 3.861e-03 |
| verlet_pd_solid_damage | 1 | 0.17222 | 0.04608 | 0.04538 | 3.74 | 3.79 | = |
| verlet_pd_solid_damage | 8 | 0.03258 | 0.00724 | 0.00811 | 4.50 | 4.02 | = |
| verlet_pd_solid_damage | 28 | 0.03409 | 0.00838 | 0.00970 | 4.07 | 3.51 | = |
| verlet_pd_solid_damage | 56 | 0.03498 | 0.01456 | 0.01642 | 2.40 | 2.13 | = |
| verlet_pd_solid_thermal | 1 | 0.14827 | 0.03642 | 0.03719 | 4.07 | 3.99 | 3.802e-03 / 3.812e-03 / 3.812e-03 |
| verlet_pd_solid_thermal | 8 | 0.04519 | 0.00712 | 0.00717 | 6.35 | 6.30 | 3.802e-03 / 3.812e-03 / 3.812e-03 |
| verlet_pd_solid_thermal | 28 | 0.03048 | 0.00868 | 0.00715 | 3.51 | 4.26 | 3.802e-03 / 3.812e-03 / 3.812e-03 |
| verlet_pd_solid_thermal | 56 | 0.04083 | 0.01056 | 0.00938 | 3.87 | 4.35 | 3.802e-03 / 3.812e-03 / 3.812e-03 |
| matrix_verlet_sparsearrays | 1 | 0.49303 | 0.03066 | 0.03062 | 16.08 | 16.10 | 1.321e+01 / 6.468e-04 / 6.468e-04 |
| linear_static_lu | 1 | 10.32842 | 9.51197 | 9.26439 | 1.09 | 1.11 | 6.294e-02 / 6.293e-02 / 6.293e-02 |
| newmark_lu | 1 | 10.00873 | 9.51033 | 9.55366 | 1.05 | 1.05 | = |
| matrix_verlet_cpu | 1 | – | – | 0.05016 | – | – | = |
| matrix_verlet_cpu | 8 | – | – | 0.00832 | – | – | = |
| matrix_verlet_cpu | 28 | – | – | 0.00690 | – | – | = |
| matrix_verlet_cpu | 56 | – | – | 0.00539 | – | – | = |
| matrix_verlet_cpu_damage | 1 | – | – | 0.09500 | – | – | = |
| matrix_verlet_cpu_damage | 8 | – | – | 0.01600 | – | – | = |
| matrix_verlet_cpu_damage | 28 | – | – | 0.04492 | – | – | = |
| matrix_verlet_cpu_damage | 56 | – | – | 0.07502 | – | – | = |
| linear_static_bicgstab | 1 | – | – | 22.26433 | – | – | = |
| linear_static_bicgstab | 8 | – | – | 3.37393 | – | – | = |
| linear_static_bicgstab | 28 | – | – | 3.12215 | – | – | = |
| linear_static_bicgstab | 56 | – | – | 2.34983 | – | – | = |
| newmark_bicgstab | 1 | – | – | 29.28211 | – | – | = |
| newmark_bicgstab | 8 | – | – | 3.93046 | – | – | = |
| newmark_bicgstab | 28 | – | – | 3.27635 | – | – | = |
| newmark_bicgstab | 56 | – | – | 2.85359 | – | – | = |
| matrix_verlet_gpu | 1 | – | – | 0.02425 | – | – | = |
| matrix_verlet_mixed | 8 | – | – | 0.02377 | – | – | = |
| matrix_verlet_mixed | 28 | – | – | 0.02765 | – | – | = |
| matrix_verlet_mixed | 56 | – | – | 0.02750 | – | – | = |

v2.2.5 vs v3:
- v3 is 4× faster per step on one rank in every explicit case, 2–11× with MPI; the
  serial LU solvers are 5–11 % faster.
- Not the same physics on these decks in v2.2.5: the correspondence Verlet cases
  (`verlet_correspondence`, `matrix_verlet_sparsearrays`) blow up (u_norm 177 and 13 instead
  of 4e-3 and 6e-4, same time step and step count as v3), the thermal case gives the
  displacements of the purely mechanical one (no thermal effect on u), and linear static
  differs by 1.4e-4. Their timings are listed but are not a like-for-like comparison.

main vs this branch:
- The existing solvers give the same displacements on both checkouts (u_norm rel. diff 0)
  and the same speed within noise.
- `Matrix Backend: CPU` scales the matrix based Verlet down to 0.0054 s/step (56 ranks),
  5.7× faster than the serial SparseArrays path (0.0307 s/step, identical u_norm).
  Its setup is slower on one rank (26 s vs 10 s) and drops with the rank count.
- BiCGStab is slower than LU on one rank, but on 56 ranks it is 4× (linear static) and
  3.3× (Newmark) faster than the serial LU solve. Linear static differs from LU by
  7e-4 in u_norm (rigid x mode of the deck is unconstrained); Newmark matches to 5e-12.
- GPU (1 rank) runs at 0.024 s/step; Mixed is limited by the GPU rank at this size
  (`GPU Weight` 16 is too high for 40k points).
- Open: `matrix_verlet_cpu_damage` gets slower beyond 8 ranks (0.016 → 0.075 s/step).
  Profiling at 40 ranks: every rank waits 40–60 ms per step at the first global collective
  of the step (`refresh_active!`); the late rank changes from step to step. It is not the
  Allreduce itself (17 µs at 56 ranks), not GC and not step allocations (~0.7 MB/step,
  mostly the existing bond field synchronisation). Root cause still to be found.

## Earlier results: Verlet Matrix Based backends, single DCB runs (2026-11-06, 2× Xeon, 112 cores, 1× RTX A6000, 50 steps)

10 k points (9 964), all backends give `max |u_y| = 2.320025825890453e-05`:

| backend | ranks | t_init [s] | t_step [s] |
|---|---:|---:|---:|
| SparseArrays | 1 | 4.57 | 0.00957 |
| CPU | 1 | 6.41 | 0.00812 |
| CPU | 8 | 0.86 | 0.00175 |
| GPU | 1 | 7.07 | 0.00143 |
| Mixed (GPU Weight 8) | 8 | 3.28 | 0.00174 |

100 k points (100 048), all backends give `max |u_y| = 4.151982590944979e-06`:

| backend | ranks | t_init [s] | t_step [s] | speedup t_step |
|---|---:|---:|---:|---:|
| SparseArrays | 1 | 24.5 | 0.0861 | 1.0 |
| CPU | 1 | 71.8 | 0.0814 | 1.1 |
| CPU | 8 | 9.87 | 0.0133 | 6.5 |
| CPU | 28 | 5.69 | 0.0098 | 8.8 |
| CPU | 56 | 3.22 | 0.0069 | 12.4 |
| GPU | 1 | 68.4 | 0.0191 | 4.5 |

Notes: the symbolic phase (sorting the element values) runs on the host and dominates
`t_init` of the 1-rank CPU and the GPU runs. The CPU runs with 28/56 ranks shared the
machine with other test jobs. Mixed on 56 ranks is measured again above (the first
measurement predates the fix that makes all ranks use the same time step): at 40 k points
it reaches 0.0275 s/step and is limited by the GPU rank (`GPU Weight` 16 too high).
