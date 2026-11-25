# Parallel Computing

PeriLab runs in parallel on two levels that can be combined:

- **MPI**: the points are distributed over several processes (ranks); every rank computes its
  own points and exchanges the values of the neighbouring points of other ranks (ghosts).
- **Threads**: inside one rank, some node loops run on several Julia threads
  (`julia -t N`).

How to start PeriLab with MPI and pin the ranks is described in
[Getting Started](getting_started.md).

## Solvers

| Solver | MPI | Notes |
|---|---|---|
| `Verlet` | yes | all material, damage, thermal, contact and additive models |
| `Verlet Matrix Based`, `Matrix Backend: SparseArrays` (default) | no | serial `SparseMatrixCSC`; also needed for `Model Reduction` |
| `Verlet Matrix Based`, `Matrix Backend: CPU`, `GPU` or `Mixed` | yes | stiffness matrix distributed over the ranks, see [Matrix backends](#Matrix-backends) |
| `Linear Static Matrix Based` / `Newmark`, `Linear Solver: LU` | no | direct solver, `Matrix Backend: SparseArrays` |
| `Linear Static Matrix Based` / `Newmark`, `Linear Solver: BiCGStab` | yes | iterative solver with Jacobi preconditioner, `Matrix Backend: CPU`, `GPU` or `Mixed` |
| `Static` | yes | nonlinear solver (NLsolve) on the force field of all ranks |

A matrix based solver with `Matrix Backend: SparseArrays` stops on more than one rank with
*"Correspondence matrix based with Matrix Backend: SparseArrays is not implemented for parallel
runs."* `Model Reduction` (Craig-Bampton, Guyan) requires `Matrix Backend: SparseArrays` and
therefore runs on one rank only.

The matrix backends `CPU`, `GPU` and `Mixed` need the licensed matrix elements
(`LICENSE_KEY`).

## Matrix backends

`Matrix Backend` (in `Verlet Matrix Based`, `Linear Static Matrix Based` and `Newmark`)
selects where the stiffness matrix lives:

- `SparseArrays`: one `SparseMatrixCSC` on one rank.
- `CPU`: every rank assembles the rows of its own points and of its ghosts. A matrix-vector
  product exchanges the ghost values (halo exchange) before the local product and, for
  materials that write into rows of their neighbours (correspondence), sends the partial
  sums of the ghost rows back to their owners afterwards.
- `GPU`: like `CPU`, but every rank keeps its matrix and vectors on a GPU. The ranks of a node
  are assigned to its GPUs round robin. Requires [CUDA.jl](https://cuda.juliagpu.org) in the
  environment and a functional GPU; PeriLab loads CUDA.jl itself.
- `Mixed`: the first rank of a node per GPU runs on the GPU, all other ranks on the CPU.
  Without a functional GPU every rank runs on the CPU (with a warning).

`GPU Weight` (`Mixed`) is the number of points of a GPU rank relative to a CPU rank. A good
value is the ratio of the time per step of one CPU rank to that of the GPU rank when all CPU
ranks are busy. It grows with the problem size and with the number of CPU ranks that share the
memory bandwidth: if the GPU rank waits for the CPU ranks, increase it; if the CPU ranks wait,
decrease it.

## Models

All models run with MPI. Two things to know:

- **Contact**: the positions of the contact points are collected from all ranks into one
  global vector in every contact search (`Allreduce`). This part grows with the total number
  of contact points, not with the points per rank, and limits the speedup of contact
  problems.
- **Threads**: only the following loops run on several threads (`Polyester.@batch`); all
  other models run on one thread per rank, so `julia -t N` only helps decks that use them:
  - deformed bond geometry (pre calculation)
  - `PD Solid Elastic` and the ordinary state-based kernels (weighted volume, dilatation,
    bond forces)
  - `Critical Energy` damage and the damage index
  - bond based `Thermal Flow`
  - bond filters (`Disk`, `Rectangular_Plane`) during the setup

  Correspondence materials, contact, boundary conditions and the matrix backends run on one
  thread per rank; use more MPI ranks instead.
- **BLAS threads**: `Verlet Matrix Based` with `Model Reduction` multiplies the dense block
  of the reduced matrices with OpenBLAS threads in every step. After each call the OpenBLAS
  threads keep spinning for about 0.1 s (2^28 cycles) and take the cores away from the
  threaded loops above, which then get slower with more threads. Start Julia with
  `OPENBLAS_THREAD_TIMEOUT=4` (the threads sleep after 2^4 cycles), e.g.
  `OPENBLAS_THREAD_TIMEOUT=4 OPENBLAS_NUM_THREADS=16 julia -t 16 ...`; it has to be set
  before Julia starts. `BLAS.set_num_threads(1)` around the loops does not stop the spinning.
  If the material point region is small, splitting the cores (`julia -t 8` with
  `OPENBLAS_NUM_THREADS=8` on 16 cores) can be faster still; compare the timers of a short
  run (`--verbose`).

## Scaling

Linear speedup with the number of ranks is the ideal only if the whole step is parallel and no
resource is shared. Amdahl's law bounds the speedup of a fixed problem (strong scaling) by its
serial part [AmdahlGM1967](@cite); for a problem that grows with the number of ranks (weak
scaling) see [GustafsonJL1988](@cite). [HagerG2010](@cite) treats both, together with the
saturation of the memory bandwidth and the cost of communication.

Most of the work per step (bond loops, sparse matrix-vector products, vector updates) reads
much more memory than it computes with, so it is limited by the **memory bandwidth**, not by
the number of cores. The roofline model [WilliamsS2009](@cite) describes this limit; for
sparse matrix-vector products and Krylov solvers it is shown in [GroppWD2000](@cite) and
[WilliamsS2007](@cite). A realistic upper bound of the speedup of such a part on ``n`` ranks is
therefore the ratio of the memory bandwidths ``B(n)/B(1)`` measured with the STREAM benchmark
[McCalpinJD1995](@cite), not ``n``:

- The speedup grows almost linearly as long as the memory bandwidth of the sockets is not
  used up (typically up to 6–10 ranks per socket) and flattens afterwards. On a machine with
  several sockets distribute the ranks over all sockets, e.g. `-map-by socket` or
  `-map-by numa` with MPICH (`--map-by socket` with OpenMPI), otherwise the ranks of the
  first socket share its bandwidth while the second socket idles.
- Keep enough points per rank. Every step has a fixed cost per rank: exchange of the ghost
  values (latency and bandwidth of the messages [HockneyRW1994](@cite)), global reductions
  such as the stable time step or the inner products of BiCGStab
  [GhyselsP2014](@cite), and waiting for the slowest rank. With too few points per rank this
  cost dominates and more ranks make a step slower. In the DCB benchmark
  (`test/benchmarks/suite`) the explicit `Verlet` stopped scaling below about 5000 points per
  rank. An asynchronous, task-based peridynamics implementation that hides part of this cost
  is described in [DiehlP2020](@cite).
- The iterative solver (`BiCGStab`) needs more iterations for larger problems; at large sizes
  the iteration count matters more than the number of ranks.

`test/benchmarks/suite/run_suite.jl` measures the time per step of all solvers for several
rank counts; see its `README.md` for MPI results.
