# Bond storage (CSR layout)

This page specifies how PeriLab stores bond data, i.e. every quantity that exists once
per bond `iID -> nlist[iID][jID]`.

## Motivation

All values of a bond field live in one contiguous buffer instead of one heap object per
node and bond. This keeps memory use and allocations low, makes whole-field operations
simple loops over flat data, and is the precondition for GPU kernels (KernelAbstractions)
and CUDA-aware MPI, both of which need one contiguous buffer per field.

## Layout

All bond fields of a rank share **one** offset vector `elem_ptr` (CSR row pointer):

```
elem_ptr :: Vector{Int64}, length nnodes + 1
elem_ptr[1] = 1
elem_ptr[i + 1] = elem_ptr[i] + number_of_neighbors[i]
bonds of node i  ->  flat indices elem_ptr[i]:(elem_ptr[i + 1] - 1)
```

Each bond field is an `ArraysOfArrays.VectorOfVectors` over a flat data vector and the
shared `elem_ptr`:

| Kind | Type | Flat element |
|---|---|---|
| neighborhood list | `BondScalarState{Int64}` | `Int64` |
| scalar bond state | `BondScalarState{T}` | `T` |
| vector bond state | `BondVectorState{T,D}` | `SVector{D,T}` |
| tensor bond state | `BondTensorState{T,D,L}` | `SMatrix{D,D,T,L}`, `L = D * D` |

The aliases name the fully concrete `VectorOfArrays{...}` type, so a local annotation such
as `bond_damage::BondScalarState{Float64} = get_field(...)` stays type stable. For vector and
tensor states the dimension `D` is part of the type; code that only knows `dof` at run time
must reach the kernel through a function barrier (`Val(dof)` or simply a call that takes the
field as argument).

`field[iID]` is a non-allocating `SubArray` view of the flat data, `field[iID][jID]` is the
bond value itself. `flatview(field)` returns the flat data vector.

### Lifecycle of `elem_ptr`

- Built from `"Number of Neighbors"` on the first bond field creation and stored in the
  Data_Manager (`get_bond_ptr()`).
- Every later bond field creation checks that `elem_ptr` still matches
  `"Number of Neighbors"` and aborts otherwise. The bond structure is fixed for the rest of
  the run.
- Cleared by `initialize_data()`.
- With MPI, `"Number of Neighbors"` is the length of the **local** neighborhood list. For
  owned nodes that is the full neighborhood; responder (ghost) nodes only keep the neighbors
  that exist on the rank. Bond data of responders is exchanged for exactly these bonds.
- Checkpoints serialize all fields and the Data_Manager state in one stream, so the shared
  identity of `elem_ptr` survives a restore.

The bond structure never changes after the neighborhood list has been distributed to the
ranks. Contact uses its own data structures, additive manufacturing and damage models only
change bond *values* (bond damage).

### What stays nested

- The **global** neighborhood list on rank 0 before distribution. It is edited in place
  (bond filters, consistent FE topology) and is converted to CSR in
  `distribute_neighborhoodlist_to_cores`. Preprocessing functions that take it are typed
  `NestedNeighborhoodList` (`= Vector{Vector{Int64}}`).
- `"FilteredNeighborhoodlist"` (bond filters with contact). It holds bond indices, not bond
  values, and its per-node lengths differ from `"Number of Neighbors"`, so it cannot share
  `elem_ptr`. It is stored as a plain `Vector{Vector{Int64}}`.
- Element fields (`"Number of Element Neighbors"`) and model-private caches such as the
  `all_CB_tensors` of the matrix based correspondence model.

## Writing bond values

`SVector` and `SMatrix` are immutable. The rules:

| Field kind | Allowed | Not allowed |
|---|---|---|
| scalar | `f[i][j] = x`, `f[i][j] += x`, `f[i] .= x` | |
| vector | `f[i][j] = SVector(...)`, `f[i][j] += v`, `f[i][j] = f[i][j] * s` | `f[i][j][m] = x`, `f[i][j] .= v`, `mul!(f[i][j], ...)`, mutating `f[i][j]` through an alias or loop variable |
| tensor | `f[i][j] = σ::SMatrix` | `f[i][j, :, :] = ...`, views into a bond tensor that are then mutated |

Component writes fail loudly (`setindex!` on an `SVector` throws), they cannot go wrong
silently.

Whole-field operations use the flat data:

```julia
copyto!(flatview(N), flatview(NP1))                  # switch_bonds!
flatview(C) .= flatview(A) .- flatview(B)            # sub_in_place! etc.
fill!(view(flatview(f), ptr[i]:(ptr[i + 1] - 1)), v) # fill_field! for active nodes
```

## Kernel signatures

Kernels receive either a whole field or a node slice `field[iID]`. Signatures use abstract
types (`AbstractVector{Int64}`, `AbstractVector{Float64}`,
`AbstractVector{<:AbstractVector{Float64}}`) instead of `Vector{...}`, so that the same
kernel works for CSR fields, their `SubArray` node slices and for small hand-built nested
vectors in unit tests.

### Function barriers for vector and tensor states

`BondVectorState{Float64}` and `BondTensorState{Float64}` leave the dimension open, so a
local annotated with them (`x::BondVectorState{Float64} = get_field(...)`) is not concrete.
Bond loops therefore live in kernels that take the fields as arguments; the model function
only fetches the fields and calls the kernel.

Model functions annotate fetched vector states with `BondVectorField{T}`
(`= Union{BondVectorState{T,2},BondVectorState{T,3}}`), not with the open
`BondVectorState{T}`. The union lets the compiler split a kernel call into static calls, so
no argument is boxed. Julia splits at most 4 combinations, i.e. two vector-state arguments
per call; kernel signatures themselves keep `BondVectorState{T}`.

### Round-off

StaticArrays evaluates small matrix products with `muladd`, which may be fused (FMA), so
results can differ in the last bit from generic or BLAS routines. Examples that amplify
perturbations (e.g. `examples/DCB/DCBmodel_ba.yaml`, roughly a factor 2.3 per step) show
this as larger relative differences after a few steps. Fullscale tests therefore compare
with tolerances, not bitwise:

- `test_thermal_flow`: `Heat Flow` (about 1e11) is compared with a relative tolerance of 1e-12.
- `test_linear_static_matrix_based_solver`: the DCB decks have no constraint in x, so the x
  displacement contains a rigid body translation that depends only on round-off (strains are
  identical). These decks use their own command file without `DisplacementsX`.

### Type stability

A CSR bond field is an immutable struct. Passing it through a dynamically dispatched call
(an argument of unknown type at the call site, e.g. after `x = invert(...)` or any other
type-unstable result) boxes it: 32 bytes per field and call. Kernels called once per bond
must therefore be reached through type-stable code; `Data_Manager` stores fields pre-boxed,
so `get_field` itself does not allocate.

## MPI and output

- `create_bond_overlap_map!` adds the flat bond indices of every rank pair to the overlap
  map: `"Bond Controller"` are the bonds of owned nodes whose neighbor exists on the other
  rank, `"Bond Responder"` are the complete responder rows. Both follow the order of the
  global neighborhood list, so `synch_controller_bonds_to_responder_flattened` sends
  `flatview(field)[controller bonds]` (reinterpreted as the scalar element type) and writes
  it to `flatview(field)[responder bonds]`. The size of a bond value comes from the field,
  the registered `dof` is not used.
- Bond norms of the contact bond filter are distributed with the same rule: a local row
  keeps the bonds to neighbors on the rank.
- Exodus bond export reads `field[node][neighbor][...]`; `SVector`/`SMatrix` bond values are
  indexed like nested vectors, so the export needs no layout specific code.
