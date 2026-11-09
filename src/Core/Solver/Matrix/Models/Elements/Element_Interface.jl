# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

module Element_Interface

using KernelAbstractions: @kernel, @index, @Const
using TimerOutputs: @timeit
using ...Sparse_CSR: launch!
using ...Matrix_Assembly: AssemblyPattern

export AbstractElement, element_name, element_dofs, compute_element!, compute_elements!,
       n_elements

"""
    AbstractElement

Supertype of all elements (README A3). An element is the integration at one integration
point: a material point (bond-based, ordinary, correspondence) or a quadrature point of
the outer integration (IGA-PD). Concrete element types must be `isbits` so they can be
passed into GPU kernels; all field data comes as a separate `data` argument
(`NamedTuple` of arrays, delivered from outside, README A0.3).

An element type implements
- `element_name(::Type{<:MyElement})::String`
- `n_elements(el::MyElement, data)::Int64`
- `element_dofs(el::MyElement, data) -> (row_dofs, col_dofs)` (host, symbolic phase)
- `compute_element!(el::MyElement, ke, start, data, e)` (device, numerical phase): writes
  the dense `|row_dofs[e]| × |col_dofs[e]|` matrix column major into `ke` from `start`
  on, without allocating.
"""
abstract type AbstractElement end

function element_name end
function element_dofs end
function compute_element! end
function n_elements end

@kernel function element_kernel!(el, ke, @Const(ke_ptr), data)
    e = @index(Global)
    @inbounds compute_element!(el, ke, ke_ptr[e], data, e)
end

"""
    compute_elements!(ke, pattern, el, data)

Fills the element buffer `ke` with one work item per element (README A3.2). The
elements are independent; the buffer positions come from the pattern.
"""
function compute_elements!(ke::AbstractVector, pattern::AssemblyPattern,
                           el::AbstractElement, data)
    @timeit "compute elements" launch!(element_kernel!, ke, length(pattern.ke_ptr) - 1, el,
                                       ke, pattern.ke_ptr, data)
    return ke
end

end
