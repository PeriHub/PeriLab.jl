# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

module Element_template

using ..Element_Interface: AbstractElement
import ..Element_Interface: element_name, element_dofs, compute_element!, n_elements
# a licensed element additionally needs: import ..Elements: register_element

"""
    ElementTemplate

Template of an element. Element implementations are licensed modules
(license-server `server/modules/Matrix_*.jl`, type `MatrixElement`): copy this file there,
rename module and type, fill the four functions and register the element in `__init__`
with `register_element("My Element", MyElement)`. The module is included into
`PD_Matrix.Elements` at runtime, so the relative imports below stay the same.

Planned elements that need rows of other integration points (README A3.6):
ordinary state-based, correspondence, IGA-PD. The assembly supports them unchanged.
"""
struct ElementTemplate <: AbstractElement end

element_name(::Type{ElementTemplate}) = "Element Template"

"Number of elements (integration points) in `data`."
n_elements(::ElementTemplate, data) = 0

"""
    element_dofs(el, data)

Host, symbolic phase: returns `(row_dofs, col_dofs)` as `VectorOfVectors{Int64}`, one
entry per element.
"""
element_dofs(::ElementTemplate, data) = error("Fill element_dofs for your element.")

"""
    compute_element!(el, ke, start, data, e)

Device, numerical phase: writes the dense element matrix of element `e` column major into
`ke[start:(start + n_rows * n_cols - 1)]`. No allocation, no dynamic types.
"""
@inline compute_element!(::ElementTemplate, ke, start, data, e) = nothing

end
