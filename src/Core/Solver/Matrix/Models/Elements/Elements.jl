# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

module Elements

include("Element_Interface.jl")
using .Element_Interface

export AbstractElement, element_name, element_dofs, compute_element!, compute_elements!,
       n_elements, register_element, element_type, registered_elements

# The element implementations are licensed modules, included into this module at runtime
# (`Matrix_Adapter.load_elements`); they call `register_element` in `__init__`. The
# template is compiled so that it stays valid.
include("Element_template/element_template.jl")

const REGISTRY = Dict{String,Type}()

"""
    register_element(name::String, T::Type{<:AbstractElement})

Makes element type `T` available under `name` (called by the element modules in
`__init__`).
"""
register_element(name::String, T::Type{<:AbstractElement}) = (REGISTRY[name] = T; nothing)

"""
    element_type(name::String) -> Union{Type,Nothing}

The registered element type of `name`, `nothing` if no such element is loaded.
"""
element_type(name::String) = get(REGISTRY, name, nothing)

"Names of all registered elements."
registered_elements() = sort!(collect(keys(REGISTRY)))

end
