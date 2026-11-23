# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

module Dof_Sets

using Adapt: Adapt, @adapt_structure

export DofSet, create_dof_set

"""
    DofSet{Vb}

Set of degrees of freedom stored as a mask with one entry per degree of freedom
(README A6.2, A4a.2). Used for prescribed displacements and for master/condensed sets of a
model reduction. Memory is O(n_dofs), independent of the number of matrix entries.

The complement is never stored; the operators read `!mask[i]`.

# Fields
- `mask::Vb`: `true` if the degree of freedom belongs to the set
"""
struct DofSet{Vb<:AbstractVector{Bool}}
    mask::Vb
end
@adapt_structure DofSet

"""
    create_dof_set(indices, n_dofs)

Creates a host `DofSet` from a list of degree of freedom indices.
"""
function create_dof_set(indices::AbstractVector{<:Integer}, n_dofs::Int64)
    mask = falses(n_dofs)
    for i in indices
        if !(1 <= i <= n_dofs)
            throw(ArgumentError("Degree of freedom $i outside 1:$n_dofs."))
        end
        mask[i] = true
    end
    return DofSet(Vector{Bool}(mask))
end

end
