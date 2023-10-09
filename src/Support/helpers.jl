# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function dof_index(index, inds)
    dest = Vector{eltype(index)}(undef, inds * length(index))
    i = 1
    @inbounds for a in index, b in 1:inds
        dest[i+=1] = a * inds + b - inds
    end
    dest
end

function find_indices(vector, what)
    return findall(item -> item == what, vector)
end

function find_active(active)
    return [i for (i, is_active) in enumerate(active) if is_active]
end

