# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export set_distribution
export set_glob_to_loc
export synch_manager

"""
    set_distribution(values::Vector{Int64})

Sets the distribution globally.

# Arguments
- `values::Vector{Int64}`: The distribution.
"""
function set_distribution(values::Vector{Int64})
    data["distribution"] = values
end

"""
    set_glob_to_loc(dict::Dict{Int64,Int64})

Sets the global-to-local mapping dict globally.

# Arguments
- `dict` (array): The dict representing the global-to-local mapping.

Example:
```julia
set_glob_to_loc([1, 3, 5])  # sets the global-to-local mapping dict
```
"""
function set_glob_to_loc(dict::Dict{Int64,Int64})
    data["glob_to_loc"] = dict
end

"""
    synch_manager(synchronise_field, direction::String)

Synchronises the fields.

# Arguments
- `synchronise_field`: The function to synchronise the field.
- `direction::String`: The direction of the synchronisation.
"""
function synch_manager(synchronise_field, direction::String)
    synch_fields = get_synch_fields()
    # @debug synch_fields
    for synch_field in keys(synch_fields)
        synchronise_field(get_comm(),
                          synch_fields,
                          get_overlap_map(),
                          get_field,
                          synch_field,
                          direction)
    end
end
