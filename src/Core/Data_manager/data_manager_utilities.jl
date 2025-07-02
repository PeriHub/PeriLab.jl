# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export has_key

"""
    has_key(field_name::String)::Bool

Control if a key exists.
"""
function has_key(field_name::String)::Bool
    return field_name in (data["field_names"]::Vector{String})
end
