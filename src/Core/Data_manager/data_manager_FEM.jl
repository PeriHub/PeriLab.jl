# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export set_fem
"""
    set_fem(value::Bool)

Activates and deactivates the FEM option in PeriLab

# Arguments
- `value::Bool`: The value to set FEM active (true) or not (false).

Example:
```julia
set_fem(true)  # sets the fem_option to true
```
"""
function set_fem(value::Bool)
    if value
        @info "FEM is enabled"
    end
    data["fem_option"] = value
end
