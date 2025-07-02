# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
export set_dof
export get_dof
export set_current_time
export get_current_time
export set_step
export get_step
export set_bc_free_dof
export get_bc_free_dof
export set_max_step
export set_iteration
export get_iteration

"""
    set_dof(n::Int64)

Sets the degree of freedom (dof) value globally.

# Arguments
- `n::Int64`: The value to set as the degree of freedom.

Example:
```julia
set_dof(3)  # sets the degree of freedom to 3
```
"""
function set_dof(n::Int64)
    if n > 3 || n < 2
        @error "Degree of freedom $n is not supported."
    end
    data["dof"] = n
end

"""
    get_dof()

Retrieves the degree of freedom (dof) value.

# Returns
- `dof` (integer): The current degree of freedom value.

Example:
```julia
get_dof()  # returns the current degree of freedom
```
"""
function get_dof()::Int64
    return data["dof"]
end

"""
    set_current_time(time::Float64)

Set the current time of the simulation.

# Arguments
- `time::Float64`: The current time of the simulation.
"""
function set_current_time(time::Float64)
    data["current_time"] = time
end

"""
    get_current_time()

Get the current time of the simulation.

# Returns
- `Float64`: The current time of the simulation.
"""
function get_current_time()::Float64
    return data["current_time"]::Float64
end

"""
    set_step(step::Int64)

Set the step of the simulation.

# Arguments
- `step::Int64`: The step of the simulation.

"""
function set_step(step::Union{Int64,Nothing})
    data["step"] = isnothing(step) ? -1 : step
end

"""
    get_step()

Get the step of the simulation.

# Returns
- `Int64`: The step of the simulation.

"""
function get_step()::Union{Int64,Nothing}
    step = data["step"]::Int64  # Type assertion
    return step == -1 ? nothing : step
end

"""
    set_max_step(max_step::Int64)

Set the max_step of the simulation.

# Arguments
- `max_step::Int64`: The max_step of the simulation.

"""
function set_max_step(max_step::Union{Int64,Nothing})
    data["max_step"] = isnothing(max_step) ? -1 : max_step
end

"""
    get_max_step()

Get the max_step of the simulation.

# Returns
- `Int64`: The max_step of the simulation.

"""
function get_max_step()::Union{Int64}
    return data["max_step"]
end

"""
    get_bc_free_dof()

Get all dof without displacment boundary conditions.

# Returns
- `Vector{Tuple{Int64, Int64}}`: The point and dof without boundary condition.

"""
function get_bc_free_dof()
    return data["BC_free_dof"]
end

"""
    set_bc_free_dof(values::Vector{Tuple{Int64, Int64}})

Set all dof without displacment boundary conditions.

# Returns
-

"""
function set_bc_free_dof(values::Vector{Int64})
    data["BC_free_dof"] = values
end

"""
    set_iteration(iteration::Int64)

Set the iteration of the simulation.

# Arguments
- `iteration::Int64`: The iteration of the simulation.

"""
function set_iteration(iteration::Int64)
    data["iteration"] = iteration
end

"""
    get_step()

Get the iteration of the simulation.

# Returns
- `Int64`: The iteration of the simulation.

"""
function get_iteration()
    return data["iteration"]
end
