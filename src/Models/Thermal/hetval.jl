# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module HETVAL
export compute_model
export init_model
export thermal_model_name
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_thermal, ParseContext, add_error!, join_path
import ......ParameterSpec: check!, key_patterns
@params struct HETVALParams
    file::String = req("File"; description = "HETVAL library, relative to the input deck")
    number_of_state_variables::Int64 = opt("Number of State Variables"; default = 1, min = 1)
    hetval_material_name::Union{Nothing,String} = opt("HETVAL Material Name";
                                                      default = nothing)
    hetval_name::String = opt("HETVAL name"; default = "HETVAL")
    predefined_field_names::Union{Nothing,String} = opt("Predefined Field Names";
                                                        default = nothing,
                                                        description = "space separated node field names")
end
key_patterns(::Type{HETVALParams}) = [r"^Property_\d+$" => Float64]
function check!(p::HETVALParams, path::String, ctx::ParseContext)
    if p.hetval_material_name !== nothing && length(p.hetval_material_name) > 80
        add_error!(ctx, join_path(path, "HETVAL Material Name"),
                   "at most 80 characters (Fortran)")
    end
    return nothing
end
__init__() = register_thermal("HETVAL", HETVALParams)

using .....Data_Manager
using .....PeriLabExceptions: @abort

global hetval_file_path = ""
global hetval_cmname::Cstring
# set to 1 to avoid a later check if the state variable field exists or not
global num_state_vars::Int64 = 1

"""
    thermal_model_name()

Gives the thermal model name. PeriLab loads the module because it defines this function; the input deck uses the name passed to `register_*` in `__init__()`.

# Arguments

# Returns
- `name::String`: The name of the thermal flow model.

Example:
```julia
println(flow_name())
"Thermal Template"
```
"""
function thermal_model_name()
    return "HETVAL"
end

"""
    compute_model(nodes, p, thermal, block, time, dt)

Calculates the thermal behavior of the material. This template has to be copied, the file renamed and edited by the user to create a new flow. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::HETVALParams`: The model parameters.
- `thermal::WithBase`: The typed thermal model of the block.
- `block::Int64`: The current block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
Example:
```julia
```
"""
function compute_model(nodes::AbstractVector{Int64}, p::HETVALParams, thermal,
                       block::Int64, time::Float64, dt::Float64)
    global hetval_cmname
    global num_state_vars
    temp_N::NodeScalarField{Float64} = Data_Manager.get_field("Temperature", "N")
    temp_NP1::NodeScalarField{Float64} = Data_Manager.get_field("Temperature", "NP1")
    # deltaT::NodeScalarField{Float64} = Data_Manager.get_field("Delta Temperature")
    flux_N::NodeScalarField{Float64} = Data_Manager.get_field("Heat Flow", "N")
    flux_NP1::NodeScalarField{Float64} = Data_Manager.get_field("Heat Flow", "NP1")
    statev::NodeField{Float64} = Data_Manager.get_field("State Variables")
    PREDEF::NodeField{Float64} = Data_Manager.get_field("Predefined Fields")
    DPRED::NodeField{Float64} = Data_Manager.get_field("Predefined Fields Increment")

    # Pre-allocated work arrays — reused every iteration
    T_in = Vector{Float64}(undef, 2)      # Temperature + ΔT
    t_in = Vector{Float64}(undef, 2)      # [time, time+dt]
    STATEV = Vector{Float64}(undef, num_state_vars)
    flux_out = Vector{Float64}(undef, 2)      # heat flow output

    for iID in nodes
        # Temperature inputs — fill pre-allocated buffer (no alloc)
        T_in[1] = temp_N[iID]
        T_in[2] = temp_NP1[iID] - temp_N[iID]

        t_in[1] = time
        t_in[2] = time + dt

        # Reuse state buffer — no copy from statev
        STATEV .= view(statev, iID, :)   # or copy!(STATEV, ...) if HETVAL modifies in-place

        flux_out[1] = flux_N[iID]
        flux_out[2] = flux_NP1[iID]
        # STATEV_temp = statev[iID, :]
        HETVAL_interface(hetval_cmname,
                         T_in,
                         t_in,
                         dt,
                         STATEV,
                         flux_out,
                         PREDEF[iID, :],
                         DPRED[iID, :])
        # HETVAL_interface(global hetval_cmname,
        #                  [temp_N[iID], temp_NP1[iID] - temp_N[iID]],
        #                  [time, time + dt],
        #                  dt,
        #                  STATEV_temp,
        #                  [flux_N[iID], flux_NP1[iID]],
        #                  PREDEF[iID, :],
        #                  DPRED[iID, :])
        copyto!(view(statev, iID, :), STATEV)
        # statev[iID, :] = STATEV_temp
    end
end

"""
    HETVAL_interface(CMNAME::Cstring, TEMP::Float64, TIME::Vector{Float64}, DTIME::Float64, STATEV::Vector{Float64}, FLUX::Float64, PREDEF::Vector{Float64}, DPRED::Vector{Float64})

HETVAL interface

# Arguments
- `CMNAME::Cstring`: Material name
- `TEMP::Float64`: Temperature
- `TIME::Vector{Float64}`: Time
- `DTIME::Float64`: Time increment
- `STATEV::Vector{Float64}`: State variables
- `FLUX::Float64`: Heat Flow
- `PREDEF::Vector{Float64}`: Predefined
- `DPRED::Vector{Float64}`: Predefined increment
"""
function HETVAL_interface(CMNAME::Cstring,
                          TEMP::Vector{Float64},
                          TIME::Vector{Float64},
                          DTIME::Float64,
                          STATEV::Vector{Float64},
                          FLUX::Vector{Float64},
                          PREDEF::Vector{Float64},
                          DPRED::Vector{Float64})
    ccall((:hetval_, hetval_file_path),
          Cvoid,
          (Cstring,
           Ptr{Float64},
           Ptr{Float64},
           Ref{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64}),
          CMNAME,
          TEMP,
          TIME,
          DTIME,
          STATEV,
          FLUX,
          PREDEF,
          DPRED)
end

"""
    init_model(nodes, p, thermal, block)

Inits the thermal model. This template has to be copied, the file renamed and edited by the user to create a new thermal. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p`: The model parameters.
- `thermal::WithBase`: The typed thermal model of the block.

"""
function init_model(nodes::AbstractVector{Int64}, p::HETVALParams, thermal, block::Int64)
    global num_state_vars
    directory = Data_Manager.get_directory()
    file = joinpath(pwd(), directory, p.file)
    global hetval_file_path = file
    if !isfile(file)
        @abort "File $file does not exist, please check name and directory."
        return
    end
    num_state_vars = p.number_of_state_variables
    # State variables are used to transfer additional information to the next step
    if num_state_vars == 1
        Data_Manager.create_constant_node_scalar_field("State Variables", Float64)
    else
        Data_Manager.create_constant_node_vector_field("State Variables", Float64,
                                                       num_state_vars)
    end

    if p.hetval_material_name === nothing
        @warn "No HETVAL Material Name is defined. Please check if you use it as method to check different material in your HETVAL."
    end
    global hetval_cmname = malloc_cstring(something(p.hetval_material_name, ""))

    dof = Data_Manager.get_dof()

    # any field if none are defined
    field_names = p.predefined_field_names === nothing ? ["Volume"] :
                  split(p.predefined_field_names, " ")
    n_fields = length(field_names)
    if n_fields == 1
        fields = Data_Manager.create_constant_node_scalar_field("Predefined Fields",
                                                                Float64)
    else
        fields = Data_Manager.create_constant_node_vector_field("Predefined Fields",
                                                                Float64,
                                                                n_fields)
    end
    for (id, field_name) in enumerate(field_names)
        if !Data_Manager.has_key(String(field_name))
            @abort "Predefined field ''$field_name'' is not defined in the mesh file."
            return
        end
        # view or copy and than deleting the old one
        # TODO check if an existing field is a bool.
        fields[:, id] = Data_Manager.get_field(String(field_name))
    end
    if n_fields == 1
        fields = Data_Manager.create_constant_node_scalar_field("Predefined Fields Increment",
                                                                Float64)
    else
        fields = Data_Manager.create_constant_node_vector_field("Predefined Fields Increment",
                                                                Float64,
                                                                n_fields)
    end
end

"""
    fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end

"""

  function is taken from here
    https://discourse.julialang.org/t/how-to-create-a-cstring-from-a-string/98566
"""
function malloc_cstring(s::String)
    n = sizeof(s) + 1 # size in bytes + NUL terminator
    return GC.@preserve s @ccall memcpy(Libc.malloc(n)::Cstring,
                                        s::Cstring,
                                        n::Csize_t)::Cstring
end

end
