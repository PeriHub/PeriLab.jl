# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module HETVAL
using TimerOutputs
export compute_model
export init_model
export thermal_model_name

"""
    thermal_model_name()

Gives the thermal model name. It is needed for comparison with the yaml input deck.

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
    compute_model(datamanager, nodes, thermal_parameter, time, dt)

Calculates the thermal behavior of the material. This template has to be copied, the file renamed and edited by the user to create a new flow. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `flow parameter::Dict(String, Any)`: Dictionary with flow parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
Example:
```julia
```
"""
function compute_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    thermal_parameter::Dict,
    block::Int64,
    time::Float64,
    dt::Float64,
)

    CMNAME::Cstring = malloc_cstring(thermal_parameter["HETVAL Material Name"])
    temp_N = datamanager.get_field("Temperature", "N")
    temp_NP1 = datamanager.get_field("Temperature", "NP1")
    deltaT = datamanager.get_field("Delta Temperature")
    statev = datamanager.get_field("State Variables")
    flux_N = datamanager.get_field("Heat Flow", "N")
    flux_NP1 = datamanager.get_field("Heat Flow", "NP1")
    PREDEF = datamanager.get_field("Predefined Fields")
    DPRED = datamanager.get_field("Predefined Fields Increment")

    for iID in nodes
        STATEV_temp = statev[iID, :]
        FLUX_temp = [flux_N[iID], flux_NP1[iID]]
        HETVAL_interface(
            thermal_parameter["File"],
            CMNAME,
            [temp_N[iID], deltaT[iID]],
            [time, time + dt],
            dt,
            STATEV_temp,
            FLUX_temp,
            PREDEF[iID, :],
            DPRED[iID, :],
        )
        statev[iID, :] = STATEV_temp
        flux_NP1[iID] = FLUX_temp[2]
    end

    return datamanager
end

"""
    HETVAL_interface(filename::String, CMNAME::Cstring, TEMP::Float64, TIME::Vector{Float64}, DTIME::Float64, STATEV::Vector{Float64}, FLUX::Float64, PREDEF::Vector{Float64}, DPRED::Vector{Float64})

UMAT interface

# Arguments
- `filename::String`: Filename
- `CMNAME::Cstring`: Material name
- `TEMP::Float64`: Temperature
- `TIME::Vector{Float64}`: Time
- `DTIME::Float64`: Time increment
- `STATEV::Vector{Float64}`: State variables
- `FLUX::Float64`: Heat Flow
- `PREDEF::Vector{Float64}`: Predefined
- `DPRED::Vector{Float64}`: Predefined increment

# Returns
- `datamanager`: Datamanager
"""
function HETVAL_interface(
    filename::String,
    CMNAME::Cstring,
    TEMP::Vector{Float64},
    TIME::Vector{Float64},
    DTIME::Float64,
    STATEV::Vector{Float64},
    FLUX::Vector{Float64},
    PREDEF::Vector{Float64},
    DPRED::Vector{Float64},
)
    expr = :(ccall(
        (:hetval_, $filename),
        Cvoid,
        (
            Cstring,
            Ptr{Float64},
            Ptr{Float64},
            Ref{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
        ),
        $CMNAME,
        $TEMP,
        $TIME,
        $DTIME,
        $STATEV,
        $FLUX,
        $PREDEF,
        $DPRED,
    ))
    eval(expr)
end

"""
    init_model(datamanager, nodes, thermal_parameter)

Inits the thermal model. This template has to be copied, the file renamed and edited by the user to create a new thermal. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `thermal parameter::Dict(String, Any)`: Dictionary with thermal parameter.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    thermal_parameter::Dict,
)
    # set to 1 to avoid a later check if the state variable field exists or not
    num_state_vars::Int64 = 1
    if !haskey(thermal_parameter, "File")
        @error "HETVAL file is not defined."
        return nothing
    end
    directory = datamanager.get_directory()
    thermal_parameter["File"] =
        joinpath(joinpath(pwd(), directory), thermal_parameter["File"])
    if !isfile(thermal_parameter["File"])
        @error "File $(thermal_parameter["File"]) does not exist, please check name and directory."
        return nothing
    end
    if haskey(thermal_parameter, "Number of State Variables")
        num_state_vars = thermal_parameter["Number of State Variables"]
    end
    # State variables are used to transfer additional information to the next step
    datamanager.create_constant_node_field("State Variables", Float64, num_state_vars)

    if !haskey(thermal_parameter, "HETVAL Material Name")
        @warn "No HETVAL Material Name is defined. Please check if you use it as method to check different material in your HETVAL."
        thermal_parameter["HETVAL Material Name"] = ""
    end
    if length(thermal_parameter["HETVAL Material Name"]) > 80
        @error "Due to old Fortran standards only a name length of 80 is supported"
        return nothing
    end

    if !haskey(thermal_parameter, "HETVAL name")
        thermal_parameter["HETVAL name"] = "HETVAL"
    end


    dof = datamanager.get_dof()

    if haskey(thermal_parameter, "Predefined Field Names")
        field_names = split(thermal_parameter["Predefined Field Names"], " ")
        fields = datamanager.create_constant_node_field(
            "Predefined Fields",
            Float64,
            length(field_names),
        )
        for (id, field_name) in enumerate(field_names)
            if !datamanager.has_key(String(field_name))
                @error "Predefined field ''$field_name'' is not defined in the mesh file."
                return nothing
            end
            # view or copy and than deleting the old one
            # TODO check if an existing field is a bool.
            fields[:, id] = datamanager.get_field(String(field_name))

        end
        datamanager.create_constant_node_field(
            "Predefined Fields Increment",
            Float64,
            length(field_names),
        )
    end

    return datamanager
end

"""
    fields_for_local_synchronization()

Returns a user developer defined local synchronization. This happens before each model.

The structure of the Dict must because

    synchfield = Dict(
        "Field name" =>
            Dict("upload_to_cores" => true, "dof" => datamanager.get_dof()),
    )

or

    synchfield = Dict(
        "Field name" =>
            Dict("download_from_cores" => true, "dof" => datamanager.get_dof()),
    )

# Arguments

"""
function fields_for_local_synchronization()
    return Dict()
end


"""

  function is taken from here
    https://discourse.julialang.org/t/how-to-create-a-cstring-from-a-string/98566
"""
function malloc_cstring(s::String)
    n = sizeof(s) + 1 # size in bytes + NUL terminator
    return GC.@preserve s @ccall memcpy(
        Libc.malloc(n)::Cstring,
        s::Cstring,
        n::Csize_t,
    )::Cstring
end

end
