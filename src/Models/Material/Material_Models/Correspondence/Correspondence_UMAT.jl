# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_UMAT
using StaticArrays
include("../../Material_Basis.jl")
using .Material_Basis: voigt_to_matrix, matrix_to_voigt, get_Hooke_matrix
include("../../../../Support/Geometry.jl")
include("../Zero_Energy_Control/global_control.jl")
using .Global_zero_energy_control: global_zero_energy_mode_stiffness
export fe_support
export init_model
export correspondence_name
export fields_for_local_synchronization

global umat_file_path = ""

# export compute_model

"""
  fe_support()

Gives the information if the material supports the FEM part of PeriLab

# Arguments

# Returns
- bool: true - for FEM support; false - for no FEM support

Example:
```julia
println(fe_support())
false
```
"""
function fe_support()
    return true
end

"""
  init_model(datamanager::Module, nodes::AbstractVector{Int64}, material_parameter::Dict)

Initializes the material model.

# Arguments
  - `datamanager::Data_manager`: Datamanager.
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.

# Returns
  - `datamanager::Data_manager`: Datamanager.
"""
function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    material_parameter::Dict)
    # set to 1 to avoid a later check if the state variable field exists or not
    num_state_vars::Int64 = 1
    if !haskey(material_parameter, "File")
        @error "UMAT file is not defined."
        return nothing
    end
    directory = datamanager.get_directory()
    material_parameter["File"] = joinpath(joinpath(pwd(), directory),
                                          material_parameter["File"])
    global umat_file_path = material_parameter["File"]
    if !isfile(material_parameter["File"])
        @error "File $(material_parameter["File"]) does not exist, please check name and directory."
        return nothing
    end
    if haskey(material_parameter, "Number of State Variables")
        num_state_vars = material_parameter["Number of State Variables"]
    end
    # State variables are used to transfer additional information to the next step
    datamanager.create_constant_node_field("State Variables", Float64, num_state_vars)

    if !haskey(material_parameter, "Number of Properties")
        @error "Number of Properties must be at least equal 1"
        return nothing
    end
    # properties include the material properties, etc.
    num_props = material_parameter["Number of Properties"]
    properties = datamanager.create_constant_free_size_field("Properties", Float64,
                                                             (num_props, 1))

    for iID in 1:num_props
        if !haskey(material_parameter, "Property_$iID")
            @warn "Property_$iID is missing. Make sure that all properties are defined."
            properties[iID] = 0.0
        else
            properties[iID] = material_parameter["Property_$iID"]
        end
    end

    if !haskey(material_parameter, "UMAT Material Name")
        @warn "No UMAT Material Name is defined. Please check if you use it as method to check different material in your UMAT."
        material_parameter["UMAT Material Name"] = ""
    end
    if length(material_parameter["UMAT Material Name"]) > 80
        @error "Due to old Fortran standards only a name length of 80 is supported"
        return nothing
    end

    if !haskey(material_parameter, "UMAT name")
        material_parameter["UMAT name"] = "UMAT"
    end

    dof = datamanager.get_dof()
    sse = datamanager.create_constant_node_field("Specific Elastic Strain Energy", Float64,
                                                 1)
    spd = datamanager.create_constant_node_field("Specific Plastic Dissipation", Float64, 1)
    scd = datamanager.create_constant_node_field("Specific Creep Dissipation Energy",
                                                 Float64,
                                                 1)
    rpl = datamanager.create_constant_node_field("Volumetric heat generation per unit time",
                                                 Float64,
                                                 1)
    DDSDDT = datamanager.create_constant_node_field("Variation of the stress increments with respect to the temperature",
                                                    Float64,
                                                    3 * dof - 3)
    DRPLDE = datamanager.create_constant_node_field("Variation of RPL with respect to the strain increment",
                                                    Float64,
                                                    3 * dof - 3)
    DRPLDT = datamanager.create_constant_node_field("Variation of RPL with respect to the temperature",
                                                    Float64,
                                                    1)
    DFGRD0 = datamanager.create_constant_node_field("DFGRD0", Float64, dof,
                                                    VectorOrMatrix = "Matrix")
    # is already initialized if thermal problems are adressed
    datamanager.create_node_field("Temperature", Float64, 1)
    deltaT = datamanager.create_constant_node_field("Delta Temperature", Float64, 1)
    if haskey(material_parameter, "Predefined Field Names")
        field_names = split(material_parameter["Predefined Field Names"], " ")
        fields = datamanager.create_constant_node_field("Predefined Fields",
                                                        Float64,
                                                        length(field_names))
        for (id, field_name) in enumerate(field_names)
            if !datamanager.has_key(String(field_name))
                @error "Predefined field ''$field_name'' is not defined in the mesh file."
                return nothing
            end
            # view or copy and than deleting the old one
            # TODO check if an existing field is a bool.
            fields[:, id] = datamanager.get_field(String(field_name))
        end
        datamanager.create_constant_node_field("Predefined Fields Increment",
                                               Float64,
                                               length(field_names))
    end

    zStiff = datamanager.create_constant_node_field("Zero Energy Stiffness",
                                                    Float64,
                                                    dof, VectorOrMatrix = "Matrix")

    return datamanager
end

"""
    correspondence_name()

Gives the correspondence material name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the material.

Example:
```julia
println(correspondence_name())
"Material Template"
```
"""
function correspondence_name()
    return "Correspondence UMAT"
end

"""
    compute_stresses(datamanager::Module, nodes::AbstractVector{Int64}, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::SubArray, stress_N::SubArray, stress_NP1::SubArray, iID_jID_nID::Tuple=())

Calculates the stresses of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `iID::Int64`: Node ID.
- `dof::Int64`: Degrees of freedom
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
- `strainInc::Union{Array{Float64,3},Array{Float64,6}}`: Strain increment.
- `stress_N::SubArray`: Stress of step N.
- `stress_NP1::SubArray`: Stress of step N+1.
- `iID_jID_nID::Tuple=(): (optional) are the index and node id information. The tuple is ordered iID as index of the point,  jID the index of the bond of iID and nID the neighborID.
# Returns
- `datamanager::Data_manager`: Datamanager.
- `stress_NP1::SubArray`: updated stresses

Example:
```julia
```
"""
function compute_stresses(datamanager::Module,
                          nodes::AbstractVector{Int64},
                          dof::Int64,
                          material_parameter::Dict,
                          time::Float64,
                          dt::Float64,
                          strain_increment::Union{SubArray,Array{Float64,3}},
                          stress_N::Union{SubArray,Array{Float64,3}},
                          stress_NP1::Union{SubArray,Array{Float64,3}})
    # the notation from the Abaqus Fortran subroutine is used.
    nstatev = material_parameter["Number of State Variables"]
    nprops = material_parameter["Number of Properties"]
    props = datamanager.get_field("Properties")

    statev = datamanager.get_field("State Variables")

    stress_temp::Vector{Float64} = zeros(Float64, 3 * dof - 3)
    DDSDDE = zeros(Float64, 3 * dof - 3, 3 * dof - 3)
    SSE = datamanager.get_field("Specific Elastic Strain Energy")
    SPD = datamanager.get_field("Specific Plastic Dissipation")
    SCD = datamanager.get_field("Specific Creep Dissipation Energy")
    RPL = datamanager.get_field("Volumetric heat generation per unit time")
    DDSDDT = datamanager.get_field("Variation of the stress increments with respect to the temperature")
    DRPLDE = datamanager.get_field("Variation of RPL with respect to the strain increment")
    DRPLDT = datamanager.get_field("Variation of RPL with respect to the temperature")
    strain_N = datamanager.get_field("Strain", "N")
    temp = datamanager.get_field("Temperature", "N")
    dtemp = datamanager.get_field("Delta Temperature")
    PREDEF = datamanager.get_field_if_exists("Predefined Fields")
    DPRED = datamanager.get_field_if_exists("Predefined Fields Increment")

    if isnothing(PREDEF)
        PREDEF = zeros(size(temp))
        DPRED = zeros(size(temp))
    end

    # only 80 characters are supported
    CMNAME::Cstring = malloc_cstring(material_parameter["UMAT Material Name"])
    coords = datamanager.get_field("Coordinates")
    zStiff = datamanager.get_field("Zero Energy Stiffness")
    Kinv = datamanager.get_field("Inverse Shape Tensor")
    # Number of normal stress components at this point
    ndi = dof
    # Number of engineering shear stress components
    nshr = 2 * dof - 3
    # Size of the stress or strain component array
    ntens = ndi + nshr
    not_supported_float::Float64 = 0.0

    rot = datamanager.get_field_if_exists("Rotation Tensor")
    # rot_NP1 = datamanager.get_field_if_exists("Rotation Tensor", "NP1")

    DROT = isnothing(rot) ? zeros(length(temp), ntens, ntens) : rot #rot_NP1 - rot_N

    DFGRD0 = datamanager.get_field("DFGRD0")
    DFGRD1 = datamanager.get_field("Deformation Gradient")
    JSTEP = datamanager.get_iteration()
    KINC::Int64 = 1
    not_supported_int::Int64 = 0
    for iID in nodes
        STATEV_temp = statev[iID, :]
        SSE_temp = SSE[iID]
        SPD_temp = SPD[iID]
        SCD_temp = SCD[iID]
        RPL_temp = RPL[iID]
        DDSDDT_temp = DDSDDT[iID, :]
        DRPLDE_temp = DRPLDE[iID, :]
        DRPLDT_temp = DRPLDT[iID]
        UMAT_interface(stress_temp,
                       STATEV_temp,
                       DDSDDE,
                       SSE_temp,
                       SPD_temp,
                       SCD_temp,
                       RPL_temp,
                       DDSDDT_temp,
                       DRPLDE_temp,
                       DRPLDT_temp,
                       matrix_to_voigt(strain_N[iID, :, :]),
                       matrix_to_voigt(strain_increment[iID, :, :]),
                       [time, time + dt],
                       dt,
                       temp[iID],
                       dtemp[iID],
                       PREDEF[iID, :],
                       DPRED[iID, :],
                       CMNAME,
                       ndi,
                       nshr,
                       ntens,
                       nstatev,
                       Vector{Float64}(props[:]),
                       nprops,
                       coords[iID, :],
                       DROT[iID, :, :],
                       not_supported_float,
                       not_supported_float,
                       DFGRD0[iID, :, :],
                       DFGRD1[iID, :, :],
                       iID,
                       not_supported_int,
                       not_supported_int,
                       not_supported_int,
                       JSTEP,
                       KINC)

        statev[iID, :] = STATEV_temp
        SSE[iID] = SSE_temp
        SPD[iID] = SPD_temp
        SCD[iID] = SCD_temp
        RPL[iID] = RPL_temp
        DDSDDT[iID, :] = DDSDDT_temp
        DRPLDE[iID, :] = DRPLDE_temp
        DRPLDT[iID] = DRPLDT_temp

        #TODO: Fix this
        # Global_zero_energy_control.global_zero_energy_mode_stiffness(
        #     iID,
        #     DDSDDE,
        #     Kinv,
        #     zStiff,
        # )
        stress_NP1[iID, :, :] = voigt_to_matrix(stress_temp)

        DFGRD0 = DFGRD1
    end
end

"""
    UMAT_interface(filename::String, STRESS::Vector{Float64}, STATEV::Vector{Float64}, DDSDDE::Matrix{Float64}, SSE::Float64, SPD::Float64, SCD::Float64, RPL::Float64, DDSDDT::Vector{Float64}, DRPLDE::Vector{Float64}, DRPLDT::Float64, STRAN::Vector{Float64}, DSTRAN::Vector{Float64}, TIME::Vector{Float64}, DTIME::Float64, TEMP::Float64, DTEMP::Float64, PREDEF::Vector{Float64}, DPRED::Vector{Float64}, CMNAME::Cstring, NDI::Int64, NSHR::Int64, NTENS::Int64, NSTATEV::Int64, PROPS::Vector{Float64}, NPROPS::Int64, COORDS::Vector{Float64}, DROT::Matrix{Float64}, PNEWDT::Float64, CELENT::Float64, DFGRD0::Matrix{Float64}, DFGRD1::Matrix{Float64}, NOEL::Int64, NPT::Int64, LAYER::Int64, KSPT::Int64, JSTEP::Int64, KINC::Int64)

UMAT interface

# Arguments
- `filename::String`: Filename
- `STRESS::Vector{Float64}`: Stress
- `STATEV::Vector{Float64}`: State variables
- `DDSDDE::Matrix{Float64}`: DDSDDE
- `SSE::Float64`: SSE
- `SPD::Float64`: SPD
- `SCD::Float64`: SCD
- `RPL::Float64`: RPL
- `DDSDDT::Vector{Float64}`: DDSDDT
- `DRPLDE::Vector{Float64}`: DRPLDE
- `DRPLDT::Float64`: DRPLDT
- `STRAN::Vector{Float64}`: Strain
- `DSTRAN::Vector{Float64}`: Strain increment
- `TIME::Vector{Float64}`: Time
- `DTIME::Float64`: Time increment
- `TEMP::Float64`: Temperature
- `DTEMP::Float64`: Temperature increment
- `PREDEF::Vector{Float64}`: Predefined
- `DPRED::Vector{Float64}`: Predefined increment
- `CMNAME::Cstring`: Material name
- `NDI::Int64`: Number of normal stress components
- `NSHR::Int64`: Number of engineering shear stress components
- `NTENS::Int64`: Size of the stress or strain component array
- `NSTATEV::Int64`: Number of state variables
- `PROPS::Vector{Float64}`: Properties
- `NPROPS::Int64`: Number of properties
- `COORDS::Vector{Float64}`: Coordinates
- `DROT::Matrix{Float64}`: Rotation
- `PNEWDT::Float64`: New time step
- `CELENT::Float64`: Thickness
- `DFGRD0::Matrix{Float64}`: Deformation gradient
- `DFGRD1::Matrix{Float64}`: Deformation gradient
- `NOEL::Int64`: Element number
- `NPT::Int64`: Point number
- `LAYER::Int64`: Layer
- `KSPT::Int64`: Partition
- `JSTEP::Int64`: Step
- `KINC::Int64`: Increment

# Returns
- `datamanager`: Datamanager
"""
function UMAT_interface(STRESS::Vector{Float64},
                        STATEV::Vector{Float64},
                        DDSDDE::Matrix{Float64},
                        SSE::Float64,
                        SPD::Float64,
                        SCD::Float64,
                        RPL::Float64,
                        DDSDDT::Vector{Float64},
                        DRPLDE::Vector{Float64},
                        DRPLDT::Float64,
                        STRAN::Union{Vector{Float64},SVector{3,Float64}},
                        DSTRAN::Union{Vector{Float64},SVector{3,Float64}},
                        TIME::Vector{Float64},
                        DTIME::Float64,
                        TEMP::Float64,
                        DTEMP::Float64,
                        PREDEF::Vector{Float64},
                        DPRED::Vector{Float64},
                        CMNAME::Cstring,
                        NDI::Int64,
                        NSHR::Int64,
                        NTENS::Int64,
                        NSTATEV::Int64,
                        PROPS::Vector{Float64},
                        NPROPS::Int64,
                        COORDS::Vector{Float64},
                        DROT::Matrix{Float64},
                        PNEWDT::Float64,
                        CELENT::Float64,
                        DFGRD0::Matrix{Float64},
                        DFGRD1::Matrix{Float64},
                        NOEL::Int64,
                        NPT::Int64,
                        LAYER::Int64,
                        KSPT::Int64,
                        JSTEP::Int64,
                        KINC::Int64)
    ccall((:umat_, umat_file_path),
          Cvoid,
          (Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ref{Float64},
           Ref{Float64},
           Ref{Float64},
           Ref{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ref{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ref{Float64},
           Ref{Float64},
           Ref{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Cstring,
           Ref{Int64},
           Ref{Int64},
           Ref{Int64},
           Ref{Int64},
           Ptr{Float64},
           Ref{Int64},
           Ptr{Float64},
           Ptr{Float64},
           Ref{Float64},
           Ref{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ref{Int64},
           Ref{Int64},
           Ref{Int64},
           Ref{Int64},
           Ref{Int64},
           Ref{Int64}),
          STRESS,
          STATEV,
          DDSDDE,
          SSE,
          SPD,
          SCD,
          RPL,
          DDSDDT,
          DRPLDE,
          DRPLDT,
          STRAN,
          DSTRAN,
          TIME,
          DTIME,
          TEMP,
          DTEMP,
          PREDEF,
          DPRED,
          CMNAME,
          NDI,
          NSHR,
          NTENS,
          NSTATEV,
          PROPS,
          NPROPS,
          COORDS,
          DROT,
          PNEWDT,
          CELENT,
          DFGRD0,
          DFGRD1,
          NOEL,
          NPT,
          LAYER,
          KSPT,
          JSTEP,
          KINC)
end

function compute_stresses_ba(datamanager::Module,
                             nodes,
                             nlist,
                             dof::Int64,
                             material_parameter::Dict,
                             time::Float64,
                             dt::Float64,
                             strain_increment::Union{SubArray,Array{Float64,3},
                                                     Vector{Float64}},
                             stress_N::Union{SubArray,Array{Float64,3},Vector{Float64}},
                             stress_NP1::Union{SubArray,Array{Float64,3},
                                               Vector{Float64}})
    @error "$(correspondence_name()) not yet implemented for bond associated."
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

"""
    fields_for_local_synchronization(datamanager::Module, model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(datamanager::Module, model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #datamanager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
    return datamanager
end

end
