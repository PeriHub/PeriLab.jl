# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_UMAT
using StaticArrays

using ......Data_Manager
using .......ParameterSpec: @params, register_material
import .......ParameterSpec: key_patterns
using ......PeriLabExceptions: @abort
using ......Helpers: voigt_to_matrix, matrix_to_voigt, matrix_to_engineering_voigt
using .....Material_Basis: get_Hooke_matrix, get_all_elastic_moduli, hooke_matrix
export fe_support
export init_model
export correspondence_name
export fields_for_local_synchronization

@params struct CorrespondenceUMATParams
    file::String = req("File"; description = "UMAT library, relative to the input deck")
    number_of_properties::Int64 = req("Number of Properties"; min = 1,
                                      description = "number of Property_N values passed to the UMAT")
    number_of_state_variables::Union{Nothing,Int64} = opt("Number of State Variables";
                                                          default = nothing, min = 0)
    predefined_field_names::Union{Nothing,String} = opt("Predefined Field Names";
                                                        default = nothing)
    umat_material_name::Union{Nothing,String} = opt("UMAT Material Name"; default = nothing)
    umat_name::Union{Nothing,String} = opt("UMAT name"; default = nothing,
                                           description = "name of the UMAT routine, default UMAT")
end
key_patterns(::Type{CorrespondenceUMATParams}) = [r"^Property_\d+$" => Float64]
__init__() = register_material("Correspondence UMAT", CorrespondenceUMATParams)

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
  init_model(nodes::AbstractVector{Int64}, material_parameter::Dict)

Initializes the material model.

# Arguments
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
"""
function init_model(nodes::AbstractVector{Int64},
                    material_parameter::Dict)
    # set to 1 to avoid a later check if the state variable field exists or not
    num_state_vars::Int64 = 1
    if !haskey(material_parameter, "File")
        @abort "UMAT file is not defined."
        return 1
    end
    directory = Data_Manager.get_directory()
    material_parameter["File"] = joinpath(pwd(), directory, material_parameter["File"])
    global umat_file_path = material_parameter["File"]
    if !isfile(material_parameter["File"])
        @abort "File $(material_parameter["File"]) does not exist, please check name and directory."
        return 1
    end
    if haskey(material_parameter, "Number of State Variables")
        num_state_vars = material_parameter["Number of State Variables"]
    end
    # State variables are used to transfer additional information to the next step
    if num_state_vars == 1
        Data_Manager.create_constant_node_scalar_field("State Variables", Float64)
    else
        Data_Manager.create_constant_node_vector_field("State Variables", Float64,
                                                       num_state_vars)
    end

    if !haskey(material_parameter, "Number of Properties")
        @abort "Number of Properties must be at least equal 1"
        return 1
    end
    # properties include the material properties, etc.
    num_props = material_parameter["Number of Properties"]
    properties = Data_Manager.create_constant_free_size_field("Properties", Float64,
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
        @abort "Due to old Fortran standards only a name length of 80 is supported"
    end

    if !haskey(material_parameter, "UMAT name")
        material_parameter["UMAT name"] = "UMAT"
    end

    _init_umat_fields!(nodes, get(material_parameter, "Predefined Field Names", nothing))
    dof = Data_Manager.get_dof()
    DDSDDE = Data_Manager.get_field("Material Gradient")
    get_all_elastic_moduli(material_parameter)
    symmetry::String = get(material_parameter, "Symmetry", "default")

    for iID in nodes
        @views DDSDDE[iID, :,
                      :] = get_Hooke_matrix(material_parameter,
                                            symmetry,
                                            dof,
                                            iID)
    end
end

# fields of a UMAT material (shared by the dict and the typed init)
function _init_umat_fields!(nodes::AbstractVector{Int64},
                            predefined_field_names::Union{Nothing,String})
    dof = Data_Manager.get_dof()
    sse = Data_Manager.create_constant_node_scalar_field("Specific Elastic Strain Energy",
                                                         Float64)
    spd = Data_Manager.create_constant_node_scalar_field("Specific Plastic Dissipation",
                                                         Float64)
    scd = Data_Manager.create_constant_node_scalar_field("Specific Creep Dissipation Energy",
                                                         Float64)
    rpl = Data_Manager.create_constant_node_scalar_field("Volumetric heat generation per unit time",
                                                         Float64)
    DDSDDT = Data_Manager.create_constant_node_vector_field("Variation of the stress increments with respect to the temperature",
                                                            Float64,
                                                            3 * dof - 3)
    DRPLDE = Data_Manager.create_constant_node_vector_field("Variation of RPL with respect to the strain increment",
                                                            Float64,
                                                            3 * dof - 3)
    DRPLDT = Data_Manager.create_constant_node_scalar_field("Variation of RPL with respect to the temperature",
                                                            Float64)
    DFGRD0 = Data_Manager.create_constant_node_tensor_field("DFGRD0", Float64, dof)

    DDSDDE::NodeTensorField{Float64} = Data_Manager.create_constant_node_tensor_field("Material Gradient",
                                                                                      Float64,
                                                                                      Int64((dof *
                                                                                             (dof +
                                                                                              1)) /
                                                                                            2))

    # is already initialized if thermal problems are adressed
    Data_Manager.create_node_scalar_field("Temperature", Float64)
    deltaT = Data_Manager.create_constant_node_scalar_field("Delta Temperature", Float64)
    if predefined_field_names !== nothing
        field_names = split(predefined_field_names, " ")
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

    zStiff = Data_Manager.create_constant_node_tensor_field("Zero Energy Stiffness",
                                                            Float64,
                                                            dof)
end

function init_model(nodes::AbstractVector{Int64}, p::CorrespondenceUMATParams, material)
    num_state_vars::Int64 = something(p.number_of_state_variables, 1)
    file = joinpath(pwd(), Data_Manager.get_directory(), p.file)
    global umat_file_path = file
    if !isfile(file)
        @abort "File $file does not exist, please check name and directory."
        return 1
    end
    if num_state_vars == 1
        Data_Manager.create_constant_node_scalar_field("State Variables", Float64)
    else
        Data_Manager.create_constant_node_vector_field("State Variables", Float64,
                                                       num_state_vars)
    end
    properties = Data_Manager.create_constant_free_size_field("Properties", Float64,
                                                              (p.number_of_properties, 1))
    for iID in 1:p.number_of_properties
        if !haskey(material.extras, "Property_$iID")
            @warn "Property_$iID is missing. Make sure that all properties are defined."
            properties[iID] = 0.0
        else
            properties[iID] = material.extras["Property_$iID"]
        end
    end
    if p.umat_material_name === nothing
        @warn "No UMAT Material Name is defined. Please check if you use it as method to check different material in your UMAT."
    elseif length(p.umat_material_name) > 80
        @abort "Due to old Fortran standards only a name length of 80 is supported"
    end
    _init_umat_fields!(nodes, p.predefined_field_names)
    for iID in nodes
        @views Data_Manager.get_field("Material Gradient")[iID, :, :] = hooke_matrix(_state_scaled(material),
                                                                                    Data_Manager.get_dof(),
                                                                                    iID)
    end
end

# moduli scaled by the state variable named by State Factor ID (the legacy init re-ran
# get_all_elastic_moduli after creating the State Variables field)
function _state_scaled(material)
    id = material.base.state_factor_id
    (id === nothing || material.moduli === nothing) && return material
    factor = Data_Manager.get_field("State Variables")[:, id]
    m = material.moduli
    moduli = (bulk_modulus = m.bulk_modulus .* factor,
              youngs_modulus = m.youngs_modulus .* factor,
              shear_modulus = m.shear_modulus .* factor,
              poissons_ratio = m.poissons_ratio)
    Data_Manager.get_field("Bulk_Modulus") .= moduli.bulk_modulus
    Data_Manager.get_field("Young's_Modulus") .= moduli.youngs_modulus
    Data_Manager.get_field("Shear_Modulus") .= moduli.shear_modulus
    return (base = material.base, moduli = moduli, hooke_symmetry = material.hooke_symmetry)
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
    compute_stresses(nodes::AbstractVector{Int64}, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::SubArray, stress_N::SubArray, stress_NP1::SubArray, iID_jID_nID::Tuple=())

Calculates the stresses of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
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
- `stress_NP1::SubArray`: updated stresses

Example:
```julia
```
"""
function _umat_stresses!(nodes::AbstractVector{Int64},
                         dof::Int64,
                         nstatev::Int64,
                         nprops::Int64,
                         cmname::String,
                         time::Float64,
                          dt::Float64,
                          strain_increment::AbstractArray{Float64,3},
                          stress_N::AbstractArray{Float64,3},
                          stress_NP1::AbstractArray{Float64,3})
    # the notation from the Abaqus Fortran subroutine is used.
    props = Data_Manager.get_field("Properties")

    statev = Data_Manager.get_field("State Variables")

    stress_temp::Vector{Float64} = zeros(Float64, 3 * dof - 3)
    # DDSDDE = zeros(Float64, 3 * dof - 3, 3 * dof - 3)
    ## TODO use C_voigt = Data_Manager.get_field("Material Gradient")
    DDSDDE_iD = Data_Manager.get_field("Material Gradient")
    #
    ##
    SSE = Data_Manager.get_field("Specific Elastic Strain Energy")
    SPD = Data_Manager.get_field("Specific Plastic Dissipation")
    SCD = Data_Manager.get_field("Specific Creep Dissipation Energy")
    RPL = Data_Manager.get_field("Volumetric heat generation per unit time")
    DDSDDT = Data_Manager.get_field("Variation of the stress increments with respect to the temperature")
    DRPLDE = Data_Manager.get_field("Variation of RPL with respect to the strain increment")
    DRPLDT = Data_Manager.get_field("Variation of RPL with respect to the temperature")
    strain_N = Data_Manager.get_field("Strain", "N")
    temp = Data_Manager.get_field("Temperature", "N")
    dtemp = Data_Manager.get_field("Delta Temperature")
    PREDEF = Data_Manager.get_field_if_exists("Predefined Fields")
    DPRED = Data_Manager.get_field_if_exists("Predefined Fields Increment")

    if isnothing(PREDEF)
        PREDEF = zeros(size(temp))
        DPRED = zeros(size(temp))
    end

    # only 80 characters are supported
    CMNAME::Cstring = malloc_cstring(cmname)
    coords = Data_Manager.get_field("Coordinates")
    zStiff = Data_Manager.get_field("Zero Energy Stiffness")
    Kinv = Data_Manager.get_field("Inverse Shape Tensor")
    # Number of normal stress components at this point
    ndi = dof
    # Number of engineering shear stress components
    nshr = 2 * dof - 3
    # Size of the stress or strain component array
    ntens = ndi + nshr
    not_supported_float::Float64 = 0.0

    rot = Data_Manager.get_field_if_exists("Rotation Tensor")
    # rot_NP1 = Data_Manager.get_field_if_exists("Rotation Tensor", "NP1")

    DROT = isnothing(rot) ? zeros(length(temp), ntens, ntens) : rot #rot_NP1 - rot_N

    DFGRD0 = Data_Manager.get_field("DFGRD0")
    DFGRD1 = Data_Manager.get_field("Deformation Gradient")
    JSTEP = Data_Manager.get_iteration()
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
        DDSDDE = DDSDDE_iD[iID, :, :]
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
                       matrix_to_engineering_voigt(strain_N[iID, :, :]),
                       matrix_to_engineering_voigt(strain_increment[iID, :, :]),
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
        # Global_Zero_Energy_Control.global_zero_energy_mode_stiffness(
        #     iID,
        #     DDSDDE,
        #     Kinv,
        #     zStiff,
        # )
        stress_NP1[iID, :, :] = voigt_to_matrix(stress_temp)

        DFGRD0 = DFGRD1
    end
end

compute_stresses(nodes::AbstractVector{Int64}, dof::Int64, material_parameter::Dict,
                 time::Float64, dt::Float64, strain_increment::AbstractArray{Float64,3},
                 stress_N::AbstractArray{Float64,3}, stress_NP1::AbstractArray{Float64,3}) = _umat_stresses!(nodes,
                                                                                                             dof,
                                                                                                             material_parameter["Number of State Variables"],
                                                                                                             material_parameter["Number of Properties"],
                                                                                                             material_parameter["UMAT Material Name"],
                                                                                                             time,
                                                                                                             dt,
                                                                                                             strain_increment,
                                                                                                             stress_N,
                                                                                                             stress_NP1)
compute_stresses(nodes::AbstractVector{Int64}, dof::Int64, p::CorrespondenceUMATParams,
                 material, time::Float64, dt::Float64,
                 strain_increment::AbstractArray{Float64,3},
                 stress_N::AbstractArray{Float64,3}, stress_NP1::AbstractArray{Float64,3}) = _umat_stresses!(nodes,
                                                                                                             dof,
                                                                                                             something(p.number_of_state_variables,
                                                                                                                       1),
                                                                                                             p.number_of_properties,
                                                                                                             something(p.umat_material_name,
                                                                                                                       ""),
                                                                                                             time,
                                                                                                             dt,
                                                                                                             strain_increment,
                                                                                                             stress_N,
                                                                                                             stress_NP1)
compute_stresses_ba(nodes, nlist, dof::Int64, p::CorrespondenceUMATParams, material,
                    time::Float64, dt::Float64, strain_increment, stress_N, stress_NP1) = @abort "$(correspondence_name()) not yet implemented for bond associated."


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
                        STRAN::Union{Vector{Float64},SVector{3,Float64},SVector{6,Float64}},
                        DSTRAN::Union{Vector{Float64},SVector{3,Float64},
                                      SVector{6,Float64}},
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

function compute_stresses_ba(nodes,
                             nlist,
                             dof::Int64,
                             material_parameter::Dict,
                             time::Float64,
                             dt::Float64,
                             strain_increment::Union{AbstractArray{Float64,3},
                                                     Vector{Float64}},
                             stress_N::Union{SubArray,Array{Float64,3},Vector{Float64}},
                             stress_NP1::Union{AbstractArray{Float64,3},
                                               Vector{Float64}})
    @abort "$(correspondence_name()) not yet implemented for bond associated."
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
    fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end

end
