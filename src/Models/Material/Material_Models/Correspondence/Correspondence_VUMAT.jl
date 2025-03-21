# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_VUMAT
using StaticArrays
include("../../Material_Basis.jl")
using .Material_Basis:
                       voigt_to_matrix, matrix_to_voigt, get_Hooke_matrix, matrix_to_vector,
                       vector_to_matrix
include("../../../../Support/Geometry.jl")
include("../Zero_Energy_Control/global_control.jl")
using .Global_zero_energy_control: global_zero_energy_mode_stiffness
export fe_support
export init_model
export correspondence_name
export fields_for_local_synchronization

global vumat_file_path = ""

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
  init_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict)

Initializes the material model.

# Arguments
  - `datamanager::Data_manager`: Datamanager.
  - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.

# Returns
  - `datamanager::Data_manager`: Datamanager.
"""
function init_model(datamanager::Module,
                    nodes::Union{SubArray,Vector{Int64}},
                    material_parameter::Dict)
    # set to 1 to avoid a later check if the state variable field exists or not
    num_state_vars::Int64 = 1
    if !haskey(material_parameter, "File")
        @error "VUMAT file is not defined."
        return nothing
    end
    directory = datamanager.get_directory()
    material_parameter["File"] = joinpath(joinpath(pwd(), directory),
                                          material_parameter["File"])
    global vumat_file_path = material_parameter["File"]
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

    if !haskey(material_parameter, "VUMAT Material Name")
        @warn "No VUMAT Material Name is defined. Please check if you use it as method to check different material in your VUMAT."
        material_parameter["VUMAT Material Name"] = ""
    end
    if length(material_parameter["VUMAT Material Name"]) > 80
        @error "Due to old Fortran standards only a name length of 80 is supported"
        return nothing
    end

    if !haskey(material_parameter, "VUMAT name")
        material_parameter["VUMAT name"] = "VUMAT"
    end

    dof = datamanager.get_dof()
    ndi = dof
    nshr = 2 * dof - 3
    ntens = ndi + nshr
    datamanager.create_node_field("Internal energy", Float64, 1)
    datamanager.create_node_field("Dissipated inelastic energy", Float64, 1)
    datamanager.create_node_field("Stretch", Float64, ntens)
    datamanager.create_node_field("defGrad", Float64, ndi + 2 * nshr)
    datamanager.create_node_field("Temperature", Float64, 1)

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
    return "Correspondence VUMAT"
end

"""
    compute_stresses(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::SubArray, stress_N::SubArray, stress_NP1::SubArray, iID_jID_nID::Tuple=())

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
                          nodes::Union{SubArray,Vector{Int64}},
                          dof::Int64,
                          material_parameter::Dict,
                          time::Float64,
                          dt::Float64,
                          strain_increment::Union{SubArray,Array{Float64,3}},
                          stress_N::Union{SubArray,Array{Float64,3}},
                          stress_NP1::Union{SubArray,Array{Float64,3}})
    nnodes = length(nodes)
    # Number of normal stress components at this point
    ndi = dof
    # Number of engineering shear stress components
    nshr = 2 * dof - 3
    nstatev = material_parameter["Number of State Variables"]
    # Size of the stress or strain component array
    ntens = ndi + nshr
    nfieldv = 1
    nprops = material_parameter["Number of Properties"]
    lanneal = 1
    # only 80 charakters are supported
    cmname::Cstring = malloc_cstring(material_parameter["VUMAT Material Name"])
    coordMp = zeros(Float64, nnodes, dof)
    charLength = zeros(Float64, nnodes)
    props = datamanager.get_field("Properties")
    density = datamanager.get_field("Density")
    strainInc = zeros(Float64, nnodes, ntens)
    relSpinInc = zeros(Float64, nnodes, nshr)
    tempOld = datamanager.get_field("Temperature", "N")
    tempNew = datamanager.get_field("Temperature", "NP1")
    stretchOld = datamanager.get_field("Stretch", "N")
    stretchNew = datamanager.get_field("Stretch", "NP1")
    defGradOld = datamanager.get_field("defGrad", "N")
    defGradNew = datamanager.get_field("defGrad", "NP1")
    deformation_gradient = datamanager.get_field("Deformation Gradient")
    fieldOld = zeros(Float64, nnodes, nfieldv)
    fieldNew = zeros(Float64, nnodes, nfieldv)
    stressOld = zeros(Float64, nnodes, ntens)
    stressNew = zeros(Float64, nnodes, ntens)
    stateNew = datamanager.get_field("State Variables")
    stateOld = copy(stateNew)
    enerInternOld = datamanager.get_field("Internal energy", "N")
    enerInternNew = datamanager.get_field("Internal energy", "NP1")
    enerInelasOld = datamanager.get_field("Dissipated inelastic energy", "N")
    enerInelasNew = datamanager.get_field("Dissipated inelastic energy", "NP1")

    for iID in nodes
        strainInc[iID, :] = matrix_to_voigt(strain_increment[iID, :, :])
        stressOld[iID, :] = matrix_to_voigt(stress_N[iID, :, :])
        defGradNew[iID, :] = matrix_to_vector(deformation_gradient[iID, :, :])
        if datamanager.get_iteration() == 1
            defGradOld = defGradNew
        end
    end
    VUMAT_interface(nnodes,
                    ndi,
                    nshr,
                    nstatev,
                    nfieldv,
                    nprops,
                    lanneal,
                    0.0, #TODO: might be wrong
                    time,
                    dt,
                    cmname,
                    coordMp,
                    charLength,
                    Vector{Float64}(props[:]),
                    density,
                    strainInc,
                    relSpinInc,
                    tempOld,
                    stretchOld,
                    defGradOld,
                    fieldOld,
                    stressOld,
                    stateOld,
                    enerInternOld,
                    enerInelasOld,
                    tempNew,
                    stretchNew,
                    defGradNew,
                    fieldNew,
                    stressNew,
                    stateNew,
                    enerInternNew,
                    enerInelasNew)
    for iID in nodes
        stress_NP1[iID, :, :] = voigt_to_matrix(stressNew[iID, :])
    end
    defGradOld = defGradNew
    return stress_NP1, datamanager
end

"""
    VUMAT_interface()

VUMAT interface

# Arguments
- `nblock::Int64`: Number of blocks
- `ndir::Int64`: Number of directions
- `nshr::Int64`: Number of engineering shear stress components
- `nstatev::Int64`: Number of state variables
- `nfieldv::Int64`: Number of field variables
- `nprops::Int64`: Number of properties
- `lanneal::Bool`: Lanneal
- `stepTime::Float64`: Step time
- `totalTime::Float64`: Total time
- `dt::Float64`: Time step
- `cmname::Cstring`: Material name
- `coordMp::Matrix{Float64}`: Coordinates
- `charLength::Vector{Float64}`: Characteristic length
- `props::Vector{Float64}`: Properties
- `density::Vector{Float64}`: Density
- `strainInc::Matrix{Float64}`: Strain increment
- `relSpinInc::Matrix{Float64}`: Relative spin increment
- `tempOld::Vector{Float64}`: Old temperature
- `stretchOld::Vector{Float64}`: Old stretch
- `defGradOld::Matrix{Float64}`: Old deformation gradient
- `fieldOld::Matrix{Float64}`: Old field
- `stressOld::Matrix{Float64}`: Old stress
- `stateOld::Matrix{Float64}`: Old state
- `enerInternOld::Vector{Float64}`: Old internal energy
- `enerInelasOld::Vector{Float64}`: Old inelastic energy
- `tempNew::Vector{Float64}`: New temperature
- `stretchNew::Matrix{Float64}`: New stretch
- `defGradNew::Matrix{Float64}`: New deformation gradient
- `fieldNew::Matrix{Float64}`: New field
- `stressNew::Matrix{Float64}`: New stress
- `stateNew::Matrix{Float64}`: New state
- `enerInternNew::Vector{Float64}`: New internal energy
- `enerInelasNew::Vector{Float64}`: New inelastic energy

# Returns
- `datamanager`: Datamanager
"""
function VUMAT_interface(nblock::Int64,
                         ndir::Int64,
                         nshr::Int64,
                         nstatev::Int64,
                         nfieldv::Int64,
                         nprops::Int64,
                         lanneal::Int64,
                         stepTime::Float64,
                         totalTime::Float64,
                         dt::Float64,
                         cmname::Cstring,
                         coordMp::Matrix{Float64},
                         charLength::Vector{Float64},
                         props::Vector{Float64},
                         density::Vector{Float64},
                         strainInc::Matrix{Float64},
                         relSpinInc::Matrix{Float64},
                         tempOld::Vector{Float64},
                         stretchOld::Matrix{Float64},
                         defGradOld::Matrix{Float64},
                         fieldOld::Matrix{Float64},
                         stressOld::Matrix{Float64},
                         stateOld::Matrix{Float64},
                         enerInternOld::Vector{Float64},
                         enerInelasOld::Vector{Float64},
                         tempNew::Vector{Float64},
                         stretchNew::Matrix{Float64},
                         defGradNew::Matrix{Float64},
                         fieldNew::Matrix{Float64},
                         stressNew::Matrix{Float64},
                         stateNew::Matrix{Float64},
                         enerInternNew::Vector{Float64},
                         enerInelasNew::Vector{Float64})
    ccall((:vumat_, vumat_file_path),
          Cvoid,
          (Ref{Int64},
           Ref{Int64},
           Ref{Int64},
           Ref{Int64},
           Ref{Int64},
           Ref{Int64},
           Ref{Int64},
           Ref{Float64},
           Ref{Float64},
           Ref{Float64},
           Cstring,
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64},
           Ptr{Float64}),
          nblock,
          ndir,
          nshr,
          nstatev,
          nfieldv,
          nprops,
          lanneal,
          stepTime,
          totalTime,
          dt,
          cmname,
          coordMp,
          charLength,
          props,
          density,
          strainInc,
          relSpinInc,
          tempOld,
          stretchOld,
          defGradOld,
          fieldOld,
          stressOld,
          stateOld,
          enerInternOld,
          enerInelasOld,
          tempNew,
          stretchNew,
          defGradNew,
          fieldNew,
          stressNew,
          stateNew,
          enerInternNew,
          enerInelasNew)
end

"""
function set_subroutine_caller(sub_name::String)
Converts each letter in the string `str` to its corresponding lowercase equivalent and transform it to the symbol needed for ccall function.

# Arguments
- `sub_name::String`: The string given by the input file;

# Returns
- `Symbol`: A symbol which is the VUMAT subroutine name
"""

function set_subroutine_caller(sub_name::String)
    return eval(Meta.parse(":" * (lowercase(sub_name)) * "_"))
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
