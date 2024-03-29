# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_UMAT
include("../material_basis.jl")
include("../../../Support/geometry.jl")
export fe_support
export init_material_model
export correspondence_name
# export compute_forces

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
  init_material_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict)

Initializes the material model.

# Arguments
  - `datamanager::Data_manager`: Datamanager.
  - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.

# Returns
  - `datamanager::Data_manager`: Datamanager.
"""
function init_material_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict{String,Any})
  # set to 1 to avoid a later check if the state variable field exists or not
  num_state_vars::Int64 = 1
  if !haskey(material_parameter, "File")
    @error "UMAT file is not defined."
    return nothing
  end

  if !isfile(material_parameter["File"])
    @error "File $(material_parameter["File"]) not there, please check name and directory."
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
  properties = datamanager.create_constant_free_size_field("Properties", Float64, (num_props, 1))

  for iID in 1:num_props
    if !haskey(material_parameter, "Property_$iID")
      @error "Property_$iID is missing. Number of properties is $num_props and properties have to be in order without a missing number."
      return nothing
    end
    properties[iID] = material_parameter["Property_$iID"]
  end

  if !haskey(material_parameter, "UMAT Material Name")
    @warn "No UMAT Material Name is defined. Please check if you use it as method to check different material in your UMAT."
    material_parameter["UMAT Material Name"] = ""
  end
  if length(material_parameter["UMAT Material Name"]) > 80
    @error "Due to old Fortran standards only a name length of 80 is supported"
    return nothing
  end
  
  if !haskey(material_parameter,"UMAT name")
    material_parameter["UMAT name"] = "UMAT" 
  end


  dof = datamanager.get_dof()
  sse = datamanager.create_constant_node_field("Specific Elastic Strain Energy", Float64, 1)
  spd = datamanager.create_constant_node_field("Specific Plastic Dissipation", Float64, 1)
  scd = datamanager.create_constant_node_field("Specific Creep Dissipation Energy", Float64, 1)
  rpl = datamanager.create_constant_node_field("Volumetric heat generation per unit time", Float64, 1)
  DDSDDT = datamanager.create_constant_node_field("Variation of the stress increments with respect to the temperature", Float64, 3 * dof - 3)
  DRPLDE = datamanager.create_constant_node_field("Variation of RPL with respect to the strain increment", Float64, 3 * dof - 3)
  DRPLDT = datamanager.create_constant_node_field("Variation of RPL with respect to the temperature", Float64, 1)
  # is already initialized if thermal problems are adressed
  temperature = datamanager.create_constant_node_field("Temperature", Float64, 1)
  deltaT = datamanager.create_constant_node_field("Delta Temperature", Float64, 1)

  if haskey(material_parameter, "Predefined Field Names")
    field_names = split(material_parameter["Predefined Field Names"], " ")
    fields = datamanager.create_constant_node_field("Predefined Fields", Float64, length(field_names))
    for (id, field_name) in enumerate(field_names)
      if !(field_name in datamanager.get_all_field_keys())
        @error "Predefined field ''$field_name'' is not defined in the mesh file."
        return nothing
      end
      # view or copy and than deleting the old one
      fields[:, id] = datamanager.get_field(String(field_name))

    end
    datamanager.create_constant_node_field("Predefined Fields Increment", Float64, length(field_names))

    rotation::Bool, angles = datamanager.rotation_data()
    rotN, rotNP1 = datamanager.create_node_field("Rotation", Float64, "Matrix", dof)
    if rotation
      angles = datamanager.get_field("Angles")
      for iID in nodes
        rotN[iID, :, :] = Geometry.rotation_tensor(angles[iID, :])
      end
    end
  end

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
    compute_stresses(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::SubArray, stress_N::SubArray, stress_NP1::SubArray)

Calculates the stresses of the material. Based on the interface explained in [WillbergC2023](@cite).

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `dof::Int64`: Degrees of freedom
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
- `strainInc::Union{Array{Float64,3},Array{Float64,6}}`: Strain increment.
- `stress_N::SubArray`: Stress of step N.
- `stress_NP1::SubArray`: Stress of step N+1.
# Returns
- `datamanager::Data_manager`: Datamanager.
- `stress_NP1::SubArray`: updated stresses
Example:
```julia
```
"""
function compute_stresses(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::SubArray, stress_N::SubArray, stress_NP1::SubArray)
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
  PREDEF = datamanager.get_field("Predefined Fields")
  DPRED = datamanager.get_field("Predefined Fields Increment")
  # only 80 charakters are supported
  CMNAME::Cstring = malloc_cstring(material_parameter["UMAT Material Name"])
  coords = datamanager.get_field("Coordinates")
  rot_N = datamanager.get_field("Rotation", "N")
  rot_NP1 = datamanager.get_field("Rotation", "NP1")

  # Number of normal stress components at this point
  ndi = dof
  # Number of engineering shear stress components
  nshr = 2 * dof - 3
  # Size of the stress or strain component array
  ntens = ndi + nshr
  not_supported_float::Float64 = 0.0
  angles = datamanager.get_field("Angles")
  DFGRD0 = datamanager.get_field("Deformation Gradient", "N")
  DFGRD1 = datamanager.get_field("Deformation Gradient", "NP1")
  not_supported_int::Int64 = 0
  for iID in nodes
    stress_NP1[iID, :, :] = UMAT_interface
    rotNP1[iID, :, :] = Geometry.rotation_tensor(angles[iID, :])
    UMAT_interface(material_parameter["file"], material_parameter["UMAT name"]  , stress_temp, statev[iID, :], DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, matrix_to_voigt(strain_N[iID, :, :]), matrix_to_voigt(strain_increment[iID, :, :]), time, dt, temperature_N[iID], temperature_increment[iID], PREDEF[iID, :], DPRED[iID, :], CMNAME, ndi, nshr, ntens, nstatev, props, nprops, coords[iID, :], rot_NP1[iID, :, :] - rot_N[iID, :, :], not_supported_float, not_supported_float, DFGRD0, DFGRD1, not_supported_int, not_supported_int, not_supported_int, not_supported_int, not_supported_int, not_supported_int)
    stress_NP1[iID, :, :] = voigt_to_matrix(stress_temp)
  end
  # CORRESPONDENCE::UMATINT(sigmaNP1LocVoigt, statev, DDSDDE, &SSE, &SPD, &SCD, &RPL,
  #   DDSDDT, DRPLDE, &DRPLDT, strainLocVoigt, depsLocVoigt, timeArray, &dtime, temp, dtemp,
  #   &PREDEF, &DPRED, matnameArray, &nnormal, &nshr, &nstresscomp, &nstatev, props,
  #   &nprops, coords, drot, &PNEWDT, &CELENT, defGradN, defGradNP1,
  #   &NOEL, &NPT, &KSLAY, &KSPT, &JSTEP, &KINC, &nname)

  return datamanager, stress_NP1
end

"""
    compute_stresses(datamanager::Module, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::SubArray, stress_N::SubArray, stress_NP1::SubArray)

Calculates the stresses of a single node. Needed for FEM. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `dof::Int64`: Degrees of freedom
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
- `strainInc::Union{Array{Float64,3},Array{Float64,6}}`: Strain increment.
- `stress_N::SubArray`: Stress of step N.
- `stress_NP1::SubArray`: Stress of step N+1.
# Returns
- `datamanager::Data_manager`: Datamanager.
- `stress_NP1::SubArray`: updated stresses
Example:
```julia
```
"""
function UMAT_interface(filename::String, function_name::String, STRESS::Vector{Float64}, STATEV::Vector{Float64}, DDSDDE::Matrix{Float64}, SSE::Float64, SPD::Float64, SCD::Float64, RPL::Float64, DDSDDT::Vector{Float64}, DRPLDE::Vector{Float64}, DRPLDT::Float64, STRAN::Vector{Float64}, DSTRAN::Vector{Float64}, TIME::Float64, DTIME::Float64, TEMP::Float64, DTEMP::Float64, PREDEF::Vector{Float64}, DPRED::Vector{Float64}, CMNAME::Cstring, NDI::Int64, NSHR::Int64, NTENS::Int64, NSTATEV::Int64, PROPS::Vector{Float64}, NPROPS::Int64, COORDS::Vector{Float64}, DROT::Matrix{Float64}, PNEWDT::Float64, CELENT::Float64, DFGRD0::Matrix{Float64}, DFGRD1::Matrix{Float64}, NOEL::Int64, NPT::Int64, LAYER::Int64, KSPT::Int64, JSTEP::Int64, KINC::Int64)
  
  
  # TODO dynamicall change folders and subroutine name
 # ccall((:umattest_, "./src/Physics/Material/UMATs/libusertest.so"), Cvoid, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
 #sub_name = set_subroutine_caller(function_name)
 #test=":umattest_"
#f=eval(:test)
 expr = :(ccall((:umattest_, $filename), Cvoid,
        (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
         Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64},
         Ref{Float64}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Cstring, Ref{Int64}, Ref{Int64},
         Ref{Int64}, Ref{Int64}, Ptr{Float64}, Ref{Int64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64},
         Ptr{Float64}, Ptr{Float64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}),
        $STRESS, $STATEV, $DDSDDE, $SSE, $SPD, $SCD, $RPL, $DDSDDT, $DRPLDE, $DRPLDT, $STRAN, $DSTRAN, $TIME, $DTIME, $TEMP, $DTEMP, $PREDEF, $DPRED, $CMNAME, $NDI,
        $NSHR, $NTENS, $NSTATEV, $PROPS, $NPROPS, $COORDS, $DROT, $PNEWDT, $CELENT, $DFGRD0, $DFGRD1, $NOEL, $NPT, $LAYER, $KSPT, $JSTEP, $KINC))
        eval(expr)
end


# TODO Subroutine not hard coded
# # sub_name::Symbol = set_subroutine_caller(function_name)    
        

"""
function set_subroutine_caller(sub_name::String)
Converts each letter in the string `str` to its corresponding lowercase equivalent and transform it to the symbol needed for ccall function.

# Arguments
- `sub_name::String`: The string given by the input file;

# Returns
- `Symbol`: A symbol which is the UMAT subroutine name
"""

function set_subroutine_caller(sub_name::String)
   return eval(Meta.parse(":"*(lowercase(sub_name))*"_"))
end



function compute_stresses(datamanager::Module, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::Vector{Float64}, stress_N::Vector{Float64}, stress_NP1::Vector{Float64})

  hookeMatrix = get_Hooke_matrix(material_parameter, material_parameter["Symmetry"], dof)

  return hookeMatrix * strain_increment + stress_N, datamanager
end
"""

  function is taken from here
    https://discourse.julialang.org/t/how-to-create-a-cstring-from-a-string/98566
"""
function malloc_cstring(s::String)
  n = sizeof(s) + 1 # size in bytes + NUL terminator
  return GC.@preserve s @ccall memcpy(Libc.malloc(n)::Cstring,
    s::Cstring, n::Csize_t)::Cstring
end

end