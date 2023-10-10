# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Material
include("../../Core/Module_inclusion/set_Modules.jl")
include("material_basis.jl")
using .Set_modules

global module_list = Set_modules.find_module_files(@__DIR__, "material_name")
Set_modules.include_files(module_list)

export compute_forces

function compute_forces(datamanager::Module, nodes::SubArray, model_param::Dict, time::Float32, dt::Float32)
    specifics = Dict{String,String}("Call Function" => "compute_forces", "Name" => "material_name")
    datamanager = Set_modules.create_module_specifics(model_param["Material Model"], module_list, specifics, (datamanager, nodes, model_param, time, dt))
    if isnothing(datamanager)
        @error "No material of name " * material["Material Model"] * " exists."
    end
    return datamanager
end
function determine_isotropic_parameter(prop)
    return get_all_elastic_moduli(prop)
end

function distribute_force_densities(datamanager, nodes)
    nlist = datamanager.get_nlist()
    bond_force = datamanager.get_field("Bond Forces")
    force_densities = datamanager.get_field("Force Densities", "NP1")
    volume = datamanager.get_field("Volume")
    force_densities[:] = distribute_forces(nodes, nlist, bond_force, volume, force_densities)
    return datamanager
end

end
