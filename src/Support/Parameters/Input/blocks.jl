# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params struct BlockParams
    block_id::Int64 = req("Block ID"; min = 1)
    density::Float64 = req("Density"; min = 0, quantity = :density)
    horizon::Float64 = req("Horizon"; min = 0, quantity = :length)
    specific_heat_capacity::Union{Nothing,Float64} = opt("Specific Heat Capacity";
                                                         default = nothing, min = 0)
    material_model::Union{Nothing,String} = opt("Material Model"; default = nothing)
    damage_model::Union{Nothing,String} = opt("Damage Model"; default = nothing)
    thermal_model::Union{Nothing,String} = opt("Thermal Model"; default = nothing)
    additive_model::Union{Nothing,String} = opt("Additive Model"; default = nothing)
    pre_calculation_model::Union{Nothing,String} = opt("Pre Calculation Model";
                                                       default = nothing)
    degradation_model::Union{Nothing,String} = opt("Degradation Model"; default = nothing)
    angle_x::Union{Nothing,Float64} = opt("Angle X"; default = nothing)
    angle_y::Union{Nothing,Float64} = opt("Angle Y"; default = nothing)
    angle_z::Union{Nothing,Float64} = opt("Angle Z"; default = nothing)
    step_id::Union{Nothing,Int64,String} = opt("Step ID"; default = nothing)
    fem::Union{Nothing,Bool} = opt("FEM"; default = nothing)
end

@params struct FEMCouplingParams
    coupling_type::String = req("Coupling Type")
    pd_weight::Union{Nothing,Float64} = opt("PD Weight"; default = nothing)
    kappa::Union{Nothing,Float64} = opt("Kappa"; default = nothing)
    coupling_block::Union{Nothing,Int64} = opt("Coupling Block"; default = nothing)
end

@params struct FEMParams
    element_type::String = req("Element Type")
    degree::Union{Int64,String} = req("Degree")
    material_model::String = req("Material Model")
    coupling::Union{Nothing,FEMCouplingParams} = opt("Coupling"; default = nothing)
end

@params struct SurfaceCorrectionParams
    type::String = req("Type")
    update::Bool = opt("Update"; default = false)
end
