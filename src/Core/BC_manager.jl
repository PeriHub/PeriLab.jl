# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Boundary_conditions
export init_BCs
export apply_bc_dirichlet
export apply_bc_neumann
include("../Support/Parameters/parameter_handling.jl")

using .Parameter_Handling: get_bc_definitions
"""
    check_valid_bcs(bcs::Dict{String,Any}, datamanager::Module

Check if the boundary conditions are valid

# Arguments
- `bcs::Dict{String,Any}`: The boundary conditions
- `datamanager::Module`: The data manager module
# Returns
- `working_bcs::Dict{String,Any}`: The valid boundary conditions
"""
function check_valid_bcs(bcs::Dict{String,Any}, datamanager::Module)
    # check bc
    working_bcs = Dict()
    for bc in keys(bcs)
        if haskey(bcs[bc], "Coordinate")
            dof = datamanager.get_dof()
            if bcs[bc]["Coordinate"] == "z" && dof < 3
                @warn "Boundary condition $bc is not possible with $dof DOF"
                break
            end
        end
        valid = false

        if !haskey(bcs[bc], "Type") ||
           bcs[bc]["Type"] != "Dirichlet" && bcs[bc]["Type"] != "Neumann"
            bcs[bc]["Type"] = "Dirichlet"
            @warn "Missing boundary condition type for $bc. Assuming Dirichlet."
        end

        for data_entry in datamanager.get_all_field_keys()
            if (occursin(bcs[bc]["Variable"], data_entry) && occursin("NP1", data_entry)) ||
               bcs[bc]["Variable"] == data_entry
                working_bcs[bc] = bcs[bc]
                bcs[bc]["Variable"] = data_entry
                bcs[bc]["Initial"] = bcs[bc]["Type"] == "Initial" ? true : false
                valid = true
                break
            end
        end
        if !valid
            @error "Boundary condition $bc is not valid."
            return nothing
        end

    end
    return working_bcs
end

"""
    init_BCs(params::Dict, datamanager)

Initialize the boundary conditions

# Arguments
- `params::Dict`: The parameters
- `datamanager::Module`: Datamanager
# Returns
- `bcs::Dict{Any,Any}`: The boundary conditions
"""
function init_BCs(params::Dict, datamanager)
    bcs = boundary_condition(params, datamanager)
    bcs = check_valid_bcs(bcs, datamanager)
    return bcs
end

"""
    boundary_condition(params::Dict, datamanager)

Initialize the boundary condition

# Arguments
- `params::Dict`: The parameters
- `datamanager::Module`: Datamanager
# Returns
- `bcs_out::Dict{Any,Any}`: The boundary conditions
"""
function boundary_condition(params::Dict, datamanager)
    bcs_in = get_bc_definitions(params)
    bcs_out = Dict{String,Any}()
    nsets = datamanager.get_nsets()

    for bc in keys(bcs_in)
        node_set_names = split(bcs_in[bc]["Node Set"], "+")
        node_set_names = map(r -> strip(r), node_set_names)
        bcs_out[bc] = Dict{String,Any}("Node Set" => [])
        for node_set_name in node_set_names
            if haskey(nsets, node_set_name)
                append!(
                    bcs_out[bc]["Node Set"],
                    datamanager.get_local_nodes(nsets[node_set_name]),
                )
                for entry in keys(bcs_in[bc])
                    if entry != "Node Set"
                        bcs_out[bc][entry] = bcs_in[bc][entry]
                    end
                end
            else
                @error "Node Set '$node_set_name' is missing"
                return nothing
            end
        end
    end
    return bcs_out
end

"""
    apply_bc_dirichlet(bcs::Dict, datamanager::Module, time::Float64)

Apply the boundary conditions

# Arguments
- `bcs::Dict{Any,Any}`: The boundary conditions
- `datamanager::Module`: Datamanager
- `time::Float64`: The current time
# Returns
- `datamanager::Module`: Datamanager
"""
function apply_bc_dirichlet(bcs::Dict, datamanager::Module, time::Float64)
    dof = datamanager.get_dof()
    dof_mapping = Dict{String,Int8}("x" => 1, "y" => 2, "z" => 3)
    coordinates = datamanager.get_field("Coordinates")
    for name in keys(bcs)
        bc = bcs[name]
        if bc["Type"] != "Dirichlet" ||
           bc["Variable"] == "Force DensitiesNP1" ||
           bc["Variable"] == "ForcesNP1"
            continue
        end
        field_to_apply_bc = datamanager.get_field(bc["Variable"])

        if ndims(field_to_apply_bc) > 1
            if haskey(dof_mapping, bc["Coordinate"])
                field_to_apply_bc[bc["Node Set"], dof_mapping[bc["Coordinate"]]] .= eval_bc(
                    field_to_apply_bc[bc["Node Set"], dof_mapping[bc["Coordinate"]]],
                    bc["Value"],
                    coordinates[bc["Node Set"], :],
                    time,
                    dof,
                    bc["Initial"],
                    name,
                )
            else
                @error "Coordinate in boundary condition must be x,y or z."
                return nothing
            end
        else
            field_to_apply_bc[bc["Node Set"]] .= eval_bc(
                field_to_apply_bc[bc["Node Set"]],
                bc["Value"],
                coordinates[bc["Node Set"], :],
                time,
                dof,
                bc["Initial"],
                name,
            )
        end

    end
    return datamanager
end

"""
apply_bc_dirichlet_force(bcs::Dict, datamanager::Module, time::Float64)

Apply the boundary conditions

# Arguments
- `bcs::Dict{Any,Any}`: The boundary conditions
- `datamanager::Module`: Datamanager
- `time::Float64`: The current time
# Returns
- `datamanager::Module`: Datamanager
"""
function apply_bc_dirichlet_force(bcs::Dict, datamanager::Module, time::Float64)
    dof = datamanager.get_dof()
    dof_mapping = Dict{String,Int8}("x" => 1, "y" => 2, "z" => 3)
    coordinates = datamanager.get_field("Coordinates")
    for name in keys(bcs)
        bc = bcs[name]
        if bc["Type"] != "Dirichlet"
            continue
        end
        if bc["Variable"] == "ForcesNP1"
            field_to_apply_bc = datamanager.get_field("External Forces")
        elseif bc["Variable"] == "Force DensitiesNP1"
            field_to_apply_bc = datamanager.get_field("External Force Densities")
        else
            continue
        end

        if ndims(field_to_apply_bc) > 1
            if haskey(dof_mapping, bc["Coordinate"])
                field_to_apply_bc[bc["Node Set"], dof_mapping[bc["Coordinate"]]] .= eval_bc(
                    field_to_apply_bc[bc["Node Set"], dof_mapping[bc["Coordinate"]]],
                    bc["Value"],
                    coordinates[bc["Node Set"], :],
                    time,
                    dof,
                    bc["Initial"],
                    name,
                )
            else
                @error "Coordinate in boundary condition must be x,y or z."
                return nothing
            end
        else
            field_to_apply_bc[bc["Node Set"]] .= eval_bc(
                field_to_apply_bc[bc["Node Set"]],
                bc["Value"],
                coordinates[bc["Node Set"], :],
                time,
                dof,
                bc["Initial"],
                name,
            )
        end

    end
    return datamanager
end

"""
    apply_bc_neumann(bcs::Dict, datamanager::Module, time::Float64)

Apply the boundary conditions

# Arguments
- `bcs::Dict{Any,Any}`: The boundary conditions
- `datamanager::Module`: Datamanager
- `time::Float64`: The current time
# Returns
- `datamanager::Module`: Datamanager
"""
function apply_bc_neumann(bcs::Dict, datamanager::Module, time::Float64)
    # Currently not supported
    dof = datamanager.get_dof()
    dof_mapping = Dict{String,Int8}("x" => 1, "y" => 2, "z" => 3)
    coordinates = datamanager.get_field("Coordinates")
    for name in keys(bcs)
        bc = bcs[name]
        if bc["Type"] != "Neumann"
            continue
        end
        field_to_apply_bc = datamanager.get_field(bc["Variable"])

        if ndims(field_to_apply_bc) > 1
            if haskey(dof_mapping, bc["Coordinate"])
                field_to_apply_bc[bc["Node Set"], dof_mapping[bc["Coordinate"]]] .+=
                    eval_bc(
                        field_to_apply_bc[bc["Node Set"], dof_mapping[bc["Coordinate"]]],
                        bc["Value"],
                        coordinates[bc["Node Set"], :],
                        time,
                        dof,
                        bc["Initial"],
                        name,
                    )
            else
                @error "Coordinate in boundary condition must be x,y or z."
                return nothing
            end
        else
            field_to_apply_bc[bc["Node Set"]] .+= eval_bc(
                field_to_apply_bc[bc["Node Set"]],
                bc["Value"],
                coordinates[bc["Node Set"], :],
                time,
                dof,
                bc["Initial"],
                name,
            )
        end

    end
    return datamanager
end

"""
    clean_up(bc::String)

Clean up the boundary condition

# Arguments
- `bc::String`: The boundary condition
# Returns
- `bc::String`: The cleaned up boundary condition
"""
function clean_up(bc::String)

    bc = replace(bc, ".*" => "*")
    bc = replace(bc, "./" => "/")
    bc = replace(bc, ".+" => "+")
    bc = replace(bc, ".-" => "-")
    # set space before the operator to avoid integer and float problems, because the dot is connected to the number and not the operator
    bc = replace(bc, "*" => " .* ")
    bc = replace(bc, "/" => " ./ ")
    bc = replace(bc, "+" => " .+ ")
    bc = replace(bc, "-" => " .- ")
    # to guarantee the scientific number notation
    bc = replace(bc, "e .- " => "e-")
    bc = replace(bc, "e .+ " => "e+")
    bc = replace(bc, "E .- " => "e-")
    bc = replace(bc, "E .+ " => "e+")
    return bc
end

"""
    eval_bc(field_values::Union{SubArray,Vector{Float64},Vector{Int64}}, bc::Union{Float64,Float64,Int64,String}, coordinates::Matrix{Float64}, time::Float64, dof::Int64)
Working with if-statements
"if t>2 0 else 20 end"
works for scalars. If you want to evaluate a vector, please use the Julia notation as input
"ifelse.(x .> y, 10, 20)"
"""
function eval_bc(
    field_values::Union{SubArray,Vector{Float64},Vector{Int64}},
    bc::Union{Float64,Int64,String},
    coordinates::Union{Matrix{Float64},Matrix{Int64}},
    time::Float64,
    dof::Int64,
    initial::Bool,
    name::String = "BC_1",
)
    # reason for global
    # https://stackoverflow.com/questions/60105828/julia-local-variable-not-defined-in-expression-eval
    # the yaml input allows multiple types. But for further use this input has to be a string
    bc = string(bc)
    bc = clean_up(bc)
    bc_value = Meta.parse(bc)

    if length(coordinates) == 0
        # @warn "Ignoring boundary condition $name.\n No nodes found, check Input Deck and or Node Sets."
        return field_values
    end

    global x = coordinates[:, 1]
    global y = coordinates[:, 2]
    global z = zeros(eltype(x), length(x))
    global t = time

    if dof > 2
        z = coordinates[:, 3]
    end

    value = eval(bc_value)

    if isnothing(value) || (initial && t != 0.0)
        return field_values
    end

    return zeros(Float64, length(x)) .+ value
end

end
