module Boundary_conditions
export init_BCs
export apply_bc

function check_valid_bcs(bcs, datamanager)
    # check bc
    working_bcs = []
    for bc in bcs
        for dataentry in datamanager.get_all_fields()
            if occursin(dataentry, bc["Type"])
                if bc["Coordinate"] == "z" && datamanager.get_dof() > 2
                    append!(working_bs, bc)
                    working_bs[end]["Type"] = [dataentry]
                    break
                end
            end

        end
    end
    return working_bcs
end


function init_BCs()
    nset, bcs = boundary_condition(params, datamager)
    bcs = check_valid_bcs(bcs, datamanager)
    return nset, bcs
end


function boundary_condition(params, datamanager)
    global_nset = get_bc_node_sets(params)
    bcs = get_bc_node_defintions(params)
    return datamanager.glob_to_loc(global_nset), bcs
end

function apply_bc(bcs, datamanager, time)
    dof = datamanager.get_dof()
    dof_mapping = Dict{String,Int8}("x" => 1, "y" => 2, "z" => 3)
    coordinates = datamanager.get_field("Coordinates")
    for bc in bcs
        field_to_apply_bc = datamanager.get_field(bc["Type"])
        field_to_apply_bc[bc["Node Set"], dof_mapping[bc["Coordinate"]]] = eval_bc(bc["Value"], coordinates[bc["Node Set"]], time, dof)
    end
end

function clean_up(bc)
    bc = replace(bc, "*" => ".*")
    bc = replace(bc, "/" => "./")
    bc = replace(bc, "+" => ".+")
    bc = replace(bc, "-" => ".-")
    return bc
end

function eval_bc(bc::String, coordinates, time, dof)
    # reason for global
    # https://stackoverflow.com/questions/60105828/julia-local-variable-not-defined-in-expression-eval
    bc = clean_up(bc)
    bc_value = Meta.parse(bc)
    """
    Working with if-statements
      "if t>2 0 else 20 end"
      works for scalars. If you want to evaluate a vector, please use the Julia notation as input
      "ifelse.(x .> y, 10, 20)"
    """
    global x = coordinates[:, 1]
    global y = coordinates[:, 2]
    global t = time

    if dof > 2
        global z = coordinates[:, 3]
    else
        global z = zeros(typeof(x[1]), length(x))
    end
    return zeros(Float32, length(x)) .+ eval(bc_value)
end

end