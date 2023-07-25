module Boundary_conditions
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

    for bc in bcs

        datamanager[bc["Type"]] = eval_bc(bc, nset[bc["Node Set"]], datamanager[bc["Type"]], datamanager["Coordinates"], time)

    end

end


function eval_bc(bc, nset, bcapply, coordinates, time)
    dof_mapping = Dict{String,Int8}("x" => 1, "y" => 2, "z" => 3)
    dof_mapping[bc["Coordinate"]]
    bc_value = bc["Value"]


    bcapply[nset] = evaluateString(bc["Value"])
    return bcapply
end

export init_BCs()
export apply_bc()

end