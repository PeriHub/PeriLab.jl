using CSV
using DataFrames

function get_bc_node_definitions(params)
    bcs = Dict{String,Any}()
    if check_element(params, "Boundary Conditions") == false
        return bcs
    end
    for entry in keys(params["Boundary Conditions"])
        bcs[entry] = params["Boundary Conditions"][entry]
    end
    return bcs
end