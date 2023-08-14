using CSV
using DataFrames

function get_bc_node_definitions(params)
    bcs = Dict{String,Any}()
    if check_element(params, "Boundary Conditions") && check_element(params, "Definitions")
        for entry in params["Boundary Conditions"]["Definitions"]
            bcs[entry] = params["Boundary Conditions"]["Definitions"][entry]
        end
    end
end