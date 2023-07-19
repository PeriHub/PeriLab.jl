using CSV
using DataFrames
function get_node_sets(params)

    nsets = Dict{String,Any}()
    if check_element(params, "Boundary Conditions") && check_element(params, "Node Sets")

        for entry in params["Boundary Conditions"]["Node Sets"]
            if occursin(".txt", params["Boundary Conditions"]["Node Sets"][entry])
                nodes = CSV.read(params["Boundary Conditions"]["Node Sets"][entry], DataFrame; delim=" ", header=false)
                nsets[entry] = nodes.Column1
            else
                nsets[entry] = params["Boundary Conditions"]["Node Sets"][entry]
            end
        end
    end
    return nsets
end

function get_node_bcs(params)
    bcs = Dict{String,Any}()
    if check_element(params, "Boundary Conditions") && check_element(params, "Definitions")
        for entry in params["Boundary Conditions"]["Definitions"]
            bcs[entry] = params["Boundary Conditions"]["Definitions"][entry]
        end
    end
end