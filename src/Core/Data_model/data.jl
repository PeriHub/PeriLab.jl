include("../../Support/parameter_handling.jl")
module DataModel

data_model = Any

function init_global_data(data)

    global n_blocks = get_number_of_blocks(data)

    #global nnodes = 
end

function add_field()

end

export init_data
export add_data
export data_model



end