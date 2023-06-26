include("MPI_receive.jl")
include("MPI_send.jl")

function init_data_field(comm, globalLen, type)
    localLen = send_single_value_from_vector(comm, 0, globalLen, type)
    local_data = zeros(type, localLen[1], 1)
    return local_data
end