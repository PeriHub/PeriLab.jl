include("../MPI_communication/MPI_send.jl")
include("../MPI_communication/MPI_receive.jl")

function field_exists(name::String)

end

function create_field(name::String, type, len::Int)
    if field_exists(name)
        return id = get_id(name)
    else
        id = set_id(name)
    end
end

function data_synchronisation()

end

function get_id(name::String)

end

function set_id(name::String)

end