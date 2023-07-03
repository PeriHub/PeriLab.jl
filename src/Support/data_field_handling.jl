# global nnodes = set_basis_length_of_vectors(comm, 0, lechunks)
# is set in the input reader
function set_basis_length_of_vectors(comm, master, globalLen)
    msg = send_single_value_from_vector(comm, master, globalLen, Int32)
    return msg[1]
end

function init_data_field(dof, type)
    local_data = zeros(type, nnodes * dof, 1)
    return local_data
end

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

function get_block_ids()
    return
end

function get_bc_ids()
    # gives you for the core the local Id of a point to apply the bc
    return
end

function get_core_ids()
    return
end