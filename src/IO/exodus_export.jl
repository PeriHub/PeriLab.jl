using Exodus

function open_result_file(filename)
    return ExodusDatabase(filename, "w")
end

function write_results(datamanager)

    if isfile("test_write.e")
        Base.Filesystem.rm("test_write.e")
    end

    # for (i, name) in enumerate(element_names)
    #     write_element_variable_name(exo, i, name)
    #     write_element_variable_values(exo, 1, 1, i, [5, 6, 7, 8])
    # end
    coords = datamanager.get_field("Coordinates")'
    block_Id = datamanager.get_field("Block_Id")
    volume = datamanager.get_field("Volume")
    # coords = [
    #     1.0 0.5 0.5 1.0 0.0 0.0 0.5 1.0 0.0
    #     1.0 1.0 0.5 0.5 1.0 0.5 0.0 0.0 0.0
    #     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    # ]

    # set the types
    maps_int_type = Int32
    ids_int_type = Int32
    bulk_int_type = Int32
    float_type = Float64

    # initialization parameters
    num_dim, num_nodes = size(coords)
    num_elems = num_nodes

    conn = reshape(collect(range(1, num_nodes)), 1, num_nodes)
    # conn = [1 2 3 4 5 6 7 8 9]

    # make some hack variables to write
    v_nodal_1 = rand(num_nodes)
    v_nodal_2 = rand(num_nodes)

    v_elem_1 = rand(num_nodes)
    v_elem_2 = rand(num_nodes)

    num_elem_blks = 1
    num_side_sets = 0
    num_node_sets = 0

    # create exodus database
    exo = ExodusDatabase(
        "test_write.e";
        maps_int_type, ids_int_type, bulk_int_type, float_type,
        num_dim, num_nodes, num_elems,
        num_elem_blks, num_node_sets, num_side_sets
    )

    @show exo

    # how to write coordinates
    write_coordinates(exo, convert(Matrix{Float64}, coords))
    # write_coordinates(exo, coords)
    # how to write a block
    write_element_block(exo, 1, "SPHERE", conn)
    # need at least one timestep to output variables
    write_time(exo, 1, 0.0)
    # write number of variables and their names
    write_number_of_nodal_variables(exo, 2)
    write_nodal_variable_names(exo, ["v_nodal_1", "v_nodal_2"])
    write_number_of_element_variables(exo, 2)
    write_element_variable_names(exo, ["v_elem_1", "v_elem_2"])
    # write variable values the 1 is for the time step
    write_nodal_variable_values(exo, 1, "v_nodal_1", convert(Vector{Float64}, v_nodal_1))
    write_nodal_variable_values(exo, 1, "v_nodal_2", convert(Vector{Float64}, v_nodal_2))
    # the first 1 is for the time step 
    # and the second 1 is for the block number
    write_element_variable_values(exo, 1, 1, "v_elem_1", v_elem_1)
    write_element_variable_values(exo, 1, 1, "v_elem_2", v_elem_2)

    # now confirm you wrote things correctly
    n_nodal_vars = read_number_of_nodal_variables(exo)
    n_elem_vars = read_number_of_element_variables(exo)
    @show n_nodal_vars
    @show n_elem_vars
    nodal_var_names = read_nodal_variable_names(exo)
    elem_var_names = read_element_variable_names(exo)
    @show nodal_var_names
    @show elem_var_names
    v_nodal_1_read = read_nodal_variable_values(exo, 1, "v_nodal_1")
    v_nodal_2_read = read_nodal_variable_values(exo, 1, "v_nodal_2")
    v_elem_1_read = read_element_variable_values(exo, 1, 1, "v_elem_1")
    v_elem_2_read = read_element_variable_values(exo, 1, 1, "v_elem_2")
    @show v_nodal_1_read ≈ v_nodal_1
    @show v_nodal_2_read ≈ v_nodal_2
    @show v_elem_1_read ≈ v_elem_1
    @show v_elem_2_read ≈ v_elem_2
    # don't forget to close the exodusdatabase, it can get corrupted otherwise if you're writing
    close(exo)
end
function close_result_file(exo)
    close(exo)
end