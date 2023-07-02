using Exodus

function open_result_file(filename)
    return ExodusDatabase(filename, "w")
end

function write_results(exo, data, export_names)

    #coords = read_coordinates(exo) # matrix of n_nodes x n_dim
    #block_ids = read_element_block_ids(exo)
    #blocks = read_element_blocks(exo, block_ids) # contains connectivity information
    #nset_ids = read_node_set_ids(exo)
    #nsets = read_node_sets(exo, nset_ids) # contains nodes on boundaries
    #nodal_var_names = read_nodal_variable_names(exo)
    #elem_var_names = read_element_variable_names(exo)


    write_time(exo, 1, 0.0)
    write_number_of_nodal_variables(exo, length(export_names))
    count::Int64 = 0
    for name in export_names
        write_nodal_variable_name(exo, count += 1, name)
    end



    write_element_variable_name(exo, 1, "stress_xx")
    write_element_variable_name(exo, 2, "stress_yy")
    write_element_variable_name(exo, 3, "stress_xy")

    write_nodal_variable_value(exo, 1, 1, randn(...))
    return exo
end
function close_result_file(exo)
    close(exo)
end