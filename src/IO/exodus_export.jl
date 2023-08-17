module Write_Exodus_Results
using Exodus
export init_result_file
export init_write_results_in_exodus
export write_results_in_exodus

function create_result_file(filename, num_nodes, num_dim, num_elem_blks, num_node_sets)
    if isfile(filename)
        rm(filename)
    end
    maps_int_type = Int32
    ids_int_type = Int32
    bulk_int_type = Int32
    float_type = Float64
    num_elems = num_nodes
    num_side_sets = 0
    num_node_sets = num_node_sets
    init = Initialization(
        Int32(num_dim), Int32(num_nodes), Int32(num_elems),
        Int32(num_elem_blks), Int32(num_node_sets), Int32(num_side_sets)
    )

    return ExodusDatabase(
        filename, "w", init,
        maps_int_type, ids_int_type, bulk_int_type, float_type
    )
end

function paraview_specifics(dof)
    convention = Dict("1" => "x", "2" => "y", "3" => "z")
    return convention[string(dof)]
end

function get_block_nodes(block_Id, block)
    conn = findall(x -> x == block, block_Id)
    return reshape(conn, 1, length(conn))
end

function init_results_in_exodus(exo, output, coords, block_Id, nsets)

    if (typeof(coords) == Matrix{Float32}) || (typeof(coords) == Matrix{Float128})
        coords = convert(Array{Float64}, coords)
    end

    write_coordinates(exo, coords)
    write_number_of_variables(exo, NodalVariable, count(x -> x == true, values(output)))
    id = 0
    for name in keys(output)
        if output[name]
            id += 1
            write_name(exo, NodalVariable, id, name)
            # how to write coordinates
            #https://github.com/cmhamel/Exodus.jl/issues/118
            #just call your variables displ_x, displ_y, displ_z. Paraview will load this up as a vector and calculate things like vector magnitude for you.
        end
    end
    write_time(exo, 1, 0.0)
    for block in unique(block_Id)
        conn = get_block_nodes(block_Id, block)# virtual elements
        write_block(exo, block, "SPHERE", conn)
        #write_name(exo, block, "block_" * string(block))

    end
    return exo
end
end