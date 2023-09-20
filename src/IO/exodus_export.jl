module Write_Exodus_Results
using Exodus
using Pkg
export get_paraviewCoordinates
export init_result_file
export init_write_results_in_exodus
export write_results_in_exodus
export merge_exodus_file

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
    init = Initialization{bulk_int_type}(
        Int32(num_dim), Int32(num_nodes), Int32(num_elems),
        Int32(num_elem_blks), Int32(num_node_sets), Int32(num_side_sets)
    )
    @info "Create output " * filename
    return ExodusDatabase(
        filename, "w", init,
        maps_int_type, ids_int_type, bulk_int_type, float_type
    )
end

function paraview_specifics(dof)
    convention = Dict(1 => "x", 2 => "y", 3 => "z")
    return convention[dof]
end

function get_paraviewCoordinates(dof, refDof)
    if dof > refDof
        @warn "Reference dof are to small and set to used dof"
        refDof = dof
    end
    if refDof < 4
        return paraview_specifics(dof)
    end
    if refDof < 10
        return paraview_specifics(Int(ceil(dof / 3))) * paraview_specifics(dof - Int(ceil(dof / 3 - 1) * 3))
    end

    @error "not exportable yet as one variable"
    return

end

function get_block_nodes(block_Id, block)
    conn = findall(x -> x == block, block_Id)
    return reshape(conn, 1, length(conn))
end

function init_results_in_exodus(exo, output, coords, block_Id, nsets)
    info = ["PeriLab Version " * string(Pkg.project().version) * ", under BSD License", "Copyright (c) 2023, Christian Willberg, Jan-Timo Hesse", "compiled with Julia Version " * string(VERSION)]
    nsetsExo = read_sets(exo, NodeSet)
    write_info(exo, info)
    # check if type of coords is int or float32
    if typeof(coords) in [Matrix{Int64}, Matrix{Float32}]
        coords = convert(Matrix{Float64}, coords)
    end

    write_coordinates(exo, coords)
    """
        id = 0
        for name in eachindex(nsets)
            id += 1
            nsetExo = NodeSet(id, convert(Array{Int32}, nsets[name]))
            write_set(exo, nsetExo)
            write_name(exo, NodeSet, name)
        end
    """
    for block in sort(unique(block_Id))
        conn = get_block_nodes(block_Id, block)# virtual elements     
        write_block(exo, block, "SPHERE", conn)
        write_name(exo, Block, block, "Block_" * string(block))
    end

    """
    output structure var_name -> [fieldname, exodus id, field dof]
    """
    names = collect(keys(sort(output)))
    write_number_of_variables(exo, NodalVariable, length(names))
    write_names(exo, NodalVariable, names)
    nnodes = exo.init.num_nodes
    for varname in keys(output)
        zero = zeros(Float64, nnodes)
        # interface does not work with Int yet 28//08//2023
        write_values(exo, NodalVariable, 1, output[varname][2], varname, zero)
    end
    return exo
end

function write_step_and_time(exo, step, time)
    write_time(exo, step, Float64(time))
    return exo
end

function write_nodal_results_in_exodus(exo, step, output, datamanager)
    #write_values
    nnodes = datamanager.get_nnodes()
    for varname in keys(output)
        field = datamanager.get_field(output[varname][1])
        #exo, timestep::Integer, id::Integer, var_index::Integer,vector
        # =>https://github.com/cmhamel/Exodus.jl/blob/master/src/Variables.jl  
        var = convert(Array{Float64}, field[1:nnodes, output[varname][3]])
        # interface does not work with Int yet 28//08//2023
        write_values(exo, NodalVariable, step, output[varname][2], varname, var)
    end
    return exo
end

function merge_exodus_file(file_name)
    epu(file_name)
end
end