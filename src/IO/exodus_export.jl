# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Write_Exodus_Results
include("csv_export.jl")
using Exodus
using Pkg
using .Write_CSV_Results
export get_paraviewCoordinates
export create_result_file
export init_results_in_exodus
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
    @debug "num_nodes: $num_nodes"
    @debug "num_elem_blks: $num_elem_blks"
    @debug "num_node_sets: $num_node_sets"
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

function init_results_in_exodus(exo::Exodus.ExodusDatabase, output::Dict{String,Vector{Any}}, computes, coords::Union{Matrix{Int64},Matrix{Float64}}, block_Id::Vector{Int64}, uniqueBlocks::Vector{Int64}, nsets::Dict{String,Vector{Int64}})
    info = ["PeriLab Version " * string(Pkg.project().version) * ", under BSD License", "Copyright (c) 2023, Christian Willberg, Jan-Timo Hesse", "compiled with Julia Version " * string(VERSION)]

    write_info(exo, info)

    # check if type of coords is int or Float64
    if typeof(coords) in [Matrix{Int64}, Matrix{Float64}]
        coords = convert(Matrix{Float64}, coords)
    end

    write_coordinates(exo, coords)

    id::Int32 = 0
    # bloecke checken
    for name in eachindex(nsets)
        id += Int32(1)
        if length(block_Id) < length(nsets[name])
            nsetExo = NodeSet(id, convert(Array{Int32}, nsets[name][1:length(block_Id)]))
        else
            nsetExo = NodeSet(id, convert(Array{Int32}, nsets[name]))
        end
        write_set(exo, nsetExo)
        write_name(exo, nsetExo, name)
    end

    for block in uniqueBlocks
        conn = get_block_nodes(block_Id, block)# virtual elements   
        write_block(exo, block, "SPHERE", conn)
        write_name(exo, Block, block, "Block_" * string(block))
    end

    """
    output structure var_name -> [fieldname, exodus id, field dof]
    """
    output_names = collect(keys(sort(output)))
    compute_names = collect(keys(sort(computes)))

    write_number_of_variables(exo, NodalVariable, length(output_names))
    write_names(exo, NodalVariable, output_names)
    nnodes = exo.init.num_nodes

    for varname in output_names
        # interface does not work with Int yet 28//08//2023
        @debug "$varname $nnodes"
        write_values(exo, NodalVariable, 1, output[varname][2], varname, zeros(Float64, nnodes))
    end
    if length(compute_names) > 0
        compute_field_names = []
        for varname in compute_names
            if haskey(computes[varname], "Mapping")
                compute_field_names = [compute_field_names; collect(keys(sort(computes[varname]["Mapping"])))]
            end
        end
        if length(compute_field_names) > 0
            write_number_of_variables(exo, GlobalVariable, length(compute_field_names))
            write_names(exo, GlobalVariable, Vector{String}(compute_field_names))
            id = 1
            for varname in compute_field_names
                write_values(exo, GlobalVariable, 1, id, varname, [Float64(0.0)])
                id += 1
            end
        end
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

function write_global_results_in_exodus(exo, csv_file, step, computes, datamanager)
    #write_values
    nnodes = datamanager.get_nnodes()
    global_values = []
    for varname in keys(computes)
        if haskey(computes[varname], "Mapping")
            compute_class = computes[varname]["Compute Class"]
            calculation_type = computes[varname]["Calculation Type"]
            field = datamanager.get_field(computes[varname]["Variable"])
            mapping = computes[varname]["Mapping"]
            for field_name in keys(mapping)
                global_value = 0
                values = field[1:nnodes, mapping[field_name]["dof"]]
                if compute_class == "Block_Data"
                    block = computes[varname]["Block"]
                    block_Id = datamanager.get_field("Block_Id")
                    block_filter = filter(x -> block_Id[x] == parse(Int, block[7:end]), 1:length(block_Id))
                    values = values[block_filter]
                end
                if length(values) == 0
                    @warn "No values for $field_name, check compute class parameters or block assignment!"
                    continue
                end
                if calculation_type == "Maximum"
                    global_value = maximum(values)
                elseif calculation_type == "Minimum"
                    global_value = minimum(values)
                elseif calculation_type == "Sum"
                    global_value = sum(values)
                end
                write_values(exo, GlobalVariable, step, mapping[field_name]["result_id"], field_name, [global_value])
                if computes[varname]["CSV Export"]
                    push!(global_values, global_value)
                end
            end
        end
    end
    if length(global_values) > 0
        Write_CSV_Results.write_global_results_in_csv(csv_file, global_values)
    end

    return exo
end

function merge_exodus_file(file_name)
    epu(file_name)
end
end