include("./Parameters/parameter_mesh_handling.jl")
include("./Parameters/parameter_handling_blocks.jl")
function check_element(data, key)
    return haskey(data, key)
end


function create_field(name::String, type::DataType, bondNode::String, dof::Int64, time_dependend::Bool)

    if !haskey(fieldnames, type)
        fieldnames[type] = Dict()
    end
    if haskey(fieldnames[type], name)
        return fieldnames[type][name]
    else
        len = field_length(bondNode, time_dependend) * dof

        if length(values(fieldnames[type])) == 0
            num = 0
        else
            # start ende -> alles in einen Vector, dann sind weniger synchros nötig
            num = maximum(values(fieldnames[type]))[2]
            fieldnames[type][name] = (num + 1, num + len + 1)
        end
        return fieldnames[type][name]
    end
end

function field_length(bondNode::String, time_dependend::Bool)
    le = 0
    if bondNode == "Bond"
        le = nbonds
    elseif bondNode == "Node"
        le = nnodes
    else
        le = 0
        @error "Not supported option $bondNode allowed options are Bond::String or Node::String"
    end
    if time_dependend
        le *= 2
    end
    return le
end
function key_magic(name::String, type, bondNode::String, dof::Int64, time_dependend::Bool)

    #integer liste wo die keys geführt werden -> nötig?
end

function get_field(name::String)
end

function synch_manager()
    # Liste mit den Daten die synchronisiert werden sollen -> 
    # upload; für init
    # download sum; für time int
    # down - up für bond information
end