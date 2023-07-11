
module Data_manager
include("./Parameters/parameter_handling.jl")
#export nnodes
#export nbonds
export set_nnodes
export get_nnodes
export set_nbonds
export get_nbonds
export create_node_field
export create_constant_node_field
export create_bond_field
export create_constant_bond_field
##########################
# Variables
##########################
nnodes = Ref(0)
nbonds = Ref(0)
fieldnames = Dict(Int64 => Dict(), Float32 => Dict(), Bool => Dict())
##########################


function set_nnodes(n)
    global nnodes = n
end
function get_nnodes()
    return nnodes
end

function set_nbonds(n)
    global nbonds = n
end
function get_nbonds()
    return nbonds
end
function create_node_field(name::String, type::DataType, dof::Int64)
    return create_field(name, type, "Node", dof, true)
end
function create_constant_node_field(name::String, type::DataType, dof::Int64)
    return create_field(name, type, "Node", dof, false)
end
function create_bond_field(name::String, type::DataType, dof::Int64)
    return create_field(name, type, "Bond", dof, true)
end
function create_constant_bond_field(name::String, type::DataType, dof::Int64)
    return create_field(name, type, "Bond", dof, false)
end
# wenn nicht existiert, wird an den gesamtvector der Teil angehängt. Dann wird der Vectorabschnitt zurückgeschickt und ist in der jeweiligen routine nutzbar
# bisher wird alles synchrononisiert
# nicht synchronisierte vectoren brauchen nicht hier rein -> Problem ist, wie mit der synch Option zu verfahren ist (eine stelle synch andere non synch)Was zählt?
# besser ist vielleicht einen synchmanager zu bauen -> dort kann man sich eintragen

function create_field(name::String, vartype::DataType, bondNode::String, dof::Int64, time_dependend::Bool)
    """
    create_field(name::String, vartype::DataType, bondNode::String, dof::Int64, time_dependend::Bool)

    Create a field with the given `name` for the specified `vartype`. If the field already exists, return the existing field. If the field does not exist, create a new field with the specified characteristics.

    # Arguments
    - `name::String`: The name of the field.
    - `vartype::DataType`: The data type of the field.
    - `bondNode::String`: length of the field (nnodes or nbonds) associated with the field.
    - `dof::Int64`: The degrees of freedom per node.
    - `time_dependend::Bool`: Indicates whether the field is time-dependent.

    # Returns
    The field with the given `name` and specified characteristics.

    """
    if haskey(fieldnames[vartype], name)
        return fieldnames[vartype][name]
    else
        len = field_length(bondNode, time_dependend)
        fieldnames[vartype][name] = zeros(vartype, len, dof)
        return fieldnames[vartype][name]
    end
end

function get_field(name::String)
    for typ in keys(fieldnames)
        if name in keys(fieldnames[typ])
            return fieldnames[typ][name]
        end
    end
    return 0
end

function field_length(bondNode::String, time_dependend::Bool)
    le = 0
    if bondNode == "Bond"
        le = get_nbonds()
    elseif bondNode == "Node"
        le = get_nnodes()
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

function synch_manager()
    # Liste mit den Daten die synchronisiert werden sollen -> 
    # upload; für init
    # download sum; für time int
    # down - up für bond information
end
end