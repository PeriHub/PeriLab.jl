
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
dof = Ref(1)
fieldnames = Dict(Int64 => Dict{String,Any}(), Float32 => Dict{String,Any}(), Bool => Dict{String,Any}())
##########################
function set_dof(n)
    global dof = n
end
function get_dof()
    return dof
end
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
    return create_field(name * "N", type, "Node_Field", dof), create_field(name * "NP1", type, "Node_Field", dof)
end
function create_constant_node_field(name::String, type::DataType, dof::Int64)
    return create_field(name, type, "Node_Field", dof)
end
function create_bond_field(name::String, type::DataType, dof::Int64)
    return create_field(name * "N", type, "Bond_Field", dof), create_field_bond_type_field(name * "NP1", type, "Bond_Field", dof)
end
function create_constant_bond_field(name::String, type::DataType, dof::Int64)
    return create_field(name, type, "Bond_Field", dof)
end
# wenn nicht existiert, wird an den gesamtvector der Teil angehängt. Dann wird der Vectorabschnitt zurückgeschickt und ist in der jeweiligen routine nutzbar
# bisher wird alles synchrononisiert
# nicht synchronisierte vectoren brauchen nicht hier rein -> Problem ist, wie mit der synch Option zu verfahren ist (eine stelle synch andere non synch)Was zählt?
# besser ist vielleicht einen synchmanager zu bauen -> dort kann man sich eintragen

function create_field(name::String, vartype::DataType, bondOrNode::String, dof::Int64)
    """
    create_field(name::String, vartype::DataType, bondNode::String, dof::Int64)

    Create a field with the given `name` for the specified `vartype`. If the field already exists, return the existing field. If the field does not exist, create a new field with the specified characteristics.

    # Arguments
    - `name::String`: The name of the field.
    - `vartype::DataType`: The data type of the field.
    - `dof::Int64`: The degrees of freedom per node.

    # Returns
    The field with the given `name` and specified characteristics.

    """
    if haskey(fieldnames, vartype) == false
        fieldnames[type] = Dict{String,Any}()
    end
    if haskey(fieldnames[vartype], name)
        return fieldnames[vartype][name]
    end

    if bondOrNode == "Node_Field"
        if dof == 1
            fieldnames[vartype][name] = zeros(vartype, nnodes)
        else
            fieldnames[vartype][name] = zeros(vartype, nnodes, dof)
        end
    elseif bondOrNode == "Bond_Field"
        nBonds = get_field("Number of Neighbors")
        fieldnames[vartype][name] = []
        for i in 1:nnodes
            if dof == 1
                append!(fieldnames[vartype][name], [Vector{vartype}(undef, nBonds[i])])
            else
                append!(fieldnames[vartype][name], [Matrix{vartype}(undef, nBonds[i], dof)])
            end
        end
    end

    return fieldnames[vartype][name]

end


function get_field(name::String, time::String)
    if time == "Constant" || time == "CONSTANT"
        return return_field(name)
    end
    return return_field(name * time)
end

function get_field(name::String)
    return return_field(name)
end

function return_field(name::String)
    """
    get_field(name::String)

    Get the field with the specified `name`. If the field exists, return the field. If the field does not exist, return 0.

    # Arguments
    - `name::String`: The name of the field.

    # Returns
    The field with the given `name` if it exists, otherwise 0.

    """
    for typ in keys(fieldnames)
        if name in keys(fieldnames[typ])
            return fieldnames[typ][name]
        end
    end
    return 0
end

function synch_manager()
    # Liste mit den Daten die synchronisiert werden sollen -> 
    # upload; für init
    # download sum; für time int
    # down - up für bond information
end
end