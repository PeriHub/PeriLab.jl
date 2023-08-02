
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
glob_to_loc = Ref([])
fields = Dict(Int64 => Dict{String,Any}(), Float32 => Dict{String,Any}(), Bool => Dict{String,Any}())
fieldnames = Dict{String,DataType}()
filtered_nodes = Ref([])
##########################
# Material information
##########################
material_type = Dict{String,Bool}("Bond-Based" => false, "Ordinary" => false, "Correspondence" => true, "Bond-Associated" => false)
##########################
function create_bond_field(name::String, type::DataType, dof::Int64)
    """
    create_bond_field(name::String, type::DataType, dof::Int64)

    Creates a bond field with the given name, data type, and degree of freedom.

    Parameters:
    - `name` (string): The name of the bond field.
    - `type` (DataType): The data type of the bond field.
    - `dof` (Int64): The degree of freedom of each bond.

    Returns:
    - `bond_field` (Field): The created bond field for the current time step.
    - `bond_field_np1` (Field): The created bond field for the next time step.

    Example:
    ```julia
    create_bond_field("stress", Float64, 6)  # creates a stress bond field with 6 degrees of freedom
    ```
    """
    return create_field(name * "N", type, "Bond_Field", dof), create_field(name * "NP1", type, "Bond_Field", dof)
end

function create_constant_bond_field(name::String, type::DataType, dof::Int64)
    """
    create_constant_bond_field(name::String, type::DataType, dof::Int64)

    Creates a constant bond field with the given name, data type, and degree of freedom.

    Parameters:
    - `name` (string): The name of the constant bond field.
    - `type` (DataType): The data type of the constant bond field.
    - `dof` (Int64): The degree of freedom for each bond.

    Returns:
    - `constant_bond_field` (Field): The created constant bond field.

    Example:
    ```julia
    create_constant_bond_field("density", Float64, 1)  # creates a density constant bond field
    ```
    """
    return create_field(name, type, "Bond_Field", dof)
end

function create_constant_node_field(name::String, type::DataType, dof::Int64)
    """
    create_constant_node_field(name::String, type::DataType, dof::Int64)

    Creates a constant node field with the given name, data type, and degree of freedom.

    Parameters:
    - `name` (string): The name of the constant node field.
    - `type` (DataType): The data type of the constant node field.
    - `dof` (Int64): The degree of freedom of each node.

    Returns:
    - `constant_node_field` (Field): The created constant node field.

    Example:
    ```julia
    create_constant_node_field("temperature", Float64, 1)  # creates a temperature constant node field
    ```
    """
    return create_field(name, type, "Node_Field", dof)
end

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
    if haskey(fields, vartype) == false
        fields[vartype] = Dict{String,Any}()
    end
    if name in get_all_field_keys()
        return get_field(name)
    end

    if bondOrNode == "Node_Field"
        if dof == 1
            fields[vartype][name] = zeros(vartype, nnodes)
        else
            fields[vartype][name] = zeros(vartype, nnodes, dof)
        end
    elseif bondOrNode == "Bond_Field"
        nBonds = get_field("Number of Neighbors")
        fields[vartype][name] = []
        for i in 1:nnodes
            if dof == 1
                append!(fields[vartype][name], [Vector{vartype}(undef, nBonds[i])])
            else
                append!(fields[vartype][name], [Matrix{vartype}(undef, nBonds[i], dof)])
            end
        end
    end
    global fieldnames[name] = vartype
    return fields[vartype][name]

end

function create_node_field(name::String, type::DataType, dof::Int64)
    """
    create_node_field(name::String, type::DataType, dof::Int64)

    Creates a node field with the given name, data type, and degree of freedom.

    Parameters:
    - `name` (string): The name of the node field.
    - `type` (DataType): The data type of the node field.
    - `dof` (Int64): The degree of freedom of each node.

    Returns:
    - `node_field` (Field): The created node field for the current time step.
    - `node_field_np1` (Field): The created node field for the next time step.

    Example:
    ```julia
    create_node_field("displacement", Float64, 3)  # creates a displacement node field with 3 degrees of freedom
    ```
    """
    return create_field(name * "N", type, "Node_Field", dof), create_field(name * "NP1", type, "Node_Field", dof)
end

function get_all_field_keys()
    return keys(fieldnames)
end

function get_dof()
    """
    get_dof()

    Retrieves the degree of freedom (dof) value.

    Returns:
    - `dof` (integer): The current degree of freedom value.

    Example:
    ```julia
    get_dof()  # returns the current degree of freedom
    ```
    """
    return dof
end

function get_local_nodes(global_nodes)
    """
    get_local_nodes()

    Determines the local node numbering.

    Returns:
    - `get_local_nodes` (array): returns local nodes.

    Example:
    ```julia
    get_local_nodes()  # returns local nodes
    ```
    """
    return [glob_to_loc[global_node] for global_node in global_nodes]

end

function get_nnodes()
    """
    get_nnodes()

    Retrieves the number of nodes.

    Returns:
    - `nnodes` (integer): The current number of nodes.

    Example:
    ```julia
    get_nnodes()  # returns the current number of nodes
    ```
    """
    return nnodes
end

function get_field(name::String, time::String)

    if time == "Constant" || time == "CONSTANT"
        return get_field(name)
    end
    return get_field(name * time)

end

function get_field(name::String)
    if name in get_all_field_keys()
        return fields[fieldnames[name]][name]
    end
    return []
end

function get_material_type(key)
    return material_type[key]
end

function get_nbonds()
    """
    get_nbonds()

    Retrieves the number of bonds.

    Returns:
    - `nbonds` (integer): The current number of bonds.

    Example:
    ```julia
    get_nbonds()  # returns the current number of bonds
    ```
    """
    return nbonds
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
    for typ in keys(fields)
        if name in keys(fields[typ])
            if length(fields[typ][name][1, :]) == 1
                # avoiding vector to matrix type transformation
                return fields[typ][name][filtered_nodes]
            end
            return fields[typ][name][filtered_nodes, :]
        end
    end
    return 0
end

function set_dof(n)
    """
    set_dof(n)

    Sets the degree of freedom (dof) value globally.

    Parameters:
    - `n` (integer): The value to set as the degree of freedom.

    Example:
    ```julia
    set_dof(3)  # sets the degree of freedom to 3
    ```
    """
    global dof = n
end

function set_filter(list_of_nodes)
    global filtered_nodes = list_of_nodes
end

function set_glob_to_loc(list)
    """
    set_glob_to_loc(list)

    Sets the global-to-local mapping list globally.

    Parameters:
    - `list` (array): The list representing the global-to-local mapping.

    Example:
    ```julia
    set_glob_to_loc([1, 3, 5])  # sets the global-to-local mapping list
    ```
    """
    global glob_to_loc = list
end

function set_material_type(key, value)
    if key in keys(material_type)
        global material_type[key] = value
    else
        @warn "Type " * key * " is not defined"
    end
end

function set_nnodes(n)
    """
    set_nnodes(n)

    Sets the number of nodes globally.

    Parameters:
    - `n` (integer): The value to set as the number of nodes.

    Example:
    ```julia
    set_nnodes(10)  # sets the number of nodes to 10
    ```
    """
    global nnodes = n
end

function set_nbonds(n)
    """
    set_nbonds(n)

    Sets the number of bonds globally.

    Parameters:
    - `n` (integer): The value to set as the number of bonds.

    Example:
    ```julia
    set_nbonds(20)  # sets the number of bonds to 20
    ```
    """
    global nbonds = n
end

function get_NP1_to_N_Dict()
    NP1_to_N = Dict{String,String}()
    for key in get_all_field_keys()
        if occursin("NP1", key)
            NP1_to_N[key] = key[1:end-2]
        end
    end
    return NP1_to_N
end


function switch_NP1_to_N()
    NP1_to_N = get_NP1_to_N_Dict()
    for NP1 in keys(NP1_to_N)
        field_NP1 = get_field(NP1)
        temp_field_name = NP1_to_N[NP1]
        field_N = get_field(temp_field_name)
        field_N[:] = field_NP1[:]
    end
end
function synch_manager()
    # Liste mit den Daten die synchronisiert werden sollen -> 
    # upload; für init
    # download sum; für time int
    # down - up für bond information
end
end