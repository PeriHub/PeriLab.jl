function get_mesh_name(params)
    check = check_element(params["Discretization"], "Input Mesh File")
    if !check
        @error "No mesh file is defined."
        return
    end
    return params["Discretization"]["Input Mesh File"]
end

function get_topology_name(params)
    check = check_element(params["Discretization"], "Input FEM Topology File")
    topoFile::String = ""
    if check
        topoFile = params["Discretization"]["Input FEM Topology File"]
    end
    return check, topoFile
end

function get_bond_filter(params)
    check = check_element(params["Discretization"], "Bond Filters")
    bfList = Any
    if check
        bfList = params["Discretization"]["Bond Filters"]
    end
    return check, bfList
end