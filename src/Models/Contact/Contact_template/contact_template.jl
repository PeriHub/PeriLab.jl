# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact_template

using .....Data_Manager
using .....ParameterSpec: @params, register_contact

export contact_model_name
export init_contact_model
export compute_contact_model

"""
    ContactTemplateParams

Declare the YAML keys your contact model needs beyond the shared contact keys
(Contact Radius, Symmetry, Contact Groups are in `contact.base`). Example:

    my_parameter::Float64 = req("My Parameter"; min = 0, description = "...")

Register it under the name the input deck uses (`Type`) by uncommenting
`__init__` (the template stays unregistered so that a copy never collides with it).
"""
@params struct ContactTemplateParams
end

# __init__() = register_contact("Contact Template", ContactTemplateParams)

"""
    contact_model_name()

Gives the contact model name. PeriLab loads the module because it defines this function; the input deck uses the name passed to `register_contact` in `__init__()`.

# Arguments

# Returns
- `name::String`: The name of the contact model.

Example:
```julia
println(contact_model_name())
"Contact Template"
```
"""
function contact_model_name()
    return "Contact Template"
end

"""
    init_contact_model(p, contact)

Inits the contact model. This template has to be copied, the file renamed and edited by the user to create a new contact. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `p::ContactTemplateParams`: The model parameters.
- `contact`: The contact model of the input deck; `contact.base` holds the shared contact keys.
"""
function init_contact_model(p::ContactTemplateParams, contact)
end

"""
    compute_contact_model(cg, p, contact, compute_master_force_density, compute_slave_force_density)

Computes the contact forces of contact group `cg`. The contact pairs of the group are in `Data_Manager.get_contact_dict(cg)`.

# Arguments
- `cg::String`: The contact group.
- `p::ContactTemplateParams`: The model parameters.
- `contact`: The contact model of the input deck; `contact.base` holds the shared contact keys.
- `compute_master_force_density::Function`: Adds a force to a master node, `(master_id, slave_id, force)`.
- `compute_slave_force_density::Function`: Adds a force to a slave node, `(slave_id, master_id, force)`.
"""
function compute_contact_model(cg, p::ContactTemplateParams, contact,
                               compute_master_force_density::Function,
                               compute_slave_force_density::Function)
    @info "Please register your contact model with register_contact in __init__()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_contact_model(cg, p, contact, compute_master_force_density, compute_slave_force_density) function."
    @info "The Data_Manager, p and contact hold all you need to solve your problem on contact level."
    @info "Add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
end

end
