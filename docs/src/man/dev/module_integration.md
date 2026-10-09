# Module integration
If you want to integrate your own model check if it suits in one of the predefined classes material, damage, additive, thermal or degradation. If so check the template folder.

!!! info "Material Template"
    Materials have multiple templates, because the correspondence formulation allows additional options.

Each template has a parameter struct, a init function, a name function and a compute function.

Copy the template and put it in the folder. Change all the functions and give the module a name.

!!! info "Automatic Integration"
    In PeriLab makros are used to automatically integrate your model.

## Parameter
Your model gets its parameters as a typed struct `p`, read and checked when the input deck is read. How parameters are declared is described [here](@ref "Parameters").

## Declare and register your parameters
Every YAML key your model reads must be declared, otherwise the input deck is rejected (unknown key). Declare the keys in an `@params` struct and register it under the name the input deck uses:

```julia
using ......ParameterSpec: @params, register_material

@params struct MyMaterialParams
    my_parameter::Float64 = req("My Parameter"; min = 0, description = "...")
end
__init__() = register_material("My Material", MyMaterialParams)
```

Material models do not repeat the shared material keys (Symmetry, Young's Modulus, Bulk Modulus, ...); they are declared once in the material factory. The material template contains the struct with the registration commented out: rename both and uncomment it. A model that is not registered is reported as "not found" when the input deck is read.

Contact models work the same way: they register with `register_contact` under the name used in `Type`; the shared contact keys (Contact Radius, Symmetry, Contact Groups) are in `contact.base`, see the contact template.

## Init function
The init function is called once before the run. The parameters are already checked (types, ranges, required keys); checks that need the mesh or other models belong here, as do the fields your model creates. The template of each category shows the arguments, e.g. `init_model(nodes, p, material)` for materials.

## Name function
PeriLab loads every file of a model folder that defines the category's name function (`material_name()`, `damage_name()`, ...). The input deck uses the name passed to `register_*`.

!!! info "Correspondence"
    If you want to integrate a correspondence model, make sure "Correspondence" occurs in the registered model name

## Compute function
This function is called from the solver. You can call whatever function you like from here. However, this function should evaluate the result needed for the solving process, e.g. heat flux or force densities.

## Module name
You can setup the module name as you like as long as it does not exist a second time in PeriLab.

!!! info "Coding Style"
    Please name the module file and the module equaly.

# Creating your own model category
!!! warn "Creating your own model category"
    This is advanced programming. Feel free to contact the developers for help.

To integrate a model category somewhere you have to do the following things. You need a main function of your modeling category. The existing ones are the factory files. These modules have a init function and a compute function. The init function find the modules of the category and the compute function calls these modules during the solving process.

The model modules are found by applying

```julia
using ...ModuleLoader: find_module_files, create_module_specifics
global module_list = find_module_files(@__DIR__, "material_name")
for mod in module_list
    include(mod["File"])
end
```

Each module registers its parameter struct in `__init__()` (`register_material(name, T)`); the factory declares the category with `register_base!` if the category has shared keys. The parsed model of a block is an instance of that struct, so the factory finds the module of the model from its type:

```julia
mod = parentmodule(typeof(p))
mod.compute_model(nodes, p, material, block, time, dt)
```
