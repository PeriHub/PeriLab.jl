# Module integration
If you want to integrate your own model check if it suits in one of the predefined classes material, damage, additive, thermal or degradation. If so check the template folder.

!!! info "Material Template"
    Materials have multiple templates, because the correspondence formulation allows additional options.

Each template has a parameter struct, an init function and a compute function.

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
The init function is called once before the run. The parameters are already checked (types, ranges, required keys); checks that need the mesh or other models belong here, as do the fields your model creates. Every category calls it the same way, `init_model(nodes, p, ctx, block)`; `ctx` is the block's shared category data (`material`, `damage`, `thermal`) or `nothing` (additive, degradation, pre-calculation).

## How PeriLab finds your model
PeriLab loads every file of a model folder that calls the category's `register_*` function (`register_material`, `register_damage`, ...), also in a comment, so the unregistered templates load too. The input deck uses the name passed to `register_*`.

!!! info "Correspondence"
    Put a correspondence model in the folder `Material_Models/Correspondence`. Every material model defined there runs in the correspondence formulation, whatever its name. Correspondence models cannot be combined with other material models by `+`.

## Compute function
This function is called from the solver, `compute_model(nodes, p, ctx, block, time, dt)`. You can call whatever function you like from here. However, this function should evaluate the result needed for the solving process, e.g. heat flux or force densities.

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
using ...ModuleLoader: find_registered_modules
for file in find_registered_modules(@__DIR__, "register_material")
    include(file)
end
```

Each module registers its parameter struct in `__init__()` (`register_material(name, T)`); the factory declares the category with `register_base!` if the category has shared keys. The parsed model of a block is an instance of that struct, so the factory finds the module of the model from its type:

```julia
using ...ParameterSpec: model_parts, model_module
for part in model_parts(material.model)   # the parts of a `+` composite
    model_module(part).compute_model(nodes, part, material, block, time, dt)
end
```
