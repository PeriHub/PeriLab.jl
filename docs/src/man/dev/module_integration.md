# Module integration
If you want to integrate your own model check if it suits in one of the predefined classes material, damage, additive, thermal or corrosion. If so check the template folder. 

!!! info "Material Template" 
    Materials have multiple templates, because the correspondence formulation allows additional options.  

Each template has a init function, a name function, a compute function. 

Copy the template and put it in the folder. Change all the functions and give the module a name.

!!! info "Automatic Integration" 
    In PeriLab makros are used to automatically integrate your model.

## Parameter
In PeriLab a field called params exists. This field provides all the material information. The structure is given [here](@ref "Parameters")

## Init function
The init function is used to read and check the  properties provided by the yaml. It should be done there, because if the compute function is used, this check is done in every time step. Also specific fields can be defined here as well.

## Name function
This function defines the name of the module. This name is used to call this model from the yaml.

!!! info "Correspondence" 
    If you want to integrate a correspondence model, make sure ''Correspodence'' occur in the material name

## Compute function
This function is called from the solver. You can call whatever function you like from here. However, this function should evaluate the result needed for the solving process, e.g. heat flux or force densities.

## Module name
You can setup the module name as you like as long as it does not exist a second time in PeriLab.

# Creating your own model category
!!! warn "Creating your own model category" 
    This is advanced programming. Feel free to contact the developers for help. 

To integrate a model category somewhere you have to do the following things. You need a main function of your modeling category. The existing ones are the factory files. These modules have a init function and a compute function. The init function find the modules of the category and the compute function calls these modules during the solving process. 

Here, the call for the init function is shown for the material factory. 

```julia
mod = Set_modules.create_module_specifics(material_model, module_list, "material_name")
datamanager.set_model_module(material_model, mod)
```

The module_list is optained, by applying 

```julia
global module_list = Set_modules.find_module_files(@__DIR__, "material_name")
Set_modules.include_files(module_list)
```

You can integrate these functions than in the compute function of the factory module.

```julia
mod = datamanager.get_model_module(material_model)
datamanager = mod.compute_forces(datamanager, nodes, model_param, time, dt, to)
```
