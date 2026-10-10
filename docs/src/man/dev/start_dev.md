## Getting started - Develop your own model

PeriLabs main development focus right now is the model section. Each model have folders with templates. Here you can find a plain .jl file with all functions you need for a minimal configuration.

This guide shows how to set up a material model. But it's the same strategy for more or less all models. The key is what you have to compute. In material it is the bond force in deformed bond direction. In damage models it is the bond damage. For thermal models it is the flux and / or the thermal deformation gradient.

### Model integration

This is how the module header looks like. You can add here additional functions you want to export or load additional packages if needed. This set up is the basis.


```julia
module Material_template
using ......Data_Manager
using ......ParameterSpec: @params, register_material
export fe_support
export init_model
export compute_model
export fields_for_local_synchronization

...

end
```
First you have to rename your module. It should be different from existing ones. The first character should be capital. The rest is your choice.

Next you declare the YAML keys your model reads in an `@params` struct and register it under the name the input deck uses (see [Module integration](module_integration.md)):

```julia
@params struct MyModelParams
    my_parameter::Float64 = req("My Parameter"; min = 0, description = "...")
end
__init__() = register_material("My model", MyModelParams)
```

PeriLab loads every file of the model folder that calls `register_material` (also in a comment, so the unregistered template loads too).

In the YAML file you can call your model like this. The id Mat_1 has to be applied to a block and your module works

```yaml
 Material Models:
      Mat_1:
        Material Model: "My model"
        My Parameter: 1.0
```
The density and the horizon depends on your problem.

!!! info "Mass density"
    The parameter density is defined in the blocks. Reason is that this parameter is used for other models as thermal as well.

```yaml
    block_1:
      Block ID: 1
      Density: 1.95e-06
      Horizon: 0.82
      Material Model: Mat_1
```
Congratulations you first model works and runs. However, it does not do anything.

`PeriLab.describe("My model")` prints the keys of your model, `PeriLab.describe("My model"; template = true)` a YAML block to copy into an input deck.

### Include features


The module TimerOutputs allows you measuring the computation speed and memory for certain functions or areas of your code.

```julia
@timeit "foo function" foo()

@timeit "more foo functions" begin
    foo1()
    foo2()
end
```

The function init_model is used to initialize the model. Here, you can set reference values. `p` is your parameter struct, read and checked when the input deck was read; `material` holds the shared material keys of the block (`material.base`), the elastic moduli (`material.moduli`) and the symmetry (`material.symmetry`).

You can and should also create fields using the [datamanger](datamanager.md).

```julia
function init_model(nodes::AbstractVector{Int64},
                    p::MyModelParams,
                    material,
                    block::Int64)
end
```
This function is used as marker if it is suited for FE support. Therefore, it must be able to compute a Cauchy stress.

```julia
function fe_support()
    return false
end
```

This function is used to specify properties which should be synchronized before the material models are called. An example is, that you need the neighboring deformation gradients to compute the bond associated correspondence formulation. Therefore, you have to synchronize.

This feature is only needed for multicore applications. The main parameters are synchronized at the beginning of each time step.

!!! warning "Cost of synchronization"
    Be aware, that synchronization cost time. Therefore, try to avoid it if possible.

```julia
function fields_for_local_synchronization(model::String)
end
```

The last module computes the model in each time step. You get all nodes and have to compute the neighorhood of it. You must call the fields you need from the datamanager. The parameters are fields of `p` (e.g. `p.my_parameter`); they are checked before the run starts.


!!! info "Separation of function"
    Please seperate the functionality. Call type stable functions for analysis. Use compute_model as function to structure your analysis (getting data, call subfunctions). It is better to read, better to test and better to optimize.

```julia
function compute_model(nodes::AbstractVector{Int64},
                       p::MyModelParams,
                       material,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
end
```

Every model category calls its models the same way:

```julia
init_model(nodes, p, ctx, block)
compute_model(nodes, p, ctx, block, time, dt)
```

`ctx` holds what the models of a category share in a block: `material` for materials, `damage` for damage models (shared keys in `damage.base`), `thermal` for thermal models (`thermal.base`). Additive, degradation and pre-calculation models have no shared keys; they get `nothing`. Two kinds of models differ: correspondence models compute stresses (`compute_stresses`, see the correspondence template), and contact models work on contact groups (see the contact template).

### Additional functions
You can add as many functions as you want. You can also create files and include them if needed. As long as the interface functions exists, the module can be called.


### Template files
#### Physical templates

[Additive](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Additive/Additive_template/additive_template.jl)

[Contact](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Contact/Contact_template/contact_template.jl)

[Damage](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Damage/Damage_template/damage_template.jl)

[Non correspondence material](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Material_template/material_template.jl)

[Correspondence material](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Material_template/correspondence_template.jl)

[Pre calculation](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Pre_calculation/Pre_calculation_template/pre_calculation_template.jl)

[Thermal](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Thermal/Thermal_template/thermal_template.jl)

#### Mesh based templates

[Bond Filter](https://github.com/PeriHub/PeriLab.jl/blob/main/src/IO/Bond_Filter/filter_template.jl)
