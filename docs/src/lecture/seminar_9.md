## Seminar 9: PeriLab - Implement your own model

[Plain Model](https://github.com/PeriHub/PeriLab.jl/tree/main/examples/Seminars/Part_09)

### Steps
- choose your model category
- take the template file, rename it and copy it in the model folder
- give you model a name, this name defines how you call your model

```julia
function damage_name()
    return "My Model"
end

```

```yaml
 Models:
    Damage Models:
      Damage:
        Damage Model: My Model
```

- Specify your model parameter

```yaml
 Models:
    Damage Models:
      Damage:
        Damage Model: My Model
        Important value: 200
```

- get this value in the code


```julia
function init_model(nodes::AbstractVector{Int64},
                    damage_parameter::Dict,
                    block::Int64)

    println(damage_parameter["Important value"])
end
```

- init_model is used to create model specific fields and to check if values, especially optional values exist

```julia
function init_model(nodes::AbstractVector{Int64},
                    damage_parameter::Dict,
                    block::Int64)
    if !inothing(get(damage_parameter, "Important value", nothing))
        println(damage_parameter["Important value"])
    else
        damage_parameter["Important value"] = 0
    end
    my_constant_node_field = Data_Manager.create_constant_node_field("my constant node field", Float64, 10)
    my_node_field_N, my_node_field_NP1 = Data_Manager.create_node_field("my node field", Float64, 2, VectorOrMatrix="Matrix")
    my_constant_bond_field = Data_Manager.create_constant_bond_field("my constant bond field", Float64, 10)
    my_bobd_field_N,  my_bobd_field_NP1 = Data_Manager.create_bond_field("my bond field", Float64, 10)

   field = Data_Manager.create_constant_free_size_field(Fieldname::String, Type_of_variable::Type, size::NTuple)
    fieldN, fieldNP1 = Data_Manager.create_free_size_field(Fieldname::String, Type_of_variable::Type, size::NTuple)

   my_free_size_field = Data_Manager.create_constant_free_size_field("my free size field", Bool, (200,1,3,4,1))
end
```

- If fields already exist, the field is returned
- Get fields and use them; if you want to now, what fields are already defined and usable

```julia
Data_Manager.get_all_field_keys()
```

- All node fields can be exported to the result file
```julia
function compute_model(nodes::AbstractVector{Int64},
                       damage_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64,
                       to::TimerOutput)
    my_constant_field = Data_Manager.get_field("my constant node field")
    my_field_N = Data_Manager.get_field("my node field","N")
    my_field_NP1 = Data_Manager.get_field("my node field","NP1")

    my_field_NP1 .= damage_parameter["Important value"]
end
```
- If N and NP1 exist only NP1 will be exported

```yaml
    Output1:
      Number of Output Steps: 100
      Output File Type: Exodus
      Output Filename: Output_file_name
      Output Variables:
        my node field: true
        my constant node field: true
```
- Activate the model class. Material is set to true by default. The rest is set to zero.


```yaml
    Solver:
        Damage Models: true
```


### Run your model
1. create your yaml with all your parameters
2. create a mesh file
3. run the model

```julia
using PeriLab
PeriLab.main("Folder where to find your yaml")
```

### Remarks
If you define parameter in your mesh, you can call them in PeriLab as well.

```ascii
header: x y block_id volume my_values my_datax my_datay
```

```julia
Data_Manager.get_field("my_values")
Data_Manager.get_field("my_data")
```
my_datax and my_datay are a 2D array.

!!! info "Mesh defined field properties"
    The variable type is given by the input (Int64, Float64 or Bool) and its a constant field.
