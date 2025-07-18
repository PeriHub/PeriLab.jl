# Output

## Output File Type
Currently `Exodus` and `CSV` are supported as output types.

!!! warning "CSV File"
    Only variables that are defined as global variables are supported, see Compute Class Parameters for more information.

!!! tip "Flush File"
    By default the output files will be flushed after each write-step, meaning you can look at the results while the simulation is still running.

## Output Frequency / Number of Output Steps
As the number of solver steps can be quite large and you don't want to buy new drives for every simulation we want to specify the number of output steps written.

You can either specify the total number of output steps via `Number of Output Steps` or the frequency of output via `Output Frequency`.

!!! tip "First time users"
    For first time users we will recommend to use `Number of Output Steps` in the range of `100` to `1000`. Depending on your discretization this is enough to get a good idea of the simulation results.

!!! tip "Only damage"
    If you want to take a closer look at damage initiation and propagation and are worried about the result file size, you can write another output file with `Write After Damage` set to `true`.

## Output Variables
Below you can find a list of all current Output Variables and the necessary Prerequisites.

| Variable | Prerequisite |
|---|:---:|
|Density|-|
|Horizon|-|
|Number of Neighbors|-|
|Number of Filtered Neighbors|Bond-Filter|
|Update List|-|
|Active|-|
|Displacements|-|
|Velocity|-|
|Acceleration|-|
|Forces|-|
|Force Densities|-|
|Cauchy Stress|-|
|von Mises Stress|-|
|Orientations|Rotated Nodes|
|Damage|Damage Model|
|Weighted Volume|PD Solid Elastic|
|Dilatation|PD Solid Elastic|
|Yield Value|PD Solid Plastic|
|Lambda Plastic|PD Solid Plastic|
|Strain|Correspondence|
|Zero Energy Stiffness|Correspondence|
|Plastic Strain|Correspondence Plastic|
|Temperature|Thermal Models|
|Delta Temperature|Thermal Models|
|Heat Flow|Thermal Models|
|Specific Heat Capacity|Thermal Models|
|Specific Volume|Thermal Models|
|Surface_Nodes|Thermal Models|
|Concentration|Degradation Models|
|Delta Concentration|Degradation Models|
|Concentration Flux|Degradation Models|
|Lumbed Mass Matrix|FEA|
|FE Nodes|FEA|

!!! info "Own Variables"
    All variables that are defined as global variables are supported as well as those defined in the mesh input f.e.: `Angles`

## Compute Class Parameters

In order to compute output values for a defined nodesets or blocks, especially for the CSV output format, you can use the following parameters:

```yaml
Compute Class Parameters:
    Block_Forces:
        Compute Class: Block_Data
        Calculation Type: Sum
        Block: block_name
        Variable: Forces
    NodeSet_Forces:
        Compute Class: Node_Set_Data
        Calculation Type: Average
        Node Set: nodeset_name
        Variable: Forces
```
Available compute classes are "Block_Data" or "Node_Set_Data".

Supported calculation types are "Sum", "Average", "Maximum", and "Minimum".

!!! tip "Referencing the compute class"
    The created variable can be used for CSV as well as Exodus Outputs by referencing it's name in the Output Variable list.
