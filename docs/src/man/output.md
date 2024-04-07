# Output

## Output File Type 
Currently `Exodus` and `CSV` are supported as output types.

!!! warning "CSV File"
    Only variables that are defined as global variables are supported.

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
|Orientations|Rotated Nodes|
|Damage|Damage Model|
|Weighted Volume|PD Solid Elastic|
|Dilatation|PD Solid Elastic|
|Yield Value|PD Solid Plastic|
|Lambda Plastic|PD Solid Plastic|
|Strain|Correspondence|
|Cauchy Stress|Correspondence|
|Zero Energy Stiffness|Correspondence|
|von Mises Stress|Correspondence Plastic|
|Plastic Strain|Correspondence Plastic|
|Temperature|Thermal Models|
|Delta Temperature|Thermal Models|
|Heat Flow|Thermal Models|
|Specific Heat Capacity|Thermal Models|
|Specific Volume|Thermal Models|
|Surface_Nodes|Thermal Models|
|Concentration|Corrosion Models|
|Delta Concentration|Corrosion Models|
|Concentration Flux|Corrosion Models|
|Lumbed Mass Matrix|FEA|
|FE Nodes|FEA|

!!! info "Own Variables"
    All variables that are defined as global variables are supported as well as those defined in the mesh input f.e.: `Angles`
