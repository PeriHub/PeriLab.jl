# Mesh Input File Structure

The mesh input file is a text file used to define the geometry of the simulation domain. The file has a header and subsequent lines representing individual nodes with their corresponding coordinates, block IDs, volumes, and optional user-defined values. Below is an example of the structure:

## Header
The header of the mesh input file specifies the format of the subsequent data columns. In this example, the header consists of the following columns:

**Variable definition**

| Parameter | Datamanager name | Header name 2D | Header name 3D | Type | 
|:---|:---|:---|:---|:---|
|x, y, z (optional) coordinate of the node | Coordinates | x, y | x, y, z | Float64, Int64|
| Definition to which block the node corresponds. Is needed in the Yaml file to define properties | Block_Id | block_id | block_id | Int64|
| Volume the node represents. | Volume | volume | volume | Float64, Int64|


**Optional**

| Parameter | Datamanager name | Header name 2D | Header name 3D | Type | 
|:---|:---|:---|:---|:---|
| Orientation of a Node | Angles | Angles | Anglesx, Anglesy, Anglesz | Float64, Int64|
| Activation time of a node, e.g. used for additive manufacturing to define when the node will be acativated | Activation_Time | Activation_Time | Activation_Time | Float64, Int64|
| Status of the node. If it is false the node is deactivated, but exists. This variable is automatically created if additive models are used and set everywhere to false, if it is not predefined | Active | Active | Active | Bool |


The difference between 2D and 3D is found automatically. If no z occurs PeriLab identifies it as 2D problem and requests a plane stress or plane strain definition.


!!! tip "Additional parameter"
    Additional parameter can be applied in the header. They will be added in the datamanager and can be used in the programm. If you add x,y,z to the parameter a multidimensional field will be created, e.g.
    MyVarx, MyVary will be created as field MyVar with 2 degrees of freedom

### Data Lines
The data lines represent individual nodes in the mesh, with values corresponding to the columns specified in the header.

Example Data Lines:

```plaintext
header: x y block_id volume
0.0 0.0 1 1.0E-02
0.1 0.0 1 1.0E-02
0.1 0.1 1 1.0E-02
...
```

# Abaqus Input

You can use the Abaqus input file to define the geometry of the simulation domain. In order to do that, refer in the input deck to your .inp file:

```yaml
PeriLab:
  Discretization:
    Input Mesh File: ABAQUS_FILE.inp
    Type: Abaqus
...
```

All elements that are defined in a element set in the Abaqus input file will be translated to PeriLab nodes. The center and volume of the elements will be calculated automatically. Have a look at the [AbaqusReader.jl](https://github.com/JuliaFEM/AbaqusReader.jl) package to see what elements are supported.

!!! warning "Supported elements"
    Currently only Hex8 elements are tested!

### How to define blocks and nodesets with Abaqus?

The element sets are defined in the Abaqus input file and can be used to define blocks and nodesets. The order of the blocks will be similar to the order in the .inp file. Nodesets can be referenced via the correspoinding element set in the Abaqus input file.

!!! tip "Block order"
    If you are not sure what order the blocks in the .inp file will be read in, you can use PeriLab to create an exodus file and check the order of the blocks in ParaView.

# Nodeset Input File Structure

The nodeset input file is a text file used to define sets of nodes in the simulation domain. The file has a header and subsequent lines representing individual node global IDs that belong to the nodeset. Below is an example of the structure:

## Header
The header of the nodeset input file specifies the format of the subsequent data columns. In this example, the header consists of the following column:

- `global_id`: Global ID of the node.

### Data Lines
The data lines represent individual nodes in the nodeset, with values corresponding to the columns specified in the header.

Example Data Lines:

```plaintext
header: global_id
1
2
11
12
```