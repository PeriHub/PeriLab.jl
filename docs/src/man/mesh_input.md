# Mesh Input File Structure

The mesh input file is a text file used to define the geometry of the simulation domain. The file has a header and subsequent lines representing individual nodes with their corresponding coordinates, block IDs, volumes, and optional user-defined values. Below is an example of the structure:

## Header
The header of the mesh input file specifies the format of the subsequent data columns. In this example, the header consists of the following columns:

- `x`: X-coordinate of the node.
- `y`: Y-coordinate of the node.
- `z` (Optional): Z-coordinate of the node (if applicable).
- `block_id`: Block ID to which the node belongs.
- `volume`: Volume associated with the node.
- `user_defined_1`, `user_defined_2`, ... (Optional): User-defined values associated with the node.

## Data Lines
The data lines represent individual nodes in the mesh, with values corresponding to the columns specified in the header.

Example Data Lines:

```plaintext
header: x y block_id volume
0.0 0.0 1 1.0E-02
0.1 0.0 1 1.0E-02
0.1 0.1 1 1.0E-02
...
```

# Nodeset Input File Structure

The nodeset input file is a text file used to define sets of nodes in the simulation domain. The file has a header and subsequent lines representing individual node global IDs that belong to the nodeset. Below is an example of the structure:

## Header
The header of the nodeset input file specifies the format of the subsequent data columns. In this example, the header consists of the following column:

- `global_id`: Global ID of the node.

## Data Lines
The data lines represent individual nodes in the nodeset, with values corresponding to the columns specified in the header.

Example Data Lines:

```plaintext
header: global_id
1
2
11
12
```