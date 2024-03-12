# PeriLab Configuration File

The PeriLab configuration file is a YAML file used to specify the parameters for running simulations in the PeriLab software. This structure represents the global configuration for PeriLab.

## `PeriLab`

- **Blocks**: List of block configurations.
  - **_block_id_**: Block-specific parameters.
    - Density: Numeric value representing density. (Float64 or Int64)
    - Horizon: Numeric value representing horizon. (Float64 or Int64)
    - Specific Heat Capacity: Numeric value representing specific heat capacity. (Float64 or Int64, optional)
    - Material Model: String representing the material model. (optional)
    - Damage Model: String representing the damage model. (optional)
    - Thermal Model: String representing the thermal model. (optional)
    - Additive Model: String representing the additive model. (optional)

## `FEM (optional)`

- **Element Type**: String representing the type of finite element.
- **Degree**: Numeric value representing the degree. (String or Int64)
- **Material Model**: String representing the material model.

## `Boundary Conditions (optional)`

- **_Own_Name_**: List of boundary condition configurations.
  - **Coordinate**: String representing the coordinate.
  - **Node Set**: String representing the node set.
  - **Type**: String representing the type of boundary condition.
  - **Value**: Numeric value, string, or a combination representing the boundary condition value. (Float64, Int64, String)

## `Compute Class Parameters (optional)`

- **_Own_Name_**: List of compute class parameters.
  - **Block**: String representing the block.
  - **Node Set**: String representing the node set.
  - **Calculation Type**: String representing the calculation type.
  - **Compute Class**: String representing the compute class.
  - **Variable**: String representing the variable.

## `Discretization`

- **Input Mesh File**: String representing the input mesh file.
- **Input External Topology**: External topology parameters. (optional)
  - **File**: String representing the external topology file.
  - **Add Neighbor Search**: Boolean indicating whether to add neighbor search. (optional)
- **Node Sets**: Dictionary of node sets. (optional)
  - **_Own_Name_**: Numeric value or string representing the node set. (Int64 or String)
- **Type**: String representing the type of discretization.
- **Distribution Type**: String representing the distribution type. (optional)
- **Surface Extrusion**: Extrusion parameters. (optional)
  - **Direction**: String representing the extrusion direction.
  - **Step**: Numeric value representing the extrusion step. (Float64 or Int64)
  - **Number**: Numeric value representing the extrusion number. (Float64 or Int64)
- **Bond Filters**: List of bond filter configurations. (optional)
  - **_Own_Name_**: Bond filter parameters.
    - **Type**: String representing the bond filter type.
    - **Normal X/Y/Z**: Numeric values representing normal components. (Float64 or Int64)
    - **Lower Left Corner X/Y/Z**: Numeric values representing lower-left corner components. (Float64 or Int64)
    - **Bottom Unit Vector X/Y/Z**: Numeric values representing bottom unit vector components. (Float64 or Int64)
    - **Center X/Y/Z**: Numeric values representing center components. (Float64 or Int64)
    - **Radius**: Numeric value representing the radius. (Float64 or Int64)
    - **Bottom Length**: Numeric value representing the bottom length. (Float64 or Int64)
    - **Side Length**: Numeric value representing the side length. (Float64 or Int64)
    - **Allow Contact**: Boolean indicating whether to allow contact. (optional)

## `Outputs (optional)`

- **_Own_Name_**: List of output configurations.
  - **Flush File**: Boolean indicating whether to flush the output file. (optional)
  - **Output Frequency**: Integer representing the output frequency.
  - **Number of Output Steps**: Integer representing the number of output steps.
  - **Output File Type**: String representing the output file type. (optional)
  - **Output Filename**: String representing the output filename.
  - **Write After Damage**: Boolean indicating whether to write after damage. (optional)
  - **Output Variables**: Dictionary of output variables.
    - **_Own_Name_**: Boolean indicating whether to output the variable.

## `Physics (optional)`

- **Damage Models**: Dictionary of damage models. (optional)
  - **_Own_Name_**: List of damage model configurations.
    - **Critical Value**: Numeric value representing the critical value. (Float64 or Int64)
    - **Damage Model**: String representing the damage model.
    - **Interblock Damage**: Dictionary of interblock damage parameters.
      - **_Own_Name_**: Numeric value representing interblock damage. (Float64 or Int64, required)
    - **Anisotropic Damage**: Dictionary of anisotropic damage parameters.
      - **Critical Value X/Y**: Numeric values representing critical values in X and Y directions. (Float64 or Int64, required)
- **Material Models**: Dictionary of material models. (optional)
  - **_Own_Name_**: List of material model configurations.
    - **Material Model**: String representing the material model.
    - **Symmetry**: String representing the symmetry. (optional)
    - **Poisson's Ratio/Young's Modulus/Bulk Modulus/Shear Modulus**: Numeric values representing material properties. (Float64 or Int64, optional)
    - **Yield Stress**: Numeric value representing the yield stress. (Float64 or Int64, optional)
    - **Zero Energy Control**: String representing zero energy control. (optional)
    - **C11/C12/.../C66**: Numeric values representing material constants. (Float64 or Int64, optional)
- **Thermal Models**: Dictionary of thermal models. (optional)
  - **_Own_Name_**: List of thermal model configurations.
    - **Thermal Model**: String representing the thermal model.
    - **Type**: String representing the type of thermal model. (optional)
    - **Heat Transfer Coefficient/Environmental Temperature/Thermal Conductivity/Thermal Expansion Coefficient/Thermal Conductivity Print Bed/Print Bed Temperature**: Numeric values representing thermal parameters. (Float64 or Int64, optional)
- **Additive Models**: Dictionary of additive models. (optional)
  - **_Own_Name_**: List of additive model configurations.
    - **Additive Model**: String representing the additive model.
    - **Print Temperature**: Numeric value representing the print temperature. (Float64 or Int64, optional)
- **Pre Calculation**: Dictionary of pre-calculation parameters.
  - **Bond Associated Deformation Gradient/Bond Associated Shape Tensor/Deformation Gradient/Deformed Bond Geometry/Shape Tensor**: Boolean values indicating whether to calculate the respective parameter. (optional)

## `Solver`

- **Solve For Displacement/Material Models/Damage Models/Thermal Models/Additive Models**: Boolean values indicating whether to solve for the respective components. (optional)
- **Maximum Damage**: Numeric value representing the maximum damage. (Float64, optional)
- **Final Time/Initial Time**: Numeric values representing the final and initial time. (Float64 or Int64, required)
- **Numerical Damping**: Numeric value representing numerical damping. (Float64 or Int64, optional)
- **Verlet**: Dictionary of Verlet solver parameters.
  - **Safety Factor/Fixed dt/Number of Steps**: Numeric values representing Verlet solver parameters. (Float64 or Int64, optional)
- **External**: Dictionary of external solver parameters.
  - **Number of Steps**: Numeric value representing the number of steps. (Int64, optional)

