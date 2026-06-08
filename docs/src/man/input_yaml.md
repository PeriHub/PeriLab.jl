# PeriLab Configuration File

The PeriLab configuration file is a YAML file used to specify the parameters for running simulations in the PeriLab software. This structure represents the global configuration for PeriLab.

The information is stored in the [params dictionary](@ref "Parameters")

## `PeriLab`

- **Blocks**: List of block configurations.
  - *__Own_Name__*: Block-specific parameters.
    - **Block ID**: Block ID corresponding to the input mesh (Int)
    - **Step ID**: Step ID for multisolver simulations (Int, optional)
    - **Density**: Numeric value representing density. (Float or Int)
    - **Horizon**: Numeric value representing horizon. (Float or Int)
    - **Specific Heat Capacity**: Numeric value representing specific heat capacity. (Float or Int, optional)
    - **Material Model**: String representing the material model. (optional)
    - **Damage Model**: String representing the damage model. (optional)
    - **Thermal Model**: String representing the thermal model. (optional)
    - **Additive Model**: String representing the additive model. (optional)
    - **Pre Calculation Model**: String representing the pre-calculation model. (optional)
    - **Angle X**: Numeric value representing angle in X direction. (Float or Int, optional)
    - **Angle Y**: Numeric value representing angle in Y direction. (Float or Int, optional)
    - **Angle Z**: Numeric value representing angle in Z direction. (Float or Int, optional)

## `FEM (optional)`

- **Element Type**: String representing the type of finite element.
- **Degree**: Numeric value representing the degree. (String or Int)
- **Material Model**: String representing the material model.
- **Coupling**: Coupling parameters. (optional)
  - **Coupling Type**: String representing the type of coupling.
  - **PD Weight**: Peridynamic Weight, default is 0.5 (Float, optional)
  - **Kappa**: Coupling stiffness, default is 1.0 (Float, optional)

## `Boundary Conditions (optional)`

- *__Own_Name__*: List of boundary condition configurations.
  - **Coordinate**: String representing the coordinate.
  - **Node Set**: String representing the node set.
  - **Variable**: String representing the variable.
  - **Type**: String representing the type of boundary condition.
  - **Value**: Numeric value, string, or a combination representing the boundary condition value. (Float, Int, String)

## `Compute Class Parameters (optional)`

- *__Own_Name__*: List of compute class parameters.
  - **Block**: String representing the block.
  - **Node Set**: String representing the node set.
  - **Calculation Type**: String representing the calculation type.
  - **Compute Class**: String representing the compute class.
  - **Variable**: String representing the variable.
  - **Equation**: Additional conversion of numeric value f.e. "-x" (optional)
  - **X**: X Value for Nearest_Point_Data. (optional)
  - **Y**: Y Value for Nearest_Point_Data. (optional)
  - **Z**: Z Value for Nearest_Point_Data. (optional)

## `Discretization`

- **Input Mesh File**: String representing the input mesh file.
- **Input External Topology**: External topology parameters. (optional)
  - **File**: String representing the external topology file.
  - **Add Neighbor Search**: Boolean indicating whether to add neighbor search. (optional)
- **Node Sets**: Dictionary of node sets. (optional)
  - *__Own_Name__*: Numeric value or string representing the node set. (Int or String)
- **Type**: String representing the type of discretization.
- **Distribution Type**: String representing the distribution type. (optional)
- **Bond Filters**: List of bond filter configurations. (optional)
  - *__Own_Name__*: Bond filter parameters.
    - **Type**: String representing the bond filter type.
    - **Normal X/Y/Z**: Numeric values representing normal components. (Float or Int)
    - **Lower Left Corner X/Y/Z**: Numeric values representing lower-left corner components. (Float or Int)
    - **Bottom Unit Vector X/Y/Z**: Numeric values representing bottom unit vector components. (Float or Int)
    - **Center X/Y/Z**: Numeric values representing center components. (Float or Int)
    - **Radius**: Numeric value representing the radius. (Float or Int)
    - **Bottom Length**: Numeric value representing the bottom length. (Float or Int)
    - **Side Length**: Numeric value representing the side length. (Float or Int)
    - **Allow Contact**: Boolean indicating whether to allow contact. (optional)
- **Surface Extrusion**: Surface extrusion parameters. (optional)
  - **Direction**: String representing surface extrusion direction (e.g., 'X', or 'Y').
  - **Step_X/Y/Z**: Numeric values representing step components for surface extrusion direction (Float or Int).

## `Outputs (optional)`

- *__Own_Name__*: List of output configurations.
  - **Flush File**: Boolean indicating whether to flush the output file. (optional)
  - **Output Frequency**: Integer representing the output frequency.
  - **Number of Output Steps**: Integer representing the number of output steps.
  - **Output File Type**: String representing the output file type. (optional)
  - **Output Filename**: String representing the output filename.
  - **Write After Damage**: Boolean indicating whether to write after damage. (optional)
  - **Output Variables**: Dictionary of output variables.
    - *__Own_Name__*: Boolean indicating whether to output the variable.

## `Models (optional)`

- **Damage Models**: Dictionary of damage models. (optional)
  - *__Own_Name__*: List of damage model configurations.
    - **Critical Value**: Numeric value representing the critical value. (Float or Int)
    - **Damage Model**: String representing the damage model.
    - **Interblock Damage**: Dictionary of interblock damage parameters.
      - **Interblock Critical Value**: Numeric value representing interblock damage. (Float or Int, required)
    - **Anisotropic Damage**: Dictionary of anisotropic damage parameters.
      - **Critical Value X/Y**: Numeric values representing critical values in X and Y directions. (Float or Int, required)
- **Material Models**: Dictionary of [material models](models/materials.md). (optional)
  - *__Own_Name__*: List of material model configurations.
    - **Material Model**: String representing the material model.
    - **Symmetry**: String representing the symmetry. (optional)
    - **Poisson's Ratio/Young's Modulus/Bulk Modulus/Shear Modulus**: Numeric values representing material properties. (Float or Int, optional)
    - **Yield Stress**: Numeric value representing the yield stress. (Float or Int, optional)
    - **Zero Energy Control**: String representing zero energy control. (optional)
    - **C11/C12/.../C66**: Numeric values representing material constants. (Float or Int, optional)
- **Thermal Models**: Dictionary of thermal models. (optional)
  - *__Own_Name__*: List of thermal model configurations.
    - **Thermal Model**: String representing the thermal model.
    - **Type**: String representing the type of thermal model. (optional)
    - **Heat Transfer Coefficient/Environmental Temperature/Thermal Conductivity/Thermal Expansion Coefficient/Thermal Conductivity Print Bed/Print Bed Temperature**: Numeric values representing thermal parameters. (Float or Int, optional)
- **Additive Models**: Dictionary of additive models. (optional)
  - *__Own_Name__*: List of additive model configurations.
    - **Additive Model**: String representing the additive model.
    - **Print Temperature**: Numeric value representing the print temperature. (Float or Int, optional)
- **Pre Calculation Models**: Dictionary of pre-calculation models. (optional)
  - *__Own_Name__* : List of pre calculation model configurations.. (optional).
    - **Bond Associated Deformation Gradient/Deformation Gradient/Deformed Bond Geometry/Shape Tensor**: Boolean values indicating whether to calculate the respective parameter. (optional)
- **Pre Calculation Global**: Dictionary of pre-calculation parameters.
  - **Bond Associated Deformation Gradient/Deformation Gradient/Deformed Bond Geometry/Shape Tensor**: Boolean values indicating whether to calculate the respective parameter. (optional)

## `Solver`

- **Additive Models/Material Models/Damage Models/Thermal Models/Pre Calculation Models**: Boolean values indicating whether to solve for the respective components. (optional)
- **Calculate Cauchy**: Boolean values indicating whether to calculate Cauchy stresses for deformation gradient and shape tensor calculations (optional).
- **Calculate von Mises**: Boolean values indicating whether to calculate von Mises stresses for deformation gradient and shape tensor calculations (optional).
- **Calculate Strain** : Boolean values indicating whether to calculate strain for deformation gradient and shape tensor calculations (optional).
- **Maximum Damage**: Numeric value representing the maximum damage. (Float, optional)
- **Final Time/Initial Time**: Numeric values representing the final and initial time. (Float or Int, required)
- **Fixed dt/Number of Steps**: Numeric values defining the step width. (Float or Int, optional)
- **Verlet**: Dictionary of Verlet solver parameters.
  - **Safety Factor**: Numeric values representing Verlet solver parameters. (Float or Int, optional)
  - **Numerical Damping**: Numeric value representing numerical damping. (Float or Int, optional)
- **Static**: Dictionary of Static solver parameters.
  - **Number of Steps**: Numeric value representing the number of steps. (Int, optional)
  - **NLsolve,Solution tolerance  | Float, Residual tolerance, Maximum number of iterations, Show solver iteration, Residual scaling, Solver Type, m, Linear Start Value** Numeric values representing Static solver parameters. (Float, Bool, String or Int, optional)

## `Surface Correction (optional)`
- **Type**: String defining what model is used, currently only Volume Correction is included
- **Update**: Bool defining if the surface is updated during the solving process, taking changes as fracture or additive processes into account

##  `Contact (optional)`
tbd
