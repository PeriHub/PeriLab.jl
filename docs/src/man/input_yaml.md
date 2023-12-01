# PeriLab Configuration File

The PeriLab configuration file is a YAML file used to specify the parameters for running simulations in the PeriLab software. Below is the expected structure of the configuration file:

```yaml
PeriLab:
  - Blocks:
    - Any:
        - Block Names: [String, true]
        - Density: [Union{Float64,Int64}, true]
        - Horizon: [Union{Float64,Int64}, true]
        - Specific Heat Capacity: [Union{Float64,Int64}, false]
        - Material Model: [String, false]
        - Damage Model: [String, false]
        - Thermal Model: [String, false]
        - Additive Model: [String, false]

  - Boundary Conditions:
    - Any:
        - Coordinate: [String, false]
        - Node Set: [String, true]
        - Type: [String, true]
        - Value: [Union{Float64,Int64,String}, true]

  - Compute Class Parameters:
    - Any:
        - Block: [String, true]
        - Calculation Type: [String, true]
        - Compute Class: [String, true]
        - Variable: [String, true]

  - Discretization:
    - Input Mesh File: [String, true]
    - Node Sets:
      - Any: [String, true]
    - Type: [String, true]
    - Bond Filters:
      - Any:
          - Type: [String, true]
          - Normal X: [Union{Float64,Int64}, true]
          - Normal Y: [Union{Float64,Int64}, true]
          - Normal Z: [Union{Float64,Int64}, true]
          - Lower Left Corner X: [Union{Float64,Int64}, false]
          - Lower Left Corner Y: [Union{Float64,Int64}, false]
          - Lower Left Corner Z: [Union{Float64,Int64}, false]
          - Bottom Unit Vector X: [Union{Float64,Int64}, false]
          - Bottom Unit Vector Y: [Union{Float64,Int64}, false]
          - Bottom Unit Vector Z: [Union{Float64,Int64}, false]
          - Center X: [Union{Float64,Int64}, false]
          - Center Y: [Union{Float64,Int64}, false]
          - Center Z: [Union{Float64,Int64}, false]
          - Radius: [Union{Float64,Int64}, false]
          - Bottom Length: [Union{Float64,Int64}, false]
          - Side Length: [Union{Float64,Int64}, false]

  - Outputs:
    - Any:
        - Flush File: [Bool, false]
        - Output Frequency: [Int64, false]
        - Number of Output Steps: [Int64, false]
        - Output File Type: [String, false]
        - Output Filename: [String, true]
        - Write After Damage: [Bool, false]
        - Output Variables:
          - Any: [Bool, true]

  - Physics:
    - Damage Models:
      - Any:
          - Critical Value: [Union{Float64,Int64}, true]
          - Damage Model: [String, true]
          - Interblock Damage:
            - Any: [Union{Float64,Int64}, true]

    - Material Models:
      - Any:
          - Material Model: [String, true]
          - Symmetry: [String, false]
          - Poisson's Ratio: [Union{Float64,Int64}, false]
          - Young's Modulus: [Union{Float64,Int64}, false]
          - Bulk Modulus: [Union{Float64,Int64}, false]
          - Shear Modulus: [Union{Float64,Int64}, false]
          - Zero Energy Control: [String, false]
          - C11: [Union{Float64,Int64}, false]
          - C12: [Union{Float64,Int64}, false]
          - ... (other material properties)

    - Thermal Models:
      - Any:
          - Thermal Model: [String, true]
          - Type: [String, false]
          - Lambda: [Union{Float64,Int64}, false]
          - Environmental Temperature: [Union{Float64,Int64}, false]
          - Kappa: [Union{Float64,Int64}, false]
          - Heat expansion: [Union{Float64,Int64}, false]

    - Additive Models:
      - Any:
          - Additive Model: [String, true]
          - Print Temperature: [Union{Float64,Int64}, false]

    - Pre Calculation:
      - Bond Associated Deformation Gradient: [Bool, false]
      - Bond Associated Shape Tensor: [Bool, false]
      - Deformation Gradient: [Bool, false]
      - Deformed Bond Geometry: [Bool, false]
      - Shape Tensor: [Bool, false]

  - Solver:
    - Solve For Displacement: [Bool, false]
    - Material Models: [Bool, false]
    - Damage Models: [Bool, false]
    - Thermal Models: [Bool, false]
    - Additive Models: [Bool, false]
    - Final Time: [Union{Float64,Int64}, true]
    - Initial Time: [Union{Float64,Int64}, true]
    - Verbose: [Bool, false]
    - Verlet:
        - Numerical Damping: [Union{Float64,Int64}, false]
        - Safety Factor: [Union{Float64,Int64}, false]
        - Fixed dt: [Union{Float64,Int64}, false]
