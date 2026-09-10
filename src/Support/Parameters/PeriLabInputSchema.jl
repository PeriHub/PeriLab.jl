# SPDX-License-Identifier: BSD-3-Clause
#
# The single source of truth for PeriLab's input YAML shape.
# Replaces `expected_structure` in parameter_handling.jl.
#
# Usage:
#   errors = Schema.validate_params(PERILAB_SCHEMA, params["PeriLab"])
#   isempty(errors) || foreach(println, errors)
#
#   json_doc = Schema.to_json_schema(PERILAB_SCHEMA, Val(:document);
#                                     id = "https://perilab.dlr.de/schema/input.json")
#   open("perilab.schema.json", "w") do io
#       JSON3.write(io, json_doc)
#   end

include("./Schema.jl")
using .Schema

const NUM = [Float64, Int64]

# --- Blocks -----------------------------------------------------------------

const BLOCK_ENTRY = SObject(Dict(
    "Block ID" => SField(Int64; required = true),
    "Step ID" => SField([Int64, String]),
    "Density" => SField(NUM; required = true, min = 0),
    "Horizon" => SField(NUM; required = true, min = 0),
    "Specific Heat Capacity" => SField(NUM, min = 0),
    "Material Model" => SField(String,
                               description = "Name referencing a block under Models.Material Models"),
    "Damage Model" => SField(String),
    "Thermal Model" => SField(String),
    "Additive Model" => SField(String),
    "Pre Calculation Model" => SField(String),
    "Degradation_template Model" => SField(String),
    "Angle X" => SField(NUM), "Angle Y" => SField(NUM), "Angle Z" => SField(NUM),
))
const BLOCKS = SAny(BLOCK_ENTRY; required = true, min_entries = 1,
                    description = "Keyed by user-chosen block name")

# --- FEM ---------------------------------------------------------------------

const FEM_COUPLING = SObject(Dict(
    "Coupling Type" => SField(String; required = true),
    "PD Weight" => SField(NUM),
    "Kappa" => SField(NUM),
))
const FEM = SObject(Dict(
    "Element Type" => SField(String; required = true),
    "Degree" => SField([String, Int64]; required = true),
    "Material Model" => SField(String; required = true),
    "Coupling" => FEM_COUPLING,
))

# --- Boundary Conditions -------------------------------------------------------

const BC_ENTRY = SObject(Dict(
    "Coordinate" => SField(String),
    "Node Set" => SField(String; required = true),
    "Variable" => SField(String; required = true),
    "Type" => SField(String, enum = ["Dirichlet", "Neumann", "Initial"]),
    "Value" => SField([Float64, Int64, String]; required = true),
))
const BOUNDARY_CONDITIONS = SAny(BC_ENTRY)

# --- Compute Class Parameters --------------------------------------------------

const COMPUTE_ENTRY = SObject(Dict(
    "Block" => SField(String),
    "Node Set" => SField(String),
    "Calculation Type" => SField(String),
    "Compute Class" => SField(String; required = true),
    "X" => SField(NUM), "Y" => SField(NUM), "Z" => SField(NUM),
    "Variable" => SField(String; required = true),
    "Equation" => SField(String),
))
const COMPUTE_CLASS_PARAMETERS = SAny(COMPUTE_ENTRY)

# --- Discretization / Mesh ------------------------------------------------------

const BOND_FILTER_ENTRY = SObject(Dict(
    "Type" => SField(String; required = true, enum = ["Rectangular_Plane", "Disk", "Wedge"]),
    "Normal X" => SField(NUM; required = true),
    "Normal Y" => SField(NUM; required = true),
    "Normal Z" => SField(NUM),
    "Lower Left Corner X" => SField(NUM), "Lower Left Corner Y" => SField(NUM),
    "Lower Left Corner Z" => SField(NUM),
    "Bottom Unit Vector X" => SField(NUM), "Bottom Unit Vector Y" => SField(NUM),
    "Bottom Unit Vector Z" => SField(NUM),
    "Center X" => SField(NUM), "Center Y" => SField(NUM), "Center Z" => SField(NUM),
    "Radius" => SField(NUM, min = 0),
    "Bottom Length" => SField(NUM, min = 0), "Side Length" => SField(NUM, min = 0),
    "Allow Contact" => SField(Bool),
))
const GCODE = SObject(Dict(
    "Overwrite Mesh" => SField(Bool; required = true),
    "Sampling" => SField(NUM; required = true),
    "Width" => SField(NUM; required = true),
    "Height" => SField(NUM; required = true),
    "Scale" => SField(NUM),
    "Start Command" => SField(String), "Stop Command" => SField(String),
    "End Command" => SField(String),
))
const SURFACE_EXTRUSION = SObject(Dict(
    "Direction" => SField(String; required = true, enum = ["X", "Y", "Z"]),
    "Step_X" => SField(NUM; required = true), "Step_Y" => SField(NUM; required = true),
    "Step_Z" => SField(NUM; required = true),
    "Number" => SField(NUM; required = true),
))
const DISCRETIZATION = SObject(Dict(
    "Input Mesh File" => SField(String; required = true),
    "Input External Topology" => SObject(Dict(
        "File" => SField(String; required = true),
        "Add Neighbor Search" => SField(Bool),
    )),
    "Node Sets" => SAny(SField([Int64, String])),
    "Type" => SField(String; required = true, enum = ["Text File", "Exodus"]),
    "Distribution Type" => SField(String),
    "Surface Extrusion" => SURFACE_EXTRUSION,
    "Bond Filters" => SAny(BOND_FILTER_ENTRY),
    "Horizon Mesh Scaling X" => SField(NUM), "Horizon Mesh Scaling Y" => SField(NUM),
    "Horizon Mesh Scaling Z" => SField(NUM),
    "Gcode" => GCODE,
); required = true)

# --- Outputs ---------------------------------------------------------------------

const OUTPUT_ENTRY = SObject(Dict(
    "Flush File" => SField(Bool),
    "Output Frequency" => SField([String, Int64]),
    "Number of Output Steps" => SField([String, Int64]),
    "Output File Type" => SField(String, enum = ["Exodus", "CSV"], default = "Exodus"),
    "Output Filename" => SField(String; required = true),
    "Write After Damage" => SField(Bool),
    "Start Time" => SField(NUM), "End Time" => SField(NUM),
    "Output Variables" => SAny(SField(Bool); required = true),
))
const OUTPUTS = SAny(OUTPUT_ENTRY)

# --- Models --------------------------------------------------------------------

const FLAW_FUNCTION = SObject(Dict(
    "Active" => SField(Bool; required = true),
    "Function" => SField(String; required = true),
    "Flaw Size" => SField(NUM), "Flaw Magnitude" => SField(NUM),
    "Flaw Location X" => SField(NUM), "Flaw Location Y" => SField(NUM),
    "Flaw Location Z" => SField(NUM),
))
const DAMAGE_MODEL_ENTRY = SObject(Dict(
    "Critical Value" => SField([Float64, Int64, String]; required = true),
    "Damage Model" => SField(String; required = true),
    "Interblock Damage" => SAny(SField(NUM)),
    "Only Tension" => SField(Bool),
    "Thickness" => SField(NUM),
    "Anisotropic Damage" => SObject(Dict(
        "Critical Value X" => SField(NUM; required = true),
        "Critical Value Y" => SField(NUM; required = true),
    )),
    "Flaw Function" => FLAW_FUNCTION,
))
const CIJ_KEYS = ["C$i$j" for (i, j) in
                  [(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(2,2),(2,3),(2,4),(2,5),(2,6),
                   (3,3),(3,4),(3,5),(3,6),(4,4),(4,5),(4,6),(5,5),(5,6),(6,6)]]

const MATERIAL_MODEL_ENTRY_BASE = SObject(Dict(
    "Material Model" => SField(String; required = true),
    "Symmetry" => SField(String, enum = ["isotropic", "anisotropic", "orthotropic"]),
    "Bond Associated" => SField(Bool),
    "Poisson's Ratio" => SField(NUM, min = -1, max = 0.5),
    "Poisson's Ratio XY" => SField([Float64, Int64, String]),
    "Poisson's Ratio YZ" => SField([Float64, Int64, String]),
    "Poisson's Ratio XZ" => SField([Float64, Int64, String]),
    "Young's Modulus" => SField(NUM, min = 0),
    "Young's Modulus X" => SField([Float64, Int64, String]),
    "Young's Modulus Y" => SField([Float64, Int64, String]),
    "Young's Modulus Z" => SField([Float64, Int64, String]),
    "Bulk Modulus" => SField(NUM, min = 0), "Shear Modulus" => SField(NUM, min = 0),
    "Shear Modulus XY" => SField([Float64, Int64, String]),
    "Shear Modulus YZ" => SField([Float64, Int64, String]),
    "Shear Modulus XZ" => SField([Float64, Int64, String]),
    "Yield Stress" => SField([Float64, Int64, String]),
    "Zero Energy Control" => SField(String),
    Dict(k => SField(NUM) for k in CIJ_KEYS)...,  # Cij stiffness matrix entries
    "File" => SField(String),
    "Number of State Variables" => SField(Int64, min = 0),
    "Number of Properties" => SField(Int64, min = 0),
    "Predefined Field Names" => SField(String),
    "State Factor ID" => SField(Int64),
    "Accuracy Order" => SField(Int64),
    "Flaw Function" => FLAW_FUNCTION,
))

# Symmetry == "anisotropic" unconditionally requires the full 21-entry
# stiffness matrix (mirrors the iID/jID 1:6 loop that currently does this
# with @abort per missing entry).
const MATERIAL_MODEL_ENTRY = SWith(MATERIAL_MODEL_ENTRY_BASE, [
    required_if("Symmetry", "anisotropic", CIJ_KEYS),
])

const THERMAL_MODEL_ENTRY = SObject(Dict(
    "Thermal Model" => SField(String; required = true),
    "Type" => SField(String),
    "Heat Transfer Coefficient" => SField(NUM),
    "Environmental Temperature" => SField([Float64, Int64, String]),
    "Allow Surface Change" => SField(Bool),
    "Thermal Conductivity" => SField(NUM, min = 0),
    "Thermal Expansion Coefficient" => SField(NUM),
    "Reference Temperature" => SField(NUM),
    "Thermal Conductivity Print Bed" => SField(NUM),
    "Print Bed Temperature" => SField(NUM),
    "Print Bed Z Coordinate" => SField(NUM),
    "File" => SField(String),
    "Number of State Variables" => SField(Int64, min = 0),
    "Predefined Field Names" => SField(String),
))
const ADDITIVE_MODEL_ENTRY = SObject(Dict(
    "Additive Model" => SField(String; required = true),
    "Print Temperature" => SField(NUM),
))
const DEGRADATION_MODEL_ENTRY = SObject(Dict(
    "Degradation Model" => SField(String; required = true),
    "Decomposition Temperature" => SField(NUM),
))
const PRE_CALC_GLOBAL = SObject(Dict(
    "Bond Associated Deformation Gradient" => SField(Bool),
    "Bond Associated Correspondence" => SField(Bool),
    "Deformation Gradient" => SField(Bool),
    "Deformed Bond Geometry" => SField(Bool),
    "Shape Tensor" => SField(Bool),
))
const PRE_CALC_MODEL_ENTRY = SObject(Dict(
    "Bond Associated Correspondence" => SField(Bool),
    "Deformation Gradient" => SField(Bool),
    "Deformed Bond Geometry" => SField(Bool),
    "Shape Tensor" => SField(Bool),
))
const MODELS = SObject(Dict(
    "Damage Models" => SAny(DAMAGE_MODEL_ENTRY),
    "Material Models" => SAny(MATERIAL_MODEL_ENTRY),
    "Thermal Models" => SAny(THERMAL_MODEL_ENTRY),
    "Additive Models" => SAny(ADDITIVE_MODEL_ENTRY),
    "Degradation Models" => SAny(DEGRADATION_MODEL_ENTRY),
    "Pre Calculation Global" => PRE_CALC_GLOBAL,
    "Pre Calculation Models" => SAny(PRE_CALC_MODEL_ENTRY),
))

# --- Contact ---------------------------------------------------------------------

const CONTACT_GROUP_ENTRY = SObject(Dict(
    "Master Block ID" => SField(Int64; required = true),
    "Slave Block ID" => SField(Int64; required = true),
    "Search Radius" => SField(NUM; required = true, min = 0),
    "Global Search Frequency" => SField(Int64, min = 1),
    "Maximum Contact Pairs" => SField(Int64, min = 1),
))
const CONTACT_ENTRY = SObject(Dict(
    "Type" => SField(String; required = true),
    "Contact Radius" => SField(NUM; required = true, min = 0),
    "Contact Stiffness" => SField(NUM; required = true, min = 0),
    "Contact Groups" => SAny(CONTACT_GROUP_ENTRY; required = true),
))
const CONTACT = SAny(CONTACT_ENTRY)

# --- Solver ------------------------------------------------------------------------
# NOTE: this is the clearest win from SOneOf. Previously "which solver is
# active" was implicit — whichever of Verlet/Static/... key happened to be
# present — and enforced only by the imperative if/elseif chain in
# get_solver_name(). Now it's declared in the schema itself.

const VERLET = SObject(Dict(
    "Safety Factor" => SField(NUM, min = 0, max = 1),
    "Fixed dt" => SField(NUM, min = 0),
    "Numerical Damping" => SField(NUM, min = 0),
); required = true)
const LINEAR_STATIC_MATRIX_BASED = SObject(Dict(
    "Safety Factor" => SField(NUM),
    "Matrix Update" => SField(Bool),
); required = true)
const NEWMARK = SObject(Dict(
    "Safety Factor" => SField(NUM),
    "Matrix Update" => SField(Bool),
); required = true)
const STATIC = SObject(Dict(
    "Maximum number of iterations" => SField(Int64, min = 1),
    "NLSolve" => SField(Bool),
    "Show solver iteration" => SField(Bool),
    "Solver Type" => SField(String),
    "Residual scaling" => SField(NUM),
    "Solution tolerance" => SField(NUM, min = 0),
    "Residual tolerance" => SField(NUM, min = 0),
    "Linear Start Value" => SField(String),
    "m" => SField(Int64),
); required = true)

const SOLVER_KIND = SOneOf(Dict(
    "Verlet" => VERLET,
    "Static" => STATIC,
    "Linear Static Matrix Based" => LINEAR_STATIC_MATRIX_BASED,
    "Newmark" => NEWMARK,
); required = true, description = "Exactly one solver algorithm must be selected")

# Common solver-level flags shared by Solver and each Multistep Solver step.
_solver_common_fields() = Dict(
    "Additive Models" => SField(Bool),
    "Degradation Models" => SField(Bool),
    "Damage Models" => SField(Bool),
    "Material Models" => SField(Bool, default = true),
    "Thermal Models" => SField(Bool),
    "Pre Calculation Models" => SField(Bool, default = true),
    "Calculate Cauchy" => SField(Bool),
    "Calculate von Mises stress" => SField(Bool),
    "Calculate Strain" => SField(Bool),
    "Maximum Damage" => SField(NUM, min = 0, max = 1),
    "Verlet" => VERLET, "Static" => STATIC, "Newmark" => NEWMARK,
    "Linear Static Matrix Based" => LINEAR_STATIC_MATRIX_BASED,
)

const SOLVER = SObject(merge(_solver_common_fields(), Dict(
    "Final Time" => SField(NUM; required = true),
    "Initial Time" => SField(NUM; required = true),
    "Number of Steps" => SField(Int64, min = 1),
)); required = true)

const MULTISTEP_ENTRY = SObject(merge(_solver_common_fields(), Dict(
    "Step ID" => SField([Int64, String]),
    "Final Time" => SField(NUM),
    "Initial Time" => SField(NUM),
    "Additional Time" => SField(NUM),
    "Number of Steps" => SField(Int64, min = 1),
)))
const MULTISTEP_SOLVER = SAny(MULTISTEP_ENTRY)

# --- Root ----------------------------------------------------------------------

const PERILAB_SCHEMA = SObject(Dict(
    "Blocks" => BLOCKS,
    "FEM" => FEM,
    "Boundary Conditions" => BOUNDARY_CONDITIONS,
    "Compute Class Parameters" => COMPUTE_CLASS_PARAMETERS,
    "Discretization" => DISCRETIZATION,
    "Outputs" => OUTPUTS,
    "Models" => MODELS,
    "Contact" => CONTACT,
    "Surface Correction" => SObject(Dict(
        "Type" => SField(String; required = true),
        "Update" => SField(Bool),
    )),
    "Solver" => SOLVER,
    "Multistep Solver" => MULTISTEP_SOLVER,
); required = true, description = "PeriLab top-level input")
