# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    InputDeck

Typed declarations of PeriLab's fixed input sections and `read_input`, which
turns the `PeriLab:` part of a YAML deck into a validated `PeriLabInput`.
Model parameters (`Models`) are declared by the model modules (phase 3).
"""
module InputDeck

using ..ParameterSpec
using ..ParameterSpec: @params, req, opt, ParseContext, add_error!, join_path,
                       parse_section
import ..ParameterSpec: check!
using ..PeriLabExceptions: @abort

include("discretization.jl")
include("blocks.jl")
include("solver.jl")
include("outputs.jl")
include("conditions.jl")
include("contact.jl")
include("input.jl")
include("generators.jl")

export read_input, PeriLabInput, SolverParams, SolverOptions, ModelReductionParams,
       active_options, solver_name, start_time, end_time, model_options, solver_steps,
       solver_step,
       BlockParams, DiscretizationParams, GcodeParams, BondFilterParams,
       SurfaceExtrusionParams, ExternalTopologyParams, block_by_id, block_angles,
       block_names_and_ids, mesh_scaling, gcode_block_ids, OutputParams, ComputeClassParams,
       check_for_duplicates, output_filenames, output_frequencies, output_fieldnames,
       compute_names, active_computes,
       BoundaryConditionParams, bc_node_set_names, bc_step_ids,
       ContactInput, ContactModelParams, ContactGroupParams, ContactGlobalsParams,
       contact_blocks, contact_search_frequency, FEMParams, FEMCouplingParams, fem_degree,
       SurfaceCorrectionParams

end
