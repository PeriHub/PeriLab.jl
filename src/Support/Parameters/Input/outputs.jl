# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params struct OutputParams
    output_filename::String = req("Output Filename")
    output_file_type::String = opt("Output File Type"; default = "Exodus",
                                   allowed = ["Exodus", "CSV"])
    output_frequency::Union{Nothing,Int64,String} = opt("Output Frequency"; default = nothing,
                                                        description = "Every n-th step; per step as \"n1 n2\"")
    number_of_output_steps::Union{Nothing,Int64,String} = opt("Number of Output Steps";
                                                              default = nothing)
    output_variables::Dict{String,Bool} = req("Output Variables")
    flush_file::Bool = opt("Flush File"; default = true)
    write_after_damage::Bool = opt("Write After Damage"; default = false)
    start_time::Float64 = opt("Start Time"; default = 0.0, quantity = :time)
    end_time::Float64 = opt("End Time"; default = Inf, quantity = :time)
    bond_export::Bool = opt("Bond Export"; default = false)
    bond_blocks::Union{Nothing,Int64,String} = opt("Bond Blocks"; default = nothing,
                                                   description = "Block id or list, e.g. \"1 3\"")
end

function check!(p::OutputParams, path::String, ctx::ParseContext)
    if p.output_frequency === nothing && p.number_of_output_steps === nothing
        add_error!(ctx, path, "\"Output Frequency\" or \"Number of Output Steps\" is required")
    end
    return nothing
end
