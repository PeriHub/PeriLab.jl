# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using CSV
using OrderedCollections: OrderedDict
export create_result_file
export write_global_results_in_csv

"""
    create_result_file(filename::String, outputs::Dict)

Creates a csv file for the results

# Arguments
- `filename::String`: The name of the file to create
- `outputs::Dict`: The outputs dictionary
# Returns
- `Dict`: The result file
"""
function create_result_file(filename::String, outputs::Dict)
    if isfile(filename)
        rm(filename)
    end
    @info "Create output " * filename
    csv_file = open(filename, "w")

    header = "Time,"
    for key in keys(sort!(OrderedDict(outputs["Fields"])))
        header = string(header, key, ",")
    end
    write(csv_file, header * "\n")

    return Dict("filename" => filename, "file" => csv_file, "type" => "CSV")
end

"""
    write_global_results_in_csv(csv_file::IOStream, time::Float64, global_values)

Writes the global results to the csv file

# Arguments
- `csv_file::IOStream`: The csv file
- `global_values`: The global values
"""
function write_global_results_in_csv(csv_file::IOStream, time::Float64, global_values)
    value_string = string(time) * ","
    for value in global_values
        value_string = string(value_string, value, ",")
    end
    value_string = value_string * "\n"
    write(csv_file, value_string)
end
