# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Write_CSV_Results
export create_result_file
export write_global_results_in_csv

function create_result_file(filename::String, outputs::Dict)

    if isfile(filename)
        rm(filename)
    end
    @info "Create output " * filename
    csv_file = open(filename, "w")

    header = ""
    for key in keys(outputs["Fields"])
        header = string(header, key, ",")
    end
    write(csv_file, header * "\n")

    return Dict("filename" => filename, "file" => csv_file)
end

function write_global_results_in_csv(csv_file::Dict, global_values)
    value_string = ""
    for value in global_values
        value_string = string(value_string, value, ",")
    end
    value_string = value_string * "\n"
    write(csv_file["file"], value_string)

end
end