# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Write_CSV_Results
export create_result_file
export write_global_results_in_csv

function create_result_file(filenames, computes)
    csv_files = []
    for (id, filename) in enumerate(filenames)

        if !haskey(computes, id)
            continue
        end
        if ".e" == filename[end-1:end]
            filename = filename[1:end-2]
        end
        filename = filename * "_globals.csv"
        if isfile(filename)
            rm(filename)
        end
        @info "Create output " * filename
        csv_file = open(filename, "w")

        header = ""
        for key in keys(computes[id])
            if computes[id][key]["CSV Export"]
                for fieldname in keys(computes[id][key]["Mapping"])
                    header = string(header, fieldname, ",")
                end
            end
        end
        write(csv_file, header * "\n")

        push!(csv_files, csv_file)
    end
    return csv_files
end

function write_global_results_in_csv(csv_file, values)

    # files = csv_files
    value_string = ""
    for value in values
        value_string = string(value_string, value, ",")
    end
    value_string = value_string * "\n"
    write(csv_file, value_string)
end
end