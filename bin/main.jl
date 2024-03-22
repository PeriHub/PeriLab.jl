# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using PeriLab
using ZipArchives: ZipWriter, zip_newfile

function zip_folder(source_folder::AbstractString, zip_filename::AbstractString)
    files = readdir(source_folder, join=true)

    ZipWriter(zip_filename) do w
        for (root, dirs, files) in walkdir(source_folder)
            for file in files
                filepath = joinpath(root, file)
                f = open(filepath, "r")
                content = read(f, String)
                close(f)
                zip_newfile(w, basename(filepath))
                write(w, content)
            end
        end
    end
end

filename = get(ENV, "filename", "examples/DCB/DCBmodel.yaml")
output_dir = "./results"
dryrun = get(ENV, "dryrun", "false") == "true"
verbose = get(ENV, "verbose", "false") == "true"
debug = get(ENV, "debug", "false") == "true"
silent = get(ENV, "silent", "false") == "true"

PeriLab.main(filename; output_dir="./results", dryrun=dryrun, verbose=verbose, debug=debug, silent=silent)

zip_filename = "results.zip"
zip_folder(output_dir, zip_filename)

ENV["RESULTS_FILE"] = zip_filename