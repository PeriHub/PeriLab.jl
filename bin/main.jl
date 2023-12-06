# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using PeriLab
using ZipFile

function zip_folder(source_folder::AbstractString, zip_filename::AbstractString)
    files = readdir(source_folder, join=true)

    # Create a Zip file
    zip = ZipFile.ZipWriter(zip_filename)

    for file in files
        # Add each file to the Zip file
        addfile(zip, file, file)
    end

    # Close the Zip file
    close(zip)
end

filename = get(ENV, "filename", "examples/DCB/DCBmodel.yaml")
output_dir = "./results"
dryrun = get(ENV, "dryrun", "false") == "true"
verbose = get(ENV, "verbose", "false") == "true"
debug = get(ENV, "debug", "false") == "true"
silent = get(ENV, "silent", "false") == "true"

PeriLab.main(filename, "./results", dryrun, verbose, debug, silent)

zip_filename = "results.zip"
zip_folder(output_dir, zip_filename)

ENV["RESULTS_FILE"] = zip_filename