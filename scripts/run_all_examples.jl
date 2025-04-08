# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using PeriLab

example_path = "../examples/"

folders = readdir(example_path)
for folder in folders
    path = joinpath(example_path, folder)
    if !isdir(path) || folder in ["Benchmark", "dev"]
        continue
    end
    for file in readdir(path)
        if split(file, ".")[end] == "yaml"
            yaml_file = joinpath(path, file)
            @info "Running example $yaml_file"
            PeriLab.main(yaml_file)
        end
    end
end
