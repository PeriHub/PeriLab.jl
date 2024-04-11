# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Exodus

function main()
    if length(ARGS) < 1
        println("Usage: julia script.jl <filename>")
        return
    end

    file_name = ARGS[1]

    epu(file_name)
end

main()