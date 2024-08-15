# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# ARGS contains `.jl` files staged for commit, if any.

if isempty(ARGS)
    # If no Julia files changed, nothing to do.
    exit(0)
end

# Install JuliaFormatter if needed.
try
    using JuliaFormatter
catch ArgumentError
    using Pkg
    Pkg.add("JuliaFormatter")
    println("Installed JuliaFormatter package.")
    using JuliaFormatter
end

# Format each changed file.
for arg in ARGS
    JuliaFormatter.format_file(arg)
end
