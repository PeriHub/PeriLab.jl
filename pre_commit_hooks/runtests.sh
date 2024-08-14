#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# check Julia is in path
if ! command which julia &>/dev/null; then
  >&2 echo 'julia not found'
  exit 1
fi

# defaults
checkbounds=yes
inline=yes
coverage=false

while :; do
    case $1 in
        coverage)
            coverage=true
            ;;
        noinline)
            inline=no
            ;;
        nocheckbounds)
            checkbounds=no
            ;;
        *)
            break
    esac

    shift
done

julia --color=yes -e 'using Pkg; VERSION >= v"1.5-" && !isdir(joinpath(DEPOT_PATH[1], "registries", "General")) && Pkg.Registry.add("General")'
julia --color=yes --check-bounds="$checkbounds" --inline="$inline" --project -e "using Pkg; Pkg.test(coverage=$coverage)"
