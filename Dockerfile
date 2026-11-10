# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

FROM julia:1.12 AS build

# Copy only necessary files for building
COPY src ./PeriLab/src
COPY Project.toml ./PeriLab/Project.toml

WORKDIR /PeriLab

# Install build dependencies
RUN apt-get update \
    && apt-get install -yq build-essential libxml2

ENV JULIA_CPU_TARGET="generic"

# JuliaC lives in the default environment, so it is not compiled into the PeriLab image
RUN julia -e 'using Pkg; Pkg.add("JuliaC")'
RUN julia --project -e 'using Pkg; Pkg.instantiate()'
RUN julia -e 'import JuliaC; JuliaC.main(["--output-exe", "PeriLab", "--bundle", "build", "."])'
# --trim=safe --experimental

#TODO: use alpine
FROM debian:trixie-slim AS main

WORKDIR /app

# Create the destination directory
RUN mkdir PeriLab

# Assuming /PeriLab/build is the build directory from previous stages
COPY --from=build /PeriLab/build /app/PeriLab
COPY Project.toml /app/Project.toml

# Move the build folder, set permissions, and delete the rest
RUN chmod +x /app/PeriLab/bin/PeriLab

ENV PATH="/app/PeriLab/bin:${PATH}"

# The bundle ships the MPICH the binary links against; this wrapper puts its
# mpiexec on PATH with the bundled libraries (hwloc, ...) it needs outside Julia.
RUN printf '%s\n' '#!/bin/sh' \
    'A=/app/PeriLab/share/julia/artifacts' \
    'export LD_LIBRARY_PATH="$(ls -d $A/*/lib | tr "\n" :)$LD_LIBRARY_PATH"' \
    'exec $(ls $A/*/bin/mpiexec) "$@"' > /usr/local/bin/mpiexec \
    && chmod +x /usr/local/bin/mpiexec

# --- API wrapper additions -------------------------------------------------

# Python runtime for the HTTP wrapper around the PeriLab binary
RUN apt-get update \
    && apt-get install -yq --no-install-recommends python3 python3-pip \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /app/requirements.txt
RUN pip3 install --no-cache-dir --break-system-packages -r /app/requirements.txt

COPY app /app/api

# Where uploaded job inputs and results are written; mount a volume here
# for persistence across container restarts.
RUN mkdir -p /app/simulations
ENV PERILAB_BIN="/app/PeriLab/bin/PeriLab"
ENV PERILAB_JOBS_DIR="/app/simulations"
ENV MAX_CONCURRENT_JOBS="2"

EXPOSE 3000

WORKDIR /app
CMD ["python3", "-m", "uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "3000"]
