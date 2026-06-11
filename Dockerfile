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

RUN julia --project -e 'using Pkg; Pkg.add("JuliaC")'
RUN julia --project -e 'import JuliaC; JuliaC.main(["--output-exe", "PeriLab", "--bundle", "build", "."])'
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

# Install SSH server and other dependencies
RUN apt-get update && apt-get install -yq openssh-server libxml2

# Configure SSH server
RUN echo 'root:root' | chpasswd \
    && sed -i'' -e's/^#PermitRootLogin prohibit-password$/PermitRootLogin yes/' /etc/ssh/sshd_config \
    && sed -i'' -e's/^#PasswordAuthentication yes$/PasswordAuthentication yes/' /etc/ssh/sshd_config \
    && sed -i'' -e's/^#PermitEmptyPasswords no$/PermitEmptyPasswords yes/' /etc/ssh/sshd_config \
    && sed -i'' -e's/^UsePAM yes/UsePAM no/' /etc/ssh/sshd_config

# Start SSH service
CMD ["/usr/sbin/sshd", "-D"]

# Expose SSH port
EXPOSE 22
