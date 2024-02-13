# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

FROM julia:latest as build

WORKDIR /env

# Copy only necessary files for building
COPY src ./src
COPY Project.toml ./Project.toml

# Install build dependencies
RUN apt-get update \
    && apt-get install -yq build-essential openssh-server \
    && julia --project=@. -e 'import Pkg; Pkg.add(url="https://github.com/JTHesse/AbaqusReader.jl"); Pkg.add("PackageCompiler"); using PackageCompiler; create_app(".", "build", executables=["PeriLab" => "main", "get_examples" => "get_examples"], incremental=true, force=true)'

FROM alpine:latest as main

WORKDIR /app

# Create the destination directory
RUN mkdir PeriLab
COPY --from=build /env/build/* /app/PeriLab

# Move the build folder, set permissions, and delete the rest
RUN chmod +x /app/PeriLab/bin/PeriLab \
    && chmod +x /app/PeriLab/bin/get_examples

ENV PATH="/app/PeriLab/bin:${PATH}"
# Allow mpirun as root, should only be used in container
ENV OMPI_ALLOW_RUN_AS_ROOT 1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM 1

# Configure SSH server
RUN mkdir /var/run/sshd \
    && echo 'root:root' | chpasswd \
    && sed -i'' -e's/^#PermitRootLogin prohibit-password$/PermitRootLogin yes/' /etc/ssh/sshd_config \
    && sed -i'' -e's/^#PasswordAuthentication yes$/PasswordAuthentication yes/' /etc/ssh/sshd_config \
    && sed -i'' -e's/^#PermitEmptyPasswords no$/PermitEmptyPasswords yes/' /etc/ssh/sshd_config \
    && sed -i'' -e's/^UsePAM yes/UsePAM no/' /etc/ssh/sshd_config

# Start SSH service
RUN service ssh start

EXPOSE 22
CMD ["/usr/sbin/sshd", "-D"]

