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
    && apt-get install -yq build-essential\
    && julia --project=@. -e 'import Pkg; Pkg.add("PackageCompiler")'

RUN julia --project=@. -e 'using PackageCompiler; create_app(".", "build", executables=["PeriLab" => "main", "get_examples" => "get_examples"], force=true)'

#TODO: use alpine
FROM debian:bookworm-slim as main 

WORKDIR /app

# Create the destination directory
RUN mkdir PeriLab

# Assuming /env/build is the build directory from previous stages
COPY --from=build /env/build /app/PeriLab

# Move the build folder, set permissions, and delete the rest
RUN chmod +x /app/PeriLab/bin/PeriLab \
    && chmod +x /app/PeriLab/bin/get_examples

ENV PATH="/app/PeriLab/bin:${PATH}"

# Install SSH server and other dependencies
RUN apt-get update && apt-get install -yq openssh-server

# Configure SSH server
RUN mkdir /var/run/sshd \
    && echo 'root:root' | chpasswd \
    && sed -i'' -e's/^#PermitRootLogin prohibit-password$/PermitRootLogin yes/' /etc/ssh/sshd_config \
    && sed -i'' -e's/^#PasswordAuthentication yes$/PasswordAuthentication yes/' /etc/ssh/sshd_config \
    && sed -i'' -e's/^#PermitEmptyPasswords no$/PermitEmptyPasswords yes/' /etc/ssh/sshd_config \
    && sed -i'' -e's/^UsePAM yes/UsePAM no/' /etc/ssh/sshd_config

# Start SSH service
CMD ["/usr/sbin/sshd", "-D"]

# Expose SSH port
EXPOSE 22

