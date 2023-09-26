# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

FROM julia:latest
WORKDIR /env
COPY . . 

RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate();'

RUN apt-get -yq update
RUN apt-get -yq install openssh-server

# Allow mpirun as root, should only be used in container
ENV OMPI_ALLOW_RUN_AS_ROOT 1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM 1

RUN mkdir /var/run/sshd

RUN  echo 'root:root' | chpasswd
RUN sed -i'' -e's/^#PermitRootLogin prohibit-password$/PermitRootLogin yes/' /etc/ssh/sshd_config \
    && sed -i'' -e's/^#PasswordAuthentication yes$/PasswordAuthentication yes/' /etc/ssh/sshd_config \
    && sed -i'' -e's/^#PermitEmptyPasswords no$/PermitEmptyPasswords yes/' /etc/ssh/sshd_config \
    && sed -i'' -e's/^UsePAM yes/UsePAM no/' /etc/ssh/sshd_config
RUN service ssh start

WORKDIR /app/

EXPOSE 22
CMD    ["/usr/sbin/sshd", "-D"]
