# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

FROM julia:1.12

WORKDIR /PeriLab

RUN apt-get update \
    && apt-get install -yq --no-install-recommends build-essential libxml2 \
    && rm -rf /var/lib/apt/lists/*

COPY Project.toml ./Project.toml

RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

COPY src ./src

# Install PeriLab itself as an app (pkg> app develop path) so `using
# PeriLab` / `main()` work from the project path at /PeriLab.
RUN julia --project=. -e 'using Pkg; Pkg.Apps.develop(path=".")'

# Ensure `julia` is on PATH inside the container; load ~/.bashrc settings
# (licenses, env vars) for interactive and `docker exec` sessions.
RUN echo 'export PATH="$HOME/.julia/bin:$PATH"' >> ~/.bashrc \
    && (source ~/.bashrc || true)

# Thin wrapper so `docker exec <container> run-perilab <args>` is the whole
# invocation -- the PeriLab App executable resolves on PATH and takes the
# passed-through CLI args (e.g. an input deck path).
RUN printf '#!/bin/bash\nexec PeriLab "$@"\n' \
    > /usr/local/bin/run-perilab \
    && chmod +x /usr/local/bin/run-perilab

COPY docker-entrypoint.sh ./docker-entrypoint.sh

# Keeps the container alive for `docker exec` -- this process does
# nothing else. Each actual simulation run is a separate `docker exec
# <container> run-perilab <args>` call, not this CMD.
ENTRYPOINT ["./docker-entrypoint.sh"]
CMD ["sleep", "infinity"]
