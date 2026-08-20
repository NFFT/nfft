FROM ubuntu:24.04
ARG TZ
ENV TZ="$TZ"

# Avoid warnings by switching to noninteractive.
ENV DEBIAN_FRONTEND=noninteractive

# Drop package extras the base image still keeps (changelogs, info pages, groff
# and lintian data). Read last, so it overrides the base image's path-include
# rules; applies to every package installed from here on.
RUN printf '%s\n' \
    'path-exclude=/usr/share/doc/*/changelog.*' \
    'path-exclude=/usr/share/info/*' \
    'path-exclude=/usr/share/groff/*' \
    'path-exclude=/usr/share/lintian/*' \
    'path-exclude=/usr/share/bug/*' \
    > /etc/dpkg/dpkg.cfg.d/zz-nfft-slim

# Create the vscode user and install dependencies. The last rm deletes mesa's software 
# rasterizer which octave-dev pulls in through libgl-dev since the container has no 
# display.
RUN groupadd --gid 2000 vscode \
    && useradd --uid 2000 --gid vscode --shell /bin/bash --create-home vscode \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
    build-essential \
    gcc-14 \
    g++-14 \
    cmake \
    libfftw3-dev \
    libfftw3-single3 \
    libfftw3-double3 \
    libfftw3-long3 \
    libcunit1-dev \
    doxygen \
    graphviz \
    octave-dev \
    clang \
    libomp-dev \
    libatomic1 \
    valgrind \
    curl \
    wget \
    ca-certificates \
    autotools-dev \
    autoconf \
    automake \
    libtool \
    pkg-config \
    clang-format \
    && rm -rf /var/lib/apt/lists/* \
    && rm -f /usr/lib/*-linux-gnu/libLLVM.so.20* \
    /usr/lib/*-linux-gnu/libLLVM-20.so \
    /usr/lib/*-linux-gnu/libgallium-*.so

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog

# Set `DEVCONTAINER` environment variable to help with orientation
ENV DEVCONTAINER=true

# Set up working directory
WORKDIR /workspace

# Change ownership to vscode user
RUN chown -R vscode:vscode /workspace

# Create Claude config directory and change ownership.
RUN mkdir -p /home/vscode/.claude && chown -R vscode:vscode /home/vscode/.claude

# Switch to vscode user
USER vscode

# Set the default editor and visual
ENV EDITOR=nano
ENV VISUAL=nano

# Install Julia, uv, Python, zizmor.
ENV PATH="/home/vscode/.local/bin:${PATH}"
# Install Julia but trim the shipped stdlib precompile cache, the
# HTML documentation and the Julia test suite.
RUN curl -fsSL https://install.julialang.org | sh -s -- -y \
    && rm -rf /home/vscode/.julia/juliaup/julia-*/share/julia/compiled \
    /home/vscode/.julia/juliaup/julia-*/share/julia/test \
    /home/vscode/.julia/juliaup/julia-*/share/doc \
    && curl -fsSL https://astral.sh/uv/install.sh | sh \
    && uv python install --default 3.13 \
    && uv tool install zizmor \
    && uv tool install ruff

# Install the CodSpeed CLI (drives the instrumented benchmarks on the CodSpeed /
# simulation track). Installs to ~/.cargo/bin and adds it to PATH.
RUN curl -fsSL https://codspeed.io/install.sh | sh
ENV PATH="/home/vscode/.cargo/bin:${PATH}"

# Install CodSpeed's Valgrind fork — the simulation/instrumentation instrument
# requires it.
USER root
RUN /home/vscode/.cargo/bin/codspeed setup --mode simulation
USER vscode
