FROM ubuntu:24.04

ARG TZ
ENV TZ="$TZ"

# Avoid warnings by switching to noninteractive.
ENV DEBIAN_FRONTEND=noninteractive

# Create the vscode user and install dependencies.
RUN groupadd --gid 2000 vscode \
    && useradd --uid 2000 --gid vscode --shell /bin/bash --create-home vscode \
    && apt update \
    && apt install -y --no-install-recommends \
    software-properties-common \
    && add-apt-repository universe \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
    build-essential \
    gcc-snapshot \
    gcc-14 \
    g++-14 \
    cmake \
    libfftw3-dev \
    libfftw3-single3 \
    libfftw3-double3 \
    libfftw3-long3 \
    libcunit1-dev \
    doxygen \
    doxygen-latex \
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
    && rm -rf /var/lib/apt/lists/*

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog

# Set `DEVCONTAINER` environment variable to help with orientation
ENV DEVCONTAINER=true

# Set up working directory
WORKDIR /workspace

# Change ownership to vscode user
RUN chown -R vscode:vscode /workspace

# Switch to vscode user
USER vscode

# Set the default editor and visual
ENV EDITOR=nano
ENV VISUAL=nano

# Install Julia, uv, Python, zizmor.
ENV PATH="/home/vscode/.local/bin:${PATH}"
RUN curl -fsSL https://install.julialang.org | sh -s -- -y \
  && curl -fsSL https://astral.sh/uv/install.sh | sh \
  && /home/vscode/.local/bin/uv python install --default 3.13 \
  && uv tool install zizmor

# Install RTK.
RUN curl -fsSL https://raw.githubusercontent.com/rtk-ai/rtk/refs/heads/master/install.sh | sh

# Install the CodSpeed CLI (drives the instrumented benchmarks on the CodSpeed /
# simulation track). Installs to ~/.cargo/bin and adds it to PATH.
RUN curl -fsSL https://codspeed.io/install.sh | sh
ENV PATH="/home/vscode/.cargo/bin:${PATH}"

# Install CodSpeed's Valgrind fork — the simulation/instrumentation instrument
# requires it.
USER root
RUN /home/vscode/.cargo/bin/codspeed setup --mode simulation
USER vscode
