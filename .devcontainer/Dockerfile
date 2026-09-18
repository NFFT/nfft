FROM ubuntu:24.04

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Avoid warnings by switching to noninteractive.
ENV DEBIAN_FRONTEND=noninteractive

# Create the vscode user and install dependencies.
RUN groupadd --gid 2000 vscode \
    && useradd --uid 2000 --gid vscode --shell /bin/bash --create-home vscode \
    && apt update \
    && apt install -y --no-install-recommends\
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
    git \
    curl \
    wget \
    ca-certificates \
    autotools-dev \
    autoconf \
    automake \
    libtool \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog

# Set up working directory
WORKDIR /workspace

# Change ownership to vscode user
RUN chown -R vscode:vscode /workspace

# Switch to vscode user
USER vscode

# Install Julia, uv, and Python.
RUN curl -fsSL https://install.julialang.org | sh -s -- -y \
  && curl -fsSL https://astral.sh/uv/install.sh | sh \
  && /home/vscode/.local/bin/uv python install --default 3.13

# Set up git safe directory (will be configured by postCreateCommand)
CMD ["/bin/bash"]