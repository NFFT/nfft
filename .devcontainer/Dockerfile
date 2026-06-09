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
    git \
    curl \
    wget \
    ca-certificates \
    sudo \
    iproute2 \
    iptables \
    ipset \
    jq \
    dnsutils \
    aggregate \
    autotools-dev \
    autoconf \
    automake \
    libtool \
    pkg-config \
    clang-format \
    zsh \
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

# Set the default shell to zsh rather than sh
ENV SHELL=/bin/zsh

# Set the default editor and visual
ENV EDITOR=nano
ENV VISUAL=nano

# Default powerline10k theme
ARG ZSH_IN_DOCKER_VERSION=1.2.0
RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v${ZSH_IN_DOCKER_VERSION}/zsh-in-docker.sh)" -- \
  -p git \
  -p fzf \
  -a "source /usr/share/doc/fzf/examples/key-bindings.zsh" \
  -a "source /usr/share/doc/fzf/examples/completion.zsh" \
  -a "export PROMPT_COMMAND='history -a' && export HISTFILE=/commandhistory/.bash_history" \
  -x

# Install Julia, uv, and Python.
RUN curl -fsSL https://install.julialang.org | sh -s -- -y \
  && curl -fsSL https://astral.sh/uv/install.sh | sh \
  && /home/vscode/.local/bin/uv python install --default 3.13

RUN curl -fsSL https://raw.githubusercontent.com/rtk-ai/rtk/refs/heads/master/install.sh | sh

# Copy and set up firewall script
COPY init-firewall.sh /usr/local/bin/
USER root
RUN chmod +x /usr/local/bin/init-firewall.sh && \
  mkdir -p /etc/sudoers.d && \
  echo "vscode ALL=(root) NOPASSWD: /usr/local/bin/init-firewall.sh" > /etc/sudoers.d/vscode-firewall && \
  chmod 0440 /etc/sudoers.d/vscode-firewall
USER vscode
