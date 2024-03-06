# Start with an Ubuntu image
FROM ubuntu:20.04

# Avoid prompts with tzdata (timezones)
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies for htslib and svaba
RUN apt update && apt install -y \
    autoconf \
    automake \
    make \
    gcc \
    g++ \
    git \
    perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# Clone and install htslib
WORKDIR /opt
RUN git clone --recursive https://github.com/samtools/htslib.git && \
    cd htslib && \
    autoheader && \
    autoconf && \
    ./configure && \
    make && \
    make install

# Ensure shared libraries are noticed
RUN ldconfig

# Clone svaba
WORKDIR /opt
RUN git clone --recursive https://github.com/walaj/svaba.git && cd svaba && mkdir build

# Compile svaba with htslib
WORKDIR /opt/svaba/build
RUN cmake .. \
    -DHTSLIB_DIR=/usr/local

# Default command can be your application run command or just an interactive shell for testing
ENV PATH "$PATH:/svaba/build"

