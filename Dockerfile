FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y git g++ zlib1g-dev cmake libbz2-dev liblzma-dev

RUN git clone --recursive https://github.com/walaj/svaba
RUN cd svaba && ./configure && make -j$(nproc) && make -j$(nproc) install
ENV PATH "$PATH:/svaba/bin"
