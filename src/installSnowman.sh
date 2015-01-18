#!/bin/bash

./autogen.sh

sparse=/xchip/gistic/Jeremiah/sparsehash-2.0.2
bamt=/broad/software/free/Linux/redhat_5_x86_64/pkgs/pezmaster31_bamtools-6708a21
seqan=/xchip/gistic/Jeremiah/seqan-trunk/core
htslib=/xchip/gistic/Jeremiah/htslib-1.1
bwalib=/xchip/gistic/Jeremiah/software/bwa
./configure --with-sparsehash=$sparse --with-bamtools=$bamt --with-seqan=$seqan --with-htslib=$htslib --with-bwalib=$bwalib

make
