# ################################################################
# Copyright (C) 2016-present Genome Research Ltd.
# Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
# ################################################################
# Base image
FROM ubuntu:16.04

################## METADATA ######################
LABEL base_image="ubuntu:16.04"
LABEL version="1.0"
LABEL software="Tomahawk"
LABEL software.version="0.7.0"
LABEL about.summary="Latest image for Tomahawk"
LABEL about.home="https://github.com/mklarqvist/tomahawk"
LABEL about.documentation="https://github.com/mklarqvist/tomahawk"
LABEL about.license_file="https://github.com/mklarqvist/tomahawk/blob/master/LICENSE"
LABEL about.license="MIT"
LABEL about.tags="Genomics,Linkage-disequilibrium,LD"

################## MAINTAINER ######################
MAINTAINER Marcus D. R. Klarqvist <mk819@cam.ac.uk>

ENV DEBIAN_FRONTEND noninteractive

RUN mv /etc/apt/sources.list /etc/apt/sources.list.bkp && \
    bash -c 'echo -e "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-backports main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list' && \
    cat /etc/apt/sources.list.bkp >> /etc/apt/sources.list && \
    cat /etc/apt/sources.list

RUN apt-get clean all &&  \
    apt-get update &&     \
    apt-get upgrade -y && \
    apt-get install -y  \
        pkg-config      \
        zip             \
        g++             \
        zlib1g-dev      \
        unzip           \
        curl            \
        git             \
        lsb-release     \
        liblz4-dev      \
        libssl-dev      \
        libcurl4-openssl-dev \
        liblz-dev       \
        libbz2-dev      \
        liblzma-dev     \
        autotools-dev   \
        automake        \
        cmake           \
        curl            \
        grep            \
        sed             \
        dpkg            \
        fuse            \
        git             \
        wget            \
        zip             \
        build-essential \
        pkg-config      \
        bzip2           \
        git             \
        zlib1g-dev &&   \
        apt-get clean && \
        apt-get purge && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir -p /opt
RUN cd /opt \
    && git clone https://github.com/samtools/htslib.git \
    && cd /opt/htslib \
    && make \
    && make lib-static \
    && make install

RUN cd /opt \
    && git clone https://github.com/facebook/zstd.git \
    && cd /opt/zstd/ \
    && make \
    && make install

RUN cd /opt \
    && git clone https://mklarqvist:CpgqYBz6-S@github.com/mklarqvist/tomahawk \
    && cd /opt/tomahawk/ \
    && make \
    && make install

ENV PATH=$PATH:/usr/local/bin
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/

VOLUME ["/data"]

CMD ["/bin/bash"]

WORKDIR /data