#!/bin/bash
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

function note_build_stage {
  echo "========== [$(date)] Stage '${1}' starting"
}

if [[ "$1" == "local" ]]; then
note_build_stage "Building zstd..."
if [ ! -d zstd ]; then
git clone https://github.com/facebook/zstd
fi
cd zstd
if [ ! -f lib/libzstd.so ]; then
    make -j$(nproc)
else
    echo "ZSTD already built! Skipping..."
fi
cd ..

note_build_stage "Building OpenSSL..."
if [ ! -d openssl ]; then
git clone https://github.com/openssl/openssl.git
fi
cd openssl
if [ ! -f libssl.so ]; then
    ./config
    make -j$(nproc)
else
    echo "OpenSSL already built! Skipping..."
fi
cd ..

note_build_stage "Building curl..."
if [ ! -d curl-7.61.0 ]; then
wget https://dl.uxnr.de/mirror/curl/curl-7.61.0.tar.gz
tar -xzvf curl-7.61.0.tar.gz
fi
cd curl-7.61.0
if [ ! -f lib/.libs/libcurl.so ]; then
    CPPFLAGS="-I${PWD}/../openssl/include" LDFLAGS="-L${PWD}/../openssl" ./configure
    make -j$(nproc)
else
    echo "curl already built! Skipping..."
fi
cd ..

note_build_stage "Building htslib"
if [ ! -d htslib ]; then
git clone https://github.com/samtools/htslib.git
fi
cd htslib
if [ ! -f htslib.so ]; then
    # Temporary fix for broken htslib compatibility with C++
    git checkout 1832d3a1b75133e55fb6abffc3f50f8a6ed5ceae
    autoheader && autoconf && ./configure CPPFLAGS="-I/usr/local/include/" LDFLAGS="-L/usr/local/lib/ -lcurl" && make -j$(nproc)
else
    echo "htslib already built! Skipping..."
fi
cd ..

else # Install with sudo
if [ "$(uname)" == "Darwin" ]; then
    # Update package list
    ################################################################################
    brew update
    # Install generic dependencies
    ################################################################################
    note_build_stage "Install openssl"
    brew install openssl
     note_build_stage "Install zstd"
    brew install zstd
     note_build_stage "Install htslib"
    brew install htslib
else
    # Update package list
    ################################################################################
    note_build_stage "Update package list"
    sudo -H apt-get -qq -y update

    # Install generic dependencies
    ################################################################################
    note_build_stage "Update misc. dependencies"
    sudo -H apt-get -y install pkg-config zip g++ zlib1g-dev unzip curl git lsb-release liblz4-dev

    # Install htslib dependencies
    ################################################################################
    note_build_stage "Install htslib dependencies"
    sudo -H apt-get -y install libssl-dev libcurl4-openssl-dev liblz-dev libbz2-dev liblzma-dev

    # Install from Github
    ################################################################################
    CURDIR=`echo $PWD`
    TMPDIR=`mktemp -d -t`
    cd ${TMPDIR}
    # Install zstd
    ################################################################################
    note_build_stage "Install zstd"
    git clone https://github.com/facebook/zstd
    cd zstd &&  make -j$(nproc) && sudo make install
    cd ..
    # Install htslib
    ################################################################################
    note_build_stage "Install htslib"
    git clone https://github.com/samtools/htslib.git
    cd htslib
    # Temporary fix for broken htslib compatibility with C++
    git checkout 1832d3a1b75133e55fb6abffc3f50f8a6ed5ceae
    # Build and install
    autoheader && autoconf && ./configure && make -j$(nproc) && sudo make install
    cd ${CURDIR}
fi
fi

# Install Tachyon
################################################################################
note_build_stage "Building Tomahawk"
make clean; make -j$(nproc)