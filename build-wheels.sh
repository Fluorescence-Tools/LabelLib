#!/bin/bash

# This script builds python/manylinux wheels

#Usage: 
#git clone --recursive https://github.com/Fluorescence-Tools/LabelLib.git
#docker run --rm -v `pwd`/LabelLib:/io:rw -w="/io" quay.io/pypa/manylinux1_x86_64 /io/build-wheels.sh [2|3]

set -e

VERSION_PREFIX=$1 # Could be `2` for python 2.*, `3` for python 3.x, `35` for python 3.5, etc; leave empy for all versions

DEFAULT_BIN=/opt/python/cp37-cp37m/bin
"${DEFAULT_BIN}/pip" install cmake
ln -s /opt/_internal/*/bin/cmake /usr/bin/cmake

#create the source package
git submodule update --init
${DEFAULT_BIN}/python3 setup.py sdist

# Compile wheels
for PYBIN in /opt/python/cp${VERSION_PREFIX}*/bin; do
    "${PYBIN}/pip" install numpy
    "${PYBIN}/pip" wheel ./ -w wheelhouse_tmp/ 
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse_tmp/LabelLib-*.whl; do
    auditwheel repair "$whl" -w ./wheelhouse/
done
