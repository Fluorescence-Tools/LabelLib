#!/bin/bash

# This script builds python/manylinux wheels

#Usage: 
#git clone --recursive https://github.com/Fluorescence-Tools/LabelLib.git
#docker run -e PYPI_USER=***** -e PYPI_PASS=***** --rm -v `pwd`/LabelLib:/io:rw quay.io/pypa/manylinux1_x86_64 /io/build-wheels.sh

VERSION_PREFIX=$1 # Could be `2` for python 2.*, `3` for python 3.x, `35` for python 3.5, etc; leave empy for all versions

PYBIN=/opt/python/cp37-cp37m/bin
"${PYBIN}/pip" install cmake
ln -s /opt/_internal/*/bin/cmake /usr/bin/cmake
"${PYBIN}/pip" install twine

# Compile wheels
for PYBIN in /opt/python/cp${VERSION_PREFIX}*/bin; do
    "${PYBIN}/pip" install numpy
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done

if [ -n "$PYPI_USER" ]; then
cat << EOF > ~/.pypirc
[distutils]
index-servers =
    pypi

[pypi]
repository: https://upload.pypi.org/legacy/
username: $PYPI_USER
password: $PYPI_PASS
EOF
/opt/_internal/*/bin/twine upload /io/wheelhouse/*
fi
