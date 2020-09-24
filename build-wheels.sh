#!/bin/bash

# This script builds python/manylinux wheels

#Usage: 
#git clone --recursive https://github.com/Fluorescence-Tools/LabelLib.git
#docker run -e PYPI_USER=***** -e PYPI_PASS=***** --rm -v `pwd`/LabelLib:/io:rw quay.io/pypa/manylinux1_x86_64 /io/build-wheels.sh

PYBIN=/opt/python/cp37-cp37m/bin
"${PYBIN}/pip" install cmake
ln -s /opt/_internal/*/bin/cmake /usr/bin/cmake
"${PYBIN}/pip" install twine

# Compile wheels
for PYBIN in /opt/python/*/bin; do
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
( cd /io; git submodule update --init; ${PYBIN}/python3 setup.py sdist)
/opt/_internal/*/bin/twine upload /io/wheelhouse/* /io/dist/*.tar.gz
fi
