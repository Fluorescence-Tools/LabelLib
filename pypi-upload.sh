#!/bin/bash

# This script builds python/manylinux wheels

#Usage: 
#git clone --recursive https://github.com/Fluorescence-Tools/LabelLib.git
#docker run -e PYPI_USER=***** -e PYPI_PASS=***** --rm -v `pwd`/LabelLib:/io:rw -w="/io" quay.io/pypa/manylinux1_x86_64 /io/pypi-upload.sh

set -e
DEFAULT_BIN=/opt/python/cp37-cp37m/bin

if [ -z "${PYPI_USER}" ] || [ -z "${PYPI_PASS}" ]; then
    echo "ERROR PYPI credintials are not set!"
    exit 1
fi

cat << EOF > ~/.pypirc
[distutils]
index-servers =
    pypi

[pypi]
repository: https://upload.pypi.org/legacy/
username = $PYPI_USER
password = $PYPI_PASS
EOF

"${DEFAULT_BIN}/pip" install twine
"${DEFAULT_BIN}/twine" upload ./wheelhouse/LabelLib-*.whl ./dist/LabelLib-*.tar.gz
