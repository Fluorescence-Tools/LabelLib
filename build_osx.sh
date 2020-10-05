#!/bin/bash 
set -e

pyenv global `pyenv versions | grep "^ *$PYTHON" | head -n1`
mkdir build && cd build
cmake -DPYTHON_EXECUTABLE:FILEPATH=`pyenv which python$PYTHON` ..
make
