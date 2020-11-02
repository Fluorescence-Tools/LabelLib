#!/bin/bash 
set -e
PY_VERSION_LIST='3.5.4 3.6.8 3.7.9 3.8.6 ""'
mkdir C:/wheelhouse
for PY_VERSION in $PY_VERSION_LIST; do
    choco install python --version=$PY_VERSION -y
    export PATH=$(cmd.exe //c "refreshenv > nul & C:\Progra~1\Git\bin\bash -c 'echo \$PATH' ")
    pip install wheel
    python setup.py bdist_wheel -d C:/wheelhouse
    choco uninstall python -y
    choco uninstall python3 -y
done
