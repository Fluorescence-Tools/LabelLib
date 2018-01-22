LabelLib
========
Library for coarse-grained simulations of probes flexibly coupled to biomolecules.

Building and installation
=========================
Ubuntu
------
### Requirements
```bash
sudo apt-get install libeigen3-dev
```
On older Ubuntu versions (<17.10) Eigen3 package can contain a bug which prevents `cmake` from seeing it. To overcome this bug newer package version can be installed from Ubuntu 17.10 repository. Luckily Eigen3 is a header-only library, so it seems to work on older Ubuntu versions as well.
```bash
wget http://mirrors.kernel.org/ubuntu/pool/universe/e/eigen3/libeigen3-dev_3.3.4-3_all.deb
sudo dpkg -i libeigen3-dev_3.3.4-3_all.deb
```
If you want to build python bindings from source, you will need pybind11-dev. If pybind11 v2.2 or newer is not available for your Ubuntu version, it can be installed from source:
```bash
git clone https://github.com/pybind/pybind11.git
mkdir pybind11/build
cd pybind11/build
cmake .. -DPYBIND11_TEST=OFF
sudo make install
```

### Compiling and installation
```bash
git clone https://github.com/Fluorescence-Tools/LabelLib.git
mkdir LabelLib/build
cd LabelLib/build
cmake ..
make package
##test [optional]
make test
sudo dpkg -i FlexLabel-*-Linux.deb
##python bindings [optional]
sudo pip install ..
##cleanup [optional]
#cd ../.. && rm -rf ./LabelLib
```

Usage
=====
### C++
C++ usage example can be found at `LabelLib/FlexLabel/test/testFlexLabel.cxx`. Your own software could be compiled like this:
```bash
cd LabelLib/FlexLabel/test
g++ -std=c++14 -O3 -o FlexLabelTest testFlexLabel.cxx -lFlexLabel
./FlexLabelTest
```
Possible output:
> AV calculation took: 20.783 ms

### Python
LabelLib can be used from python code, once python bindings are installed. Usage example is available at `LabelLib/FlexLabel/python/usage.py`
```python
import LabelLib as ll
av1 = ll.dyeDensityAV1(atoms, source, 20.0, 2.0, 3.5, 0.9)
```

Citation
========
If you have used LabelLib in a scientific publication, we would appreciate citations to the following paper: [![DOI for citing LabelLib](https://img.shields.io/badge/DOI-10.1016%2Fj.sbi.2016.11.012-blue.svg)](https://doi.org/10.1016/j.sbi.2016.11.012)
> Dimura, M., Peulen, T.O., Hanke, C.A., Prakash, A., Gohlke, H. and Seidel, C.A., 2016. Quantitative FRET studies and integrative modeling unravel the structure and dynamics of biomolecular systems. Current opinion in structural biology, 40, pp.163-185.

Additional information is available in FPS toolkit paper: [![DOI for citing FPS](https://img.shields.io/badge/DOI-10.1038%2Fnmeth.2222-blue.svg)](https://doi.org/10.1038/nmeth.2222)
> Kalinin, S., Peulen, T., Sindbert, S., Rothwell, P.J., Berger, S., Restle, T., Goody, R.S., Gohlke, H. and Seidel, C.A., 2012. A toolkit and benchmark study for FRET-restrained high-precision structural modeling. Nature methods, 9(12), pp.1218-1225.
