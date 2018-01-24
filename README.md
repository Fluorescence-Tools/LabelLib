LabelLib
========
Library for coarse-grained simulations of probes flexibly coupled to biomolecules.

Building and installation
=========================
C++ shared library
------------------
C++ shared library can be installed from source with cmake:
```bash
git clone --recursive https://github.com/Fluorescence-Tools/LabelLib.git
mkdir LabelLib/build
cd LabelLib/build
cmake ..
sudo make install
```
On Linux you can build and install a package instead (prefered):
```bash
...
cmake .. && make package
sudo dpkg -i FlexLabel-*-Linux.deb
```

Python bindings
---------------
Python bindings can be installed via pip. Installation of the C++ library is not necessary for this.
```bash
sudo pip install LabelLib
```

Usage
=====
C++
---
C++ usage example can be found at `LabelLib/FlexLabel/test/testFlexLabel.cxx`. Your own software could be compiled like this:
```bash
cd LabelLib/FlexLabel/test
g++ -std=c++14 -O3 -o FlexLabelTest testFlexLabel.cxx -lFlexLabel
./FlexLabelTest
```
Possible output:
> AV calculation took: 20.783 ms

Python
------
LabelLib can be used from python code. Usage example is available at `LabelLib/FlexLabel/python/usage.py`
```python
import LabelLib as ll
import numpy as np
av1 = ll.dyeDensityAV1(np.zeros((4,11)), np.zeros(3), 20.0, 2.0, 3.5, 0.9)
```

Citation
========
If you have used LabelLib in a scientific publication, we would appreciate citations to the following paper: [![DOI for citing LabelLib](https://img.shields.io/badge/DOI-10.1016%2Fj.sbi.2016.11.012-blue.svg)](https://doi.org/10.1016/j.sbi.2016.11.012)
> Dimura, M., Peulen, T.O., Hanke, C.A., Prakash, A., Gohlke, H. and Seidel, C.A., 2016. Quantitative FRET studies and integrative modeling unravel the structure and dynamics of biomolecular systems. Current opinion in structural biology, 40, pp.163-185.

Additional information is available in FPS toolkit paper: [![DOI for citing FPS](https://img.shields.io/badge/DOI-10.1038%2Fnmeth.2222-blue.svg)](https://doi.org/10.1038/nmeth.2222)
> Kalinin, S., Peulen, T., Sindbert, S., Rothwell, P.J., Berger, S., Restle, T., Goody, R.S., Gohlke, H. and Seidel, C.A., 2012. A toolkit and benchmark study for FRET-restrained high-precision structural modeling. Nature methods, 9(12), pp.1218-1225.
