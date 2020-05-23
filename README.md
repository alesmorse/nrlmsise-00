NRLMSISE-00 ATMOSPHERE MODEL
[![License: LGPL-3.0](https://img.shields.io/badge/License-LGPL--3.0-blue)](https://www.gnu.org/licenses/lgpl-3.0)
[![Actions Status](https://github.com/alesmorse/nrlmsise-00/workflows/CI/badge.svg)](https://github.com/alesmorse/nrlmsise-00/actions)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/alesmorse/nrlmsise-00.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/alesmorse/nrlmsise-00/context:cpp)
===========

Unofficial C++ port of the NRLMSISE-00 Model 2001 empirical atmosphere model.

## Version
This version of the code has been derived from the C source code of NRLMSISE-00 (version 2019-07-09), which in turn has been derived from the official NRLMSISE-00 version 2.0 Fortran release.

### Changes
The following changes have been made:

* A NRLMSISE-00 class is used to interface the model calls. Flags are passed when calling the constructor, which then sets the internal switches.
  * Pointer-to-implementation is used: all implementation details are moved to a dedicated class.
  * Moved data arrays and math routines (splin*) to nested namespaces.
  * Input and output structure have been removed. Calls to routines are more similar to the original Fortran version.
  * Added a function that returns directly the total atmospheric density (including anomalous oxygen contribution). This function shall be used for satellite orbit propagations.
* Added unit-testing suite (CppUnit).
* The library is making use of C++11 standard.

## Building and installation
0. Download or clone the repository:
```bash
git clone "https://github.com/alesmorse/nrlmsise-00.git"
```
1. Navigate to the repository folder:
```bash
cd nrlmsise-00
```
2. Run cmake and create a folder where to build the repository, e.g.:
```bash
cmake -S . -B ./_build
```
2. Go to the build directory and run
```bash
cd _build
make
```
3. Run the unit tests to see that everything is fine!
```bash
make test
```
4. If everything is fine, you can install the library
```bash
sudo make install
```
The ```sudo``` is there to give you the required permissions to install into your system directories (usually ```/usr/local```).


## Resources
* [Original C source code](http://www.brodo.de/space/nrlmsise/)
* [Official Fortran source code](https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/)
* [CCMC Modelweb description](http://ccmc.gsfc.nasa.gov/modelweb/atmos/nrlmsise00.html)
