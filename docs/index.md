# NRLMSISE-00 ATMOSPHERE MODEL
Unofficial C++ port of the NRLMSISE-00 Model 2001 empirical atmosphere model.

## Introduction
The NRLMSISE-00 model was developed by Mike Picone, Alan Hedin, and Doug Drob. They also wrote a NRLMSISE-00 distribution package in FORTRAN which is available at [CCMC Modelweb description](http://ccmc.gsfc.nasa.gov/modelweb/atmos/nrlmsise00.html).

## Library

* A NRLMSISE-00 class (atmos::CNrlmsise00) is used to interface the model calls. Flags are passed when calling the constructor, which then sets the internal switches.
  * Pointer-to-implementation is used: all implementation details are moved to a dedicated class, atmos::CNrlmsise00_p).
  * Moved math routines (splin*) to nested namespace, atmos::math.
  * Input and output structure have been removed. Calls to routines are more similar to the original Fortran version.
  * Added a function that returns directly the total atmospheric density (including anomalous oxygen contribution), atmos::CNrlmsise00::density. This function shall be used for satellite orbit propagations.
* Added unit-testing suite (CppUnit).
* The library is making use of C++11 standard.