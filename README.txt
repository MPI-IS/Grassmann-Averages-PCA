

RobustPCA C++ library
=====================

Homepage : http://ps.is.tuebingen.mpg.de/project/Robust_PCA
 
The library is an implementation of the Robust PCA using Grassmann averages and presented in the paper:

  "Grassmann averages for scalable robust PCA", Søren Hauberg, Aasa Feragen and Michael J. Black, CVPR 2014.



Table of content
----------------

0- License and terms of use
1- Preparation
2- Compilation & Installation
3- Tests


----------------------------------------------------------------
0- License and terms of use
----------------------------------------------------------------
This program is an open-source project of the Max Planck Society. It is distributed
under the terms of the BSD-3 Clause license. See accompanying LICENSE.txt file. 






----------------------------------------------------------------
1- Preparation
----------------------------------------------------------------

The Robust PCA C++ library depends on 
- Boost C++ (www.boost.org) 
- cmake 2.8.11+
- Matlab (the bindings are optional).


1.1 Boost 
----------------------
Usually boost can be installed by the system package manager. For instance:
- on Ubuntu: sudo apt-get install boost-alldev
- on Mac: brew install boost
- on Windows: you have to compile boost as the binaries are not distributed. After untarring boost, you can build boost this way:


~~~~~~~~~~~~~~
cd $BOOST_DIR
./boostrap.bat
./b2 --prefix=$BOOST_INSTALL_PREFIX --layout=versioned --ignore-config --with-test --with-thread --with-chrono --with-date_time install
~~~~~~~~~~~~~~

where BOOST_INSTALL_PREFIX is the location where the library should be installed.



1.2 Cmake
----------------------
CMake (http://www.cmake.org) is usually available with an installer or as a package for your system. A version 2.8.11 or above is needed.










----------------------------------------------------------------
2- Compilation & installation of the Robust PCA
----------------------------------------------------------------
Simple create a "build" directory under the path you have untarred the library

~~~~~~~~~~~~~~
cd $UNTAR_PATH
mkdir build
cmake ..
make 
~~~~~~~~~~~~~~

2.1 Indicating the location of Boost
----------------------
If Boost is not found by the cmake scripts, an error will be thrown and indicated. Usually it means that boost is not
installed in the default library path of the system. In that case, the location of the installation directory of boost
should be provided to the cmake script. This location is the BOOST_INSTALL_PREFIX indicated in the Boost section.

~~~~~~~~~~~~~~
cd $UNTAR_PATH
mkdir build
cmake -DBOOST_ROOT=$BOOST_INSTALL_PREFIX ..
make 
~~~~~~~~~~~~~~


2.2 Matlab bindings
----------------------
The Matlab bindings are on by default. It is possible to disable them by providing a variable to the cmake script:
~~~~~~~~~~~~~~
cd $UNTAR_PATH
mkdir build
cmake -DWITHOUT_MATLAB=1 ..
make 
~~~~~~~~~~~~~~

The cmake script tries to find Matlab automatically. If it fails to do so, the Matlab root location can be provided:
~~~~~~~~~~~~~~
cd $UNTAR_PATH
mkdir build
cmake -DMATLAB_USER_ROOT="C:\Program Files\MATLAB\R2013b" ..
make 
~~~~~~~~~~~~~~

In order to be able to find Matlab, a valid license should be available. In particular, if the license is a floating network
license, the license server should be accessible during the compilation. 


2.3 Compilation options
----------------------
Depending on the platform and the available instruction set (SSE, SSE2, SSE4, AVX, ...), the compilation options may be changed in order to produce 
a final binary that is more efficient.








----------------------------------------------------------------
3- Test
----------------------------------------------------------------
