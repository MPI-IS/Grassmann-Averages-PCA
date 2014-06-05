

RobustPCA C++ library
=====================

Homepage : http://ps.is.tuebingen.mpg.de/project/Robust_PCA
 
The library is an implementation of the Robust PCA using Grassmann averages and presented in the paper:

  "Grassmann averages for scalable robust PCA", Søren Hauberg, Aasa Feragen and Michael J. Black, CVPR 2014.

  
Table of content
----------------

1- Preparation
2- Compilation & Installation
3- Tests


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



2.2 Compilation options
----------------------
Depending on the platform and the available instruction set (SSE, SSE2, SSE4, AVX, ...), the compilation options may be changed in order to produce a more computationally efficient binary. 



----------------------------------------------------------------
3- Test
----------------------------------------------------------------