

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
3- Matlab
4- Tests
5- Known issues


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
- Matlab 2013a+ (optional bindings). A previous version might be used, see §Matlab.

You should be able to compile the code with a free development environment. The code has been tested on:
- Linux: from a regular command line, using gcc 4.6 and 4.8
- OSX: regular command line (Clang 5.1) or XCode (5.0, 5.1), XCode should be installed on the system in order to have the compiler
- Windows: Visual Studio Express 2013

1.1 Boost 
----------------------
Usually boost can be installed through the system package manager. For instance:
- on Ubuntu: sudo apt-get install libboost-all-dev
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




1.3 Matlab (2013a+, optional)
----------------------
The library contains a target that enables the Robust PCA to be run from Matlab (MEX files). 
Matlab 2013a+ is required for these MEX files, but this requirement is mainly because of the use of the Unit Testing Framework of Matlab. 
If you are familiar with CMake, it should be easy for you to disable the Matlab unit tests while still being able to compile the 
MEX extensions. 






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

On Windows, it would be 
~~~~~~~~~~~~~~
cd $UNTAR_PATH
mkdir build
cmake -G "Visual Studio 12 Win64" ..
robustpca.sln
~~~~~~~~~~~~~~
And then you build the Release/Debug version directly from Visual Studio.


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
The library contains a target that enables the Robust PCA to be run from Matlab. 
The Matlab bindings are on by default. It is possible to disable them by providing a variable to the cmake script:

~~~~~~~~~~~~~~
cd $UNTAR_PATH
mkdir build
cmake -DWITHOUT_MATLAB=1 ..
make 
~~~~~~~~~~~~~~

The cmake script tries to find Matlab automatically. If it fails to do so, the Matlab root location can be provided 
with the following commands:
~~~~~~~~~~~~~~
cd $UNTAR_PATH
mkdir build
cmake -DMATLAB_USER_ROOT="C:\Program Files\MATLAB\R2013b" ..
make 
~~~~~~~~~~~~~~

A valid license should be available in order to be able to find Matlab. In particular, if the license is a floating network
license, the license server should be accessible during the compilation. 


2.3 Compilation options
----------------------
Depending on the platform and the available instruction set (SSE, SSE2, SSE4, AVX, ...), the compilation options may be changed in order to produce 
a final binary that is more efficient.










----------------------------------------------------------------
3- Matlab
----------------------------------------------------------------
The library contains a target that enables the Robust PCA to be run from Matlab (MEX files). 
Two files are actually needed to make the bindings work:
- robustpca_m.m: this file contains only the help
- robustpca_m.mexXXX : the extension depends on the platform. This is the binary MEX file

If you followed the previous instructions, after the build these two files can be found under $UNTAR_PATH/build, 
or $UNTAR_PATH/build/[Release|Debug] if you are using an IDE (Xcode, Visual). Several files are accompanying the 
mex file: these are the dependencies on Boost. Note that there is no such dependencies if Boost is linked against
the static library version. The dependencies are the following:
- boost_chrono
- boost_date_time
- boost_system
- boost_thread

These files are copied automatically during the Matlab mex file generation into the build directory. Note that the file
robustpca_m.m does not contain any code but only the documentation for the mex file. This file should be in the same path as 
robustpca_m.mexXXX and is also copied automatically into the build directory. Inside Matlab, typing 

  help robustpca_m 
  
prints the help for using the Robust PCA MEX file.






----------------------------------------------------------------
4- Test
----------------------------------------------------------------
Several tests exists for the library. In order to run the tests, you could run:

~~~~~~~~~~~~~~
cd $UNTAR_PATH
mkdir build
cmake <cmake options> ..
make 
make test
~~~~~~~~~~~~~~

or in case you are using a IDE (XCode, Visual Studio), "building" the target "RUN_TESTS" will also run the tests. Some tests
are running Matlab, and need a valid Matlab installation (with a valid license).






----------------------------------------------------------------
5- Known issues
----------------------------------------------------------------
- It is known that the MEX files on OSX does not honour its dependencies (boost) because of a troubleshooting with Matlab. On
  OSX currently the only possible way to run the MEX file is to use the static boost libraries.
