
Grassmann Averages PCA C++ library
==================================

**Homepage** : http://ps.is.tuebingen.mpg.de/project/Robust_PCA

The library is an implementation of the Grassmann Averages for computing a PCA in robust and linear manner, presented in the paper:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**"Grassmann averages for scalable robust PCA", Soren Hauberg, Aasa Feragen and Michael J. Black, CVPR 2014.**

The *Grassmann Averages PCA* is a method for extracting the principal components from a sets of vectors, with the nice following properties

- it is of **linear complexity** wrt. the dimension of the vectors and the size of the data, which makes the method highly scalable,
- it is more **robust** to outliers than PCA in the sense that it minimizes an L1 norm instead of the L2 norm of the standard PCA.

It comes with two variants:

- the standard computation, that coincides with the PCA for normally distributed data, also referred to as the **GA**
- a trimmed variant, that is more robust to outliers, referred to the **TGA**

The paper is available from the homepage above and the details of the computations are given there.

**Content of the repository**

The repository contains the following:

- a C++ multi-threaded implementation of the *GA* and *TGA*
- a C++ multi-threaded implementation of the *EM-PCA* (for comparisons)
- binaries that computes the *GA*, *TGA* and *EM-PCA* on a set of images (frames of a video)
- Matlab bindings
- Documentation of the C++ API

Table of content
----------------

- License and terms of use
- Programs
- Preparation
- Compilation & Installation
- Matlab
- Tests
- Known issues

## 0 - License and terms of use
This program is an open-source project of the Max Planck Society. It is distributed
under the terms of the *BSD-3 Clause license*. See accompanying `LICENSE.txt` file.

## 1 - Preparation

The GrassmannAveragesPCA C++ library depends on
- Boost C++ (www.boost.org)
- cmake 2.8.11+
- Matlab 2013a+ (optional bindings). A previous version might be used, see §Matlab.

You should be able to compile the code with a free development environment. The code has been tested on:
- Linux: from a regular command line, using gcc 4.6 and 4.8
- OSX: regular command line (Clang 5.1) or XCode (5.0, 5.1), XCode should be installed on the system in order to have the compiler
- Windows: Visual Studio Express 2013

### 1.1 Boost
Usually boost can be installed through the system package manager. For instance:
- on Ubuntu: ```sudo apt-get install libboost-all-dev```
- on Mac: ```brew install boost```
- on Windows: binaries might be available (see http://www.boost.org/users/download/)

However, there are some issues using Boost when you want to compile the MEX file extension, as Matlab is also
shipped with a version of Boost that may conflict with the version you are using for the compilation.

#### Compiling boost on Windows
```
cd $BOOST_DIR
boostrap.bat
b2 --prefix=$BOOST_INSTALL_PREFIX --layout=versioned --ignore-config \
   --with-test --with-thread --with-system \
   --with-chrono --with-date_time --with-program_options install
```

where `BOOST_INSTALL_PREFIX` is the location where the library should be installed.

#### Compiling boost on Linux

If you wish to link with the static version of Boost on Linux, Boost should be compiled the following way
(by default `-fPIC` is missing for static libraries):

```
cd $BOOST_DIR
./boostrap.sh
./b2 --prefix=$BOOST_INSTALL_PREFIX --layout=versioned address-model=64 cflags=-fPIC cxxflags=-fPIC \
       --with-test --with-thread --with-system \
       --with-chrono --with-date_time --with-program_options install
```

#### Compiling boost on OSX
The command line is the following:
```
cd $BOOST_DIR
./boostrap.sh
./b2 --prefix=$BOOST_INSTALL_PREFIX --layout=versioned address-model=32_64 --ignore-config \
       --with-test --with-thread --with-system \
       --with-chrono --with-date_time --with-program_options install
```

**Note**
- On OSX, the version of boost 1.55 should be patched for Clang 5.1+ as the original release does not compile. See the patches here:
  https://github.com/boostorg/atomic/commit/e4bde20f2eec0a51be14533871d2123bd2ab9cf3.diff
  https://github.com/boostorg/atomic/commit/6bb71fdd8f7cc346d90fb14beb38b7297fc1ffd9.diff
  Boost 1.56 should include these fixes.




### 1.2 Cmake
[CMake](http://www.cmake.org) is usually available with an installer or as a package for your system.
Version 2.8.11 or above is needed.




### 1.3 Matlab (2013a+, optional)
The library contains a target that enables the GrassmannAveragesPCA to be run from Matlab (MEX files).
Matlab 2013a+ is required for these MEX files, but this requirement is mainly because of the use of the Unit Testing Framework of Matlab.
If you are familiar with CMake, it should be easy for you to disable the Matlab unit tests while still being able to compile the MEX extensions.


----------------------------------------------------------------




## 2 - Compilation & installation of the GrassmannAveragesPCA
Create a `build` directory in the directory you have cloned the project (`$source_folder`)

```
cd $source_folder
mkdir build
cmake ..
make
```

On Windows, it would be
```
cd $source_folder
mkdir build
cmake -G "Visual Studio 12 Win64" ..
GrassmannAveragesPCA.sln
```
And then you build the Release/Debug version directly from Visual Studio.


### 2.1 Indicating the location of Boost
If Boost is not found by the `cmake` scripts, an error will be thrown and indicated. Usually it means that boost is not
installed in the default library path of the system. In that case, the location of the installation directory of boost
should be provided to the `cmake` script. This location is the `BOOST_INSTALL_PREFIX` indicated in the Boost section.

```
cd $source_folder
mkdir build
cmake -DBOOST_ROOT=$BOOST_INSTALL_PREFIX ..
make
```

If Boost is not found, it may be because of the compiler id that is added to the file names (OSX). In that case, the correct
identifier may be given on the command line:

```
cd $source_folder
mkdir build
cmake -DBOOST_ROOT=$BOOST_INSTALL_PREFIX -DBoost_COMPILER=-xgcc42 ..
make
```


### 2.2 Matlab bindings
The library contains a target that enables the GrassmannAveragesPCA to be run from Matlab.
The Matlab bindings are on by default. It is possible to disable them by providing a variable to the `cmake` script:

```
cd $source_folder
mkdir build
cmake -DWITHOUT_MATLAB=1 ..
make
```

The `cmake` script tries to find Matlab automatically. If it fails to do so, the Matlab root location can be provided
with the following commands:
```
cd $source_folder
mkdir build
cmake -DMatlab_ROOT_DIR="C:\Program Files\MATLAB\R2013b" ..
make
```

A valid license should be available in order to be able to find Matlab. In particular, if the license is a
floating network license, the license server should be accessible during the compilation.


### 2.3 Compilation options
Depending on the platform and the available instruction set (SSE, SSE2, SSE4, AVX, ...), the compilation
options may be changed in order to produce a final binary that is more efficient.

----------------------------------------------------------------

## 3 - Programs

The programs were used to compute the *GA*, *TGA* and *EM-PCA* on a (large) set of images: each image is considered as a vector
and the principal components on those vectors are extracted.

They are basic applications/implementation of the C++
API where some of the parameters may be specified from the command line. They all depend on **OpenCV** for loading the images
from the files.

Their compilation is disabled automatically if *OpenCV* is not found. The location of *OpenCV* can be specified, if needed, with the
`-DOpenCVRoot=XXX` option to `cmake`.

----------------------------------------------------------------

## 4 - Matlab extensions

The library contains a target that enables the GrassmannAveragesPCA to be run from Matlab (MEX files).
Two files are actually needed to make the bindings work:
- GrassmannAveragesPCA.m: this file contains only the help
- GrassmannAveragesPCA.mexXXX : the extension depends on the platform. This is the binary MEX file

If you followed the previous instructions, after the build these two files can be found under `$source_folder/build`,
or `$source_folder/build/[Release|Debug]` if you are using an IDE (Xcode, Visual). Several files are accompanying the
mex file: these are the dependencies on Boost. Note that there is no such dependencies if Boost is linked against
the static library version. The dependencies are the following:
- boost_chrono
- boost_date_time
- boost_system
- boost_thread

These files are copied automatically during the Matlab mex file generation into the build directory.
Note that the file GrassmannAveragesPCA.m does not contain any code but only the documentation for
the mex file. This file should be in the same path as GrassmannAveragesPCA.mexXXX and is also
copied automatically into the build directory. Inside Matlab, typing

```
  help GrassmannAveragesPCA
```
prints the help for using the GrassmannAveragesPCA MEX file.

----------------------------------------------------------------


## 5 - Tests

Several tests exists for the library. In order to run the tests, you could run:

```
cd $source_folder
mkdir build
cmake <cmake options> ..
make
make test
```

or in case you are using a IDE (XCode, Visual Studio), "building" the target "RUN_TESTS" will also run the tests. Some tests are running Matlab, and need a valid Matlab installation (with a valid license).



----------------------------------------------------------------


## 6 - Known issues
- It is known that the MEX files on OSX do not honor their dependencies (boost) because of a troubleshooting with Matlab. On OSX currently the only possible way to run the MEX file is to use the static boost libraries.
- On Linux, some symbol may clash with Matlab. This is mainly because Matlab is shipped with its own `libstdc++`, that may be older than the one used for compiling the MEX extensions. In that case, you can safely run Matlab with the new version of `libstdc++` that is shipped with your system, using the following command

   ```
   LD_PRELOAD=/path/to/your/libstdc++.so matlab
   ```
 This will load the new version of `libstdc++.so` in memory prior to launching Matlab. Since those are ABI-compatible, you should not experience any problem with that
 hack.
