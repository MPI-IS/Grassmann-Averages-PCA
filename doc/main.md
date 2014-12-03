GrassmannAveragesPCA C++ library                         {#mainpage}
============

@tableofcontents


Welcome to the "Grassmann averages for scalable robust PCA C++" library.

Content
-------
This library provides a method for computing a PCA like decomposition using Grassmann averages. The method scales particularly well
for high dimensional problems. The full details can be found in the following article 

> **Grassmann averages for Scalable Robust PCA**, Søren Hauberg, Aasa Feragen and Michael J. Black
> Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2014

or on the web page of [the Robust PCA](http://ps.is.tuebingen.mpg.de/project/Robust_PCA). Two algorithms are actually implemented:
- the Grassmann averages PCA
- the trimmed Grassmann averages PCA

The code is written in C++ that should be compatible with plain C++03 compilers. The implementations is done in the following two classes:
- @c grassmann_averages_pca::grassmann_pca for the Grassmann average PCA
- and @c grassmann_averages_pca::grassmann_pca_with_trimming for the trimmed version

These classes are templates and you should be able to run the algorithms on different types of data quite easily. 
The implementation uses several threads in order to do the processing. 

@note an C++ multithreaded implementation of EM-PCA is also provided, mainly for comparison

License
-------
All the code is licensed under the terms of the BSD-3 Clause license. 
A copy of this license can be found [here](http://opensource.org/licenses/BSD-3-Clause) and in the file LICENSE.txt 
distributed with the source code.


How to build and install
--------------
The complete build and installation instructions can be found in the README.txt file accompanying the source code.


Example code
---------------
The code contains several examples of use:
- the unit test files in C++ and in Matlab may give valuable information
- the examples attached to this documentation
- an application of the algorithm to the computation of the basis vectors of a full movie. 


@example example_grassmannpca.cpp
This is an example on how to use the non-trimmed version of the GrassmannAveragPCA. The use of the trimmed version is 
similar to this one except for the instanciation of the class @c grassmann_averages_pca::grassmann_pca_with_trimming, where the 
trimming percentage should be given.

@example video_processing.cpp
This is an example on how to adapt the GrassmannAveragPCA for subspace computation of a video sequence:

* An additional iterator class is implemented for reading the image files. This avoids having everything in memory before the algorithm starts
* An observer class that saves the subspaces as they are computed. 
