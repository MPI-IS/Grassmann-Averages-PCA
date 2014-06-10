Robust PCA C++ library                         {#mainpage}
============

@tableofcontents


Welcome to the "Grassmann averages for scalable robust PCA C++" library.

Content
-------
This library provides a method for computing a PCA in a robust manner, using Grassmann averages. The method scales particularly well
for high dimensional problems. The full details can be found in the following article 

> **Grassmann averages for Scalable Robust PCA**, Søren Hauberg, Aasa Feragen and Michael J. Black
> Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2014

or on the web page of [the Robust PCA](http://ps.is.tuebingen.mpg.de/project/Robust_PCA). Two algorithms are actually implemented:
- the robust PCA
- the trimmed robust PCA

The code is written in C++ that should be compatible with plain C++03 compilers. The implementations is done in two classes:
- @c robust_pca::robust_pca_impl for the robust PCA
- and @c robust_pca::robust_pca_with_trimming_impl for the trimmed version

These classes are templates and you should be able to run the algorithms on different types of data quite easily. 
The implementation uses several threads in order to do the processing. 

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


@example example_robustpca.cpp
This is an example of how to use the non-trimmed version of the Robust PCA. The use of the trimmed version is exacly
similar to this one except for the instanciation of the class @c robust_pca::robust_pca_with_trimming_impl, where the 
trimming percentage should be given.