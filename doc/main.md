Robust PCA C++ library                         {#mainpage}
============

@tableofcontents


Welcome to the "Grassmann averages for scalable robust PCA C++" library.

Content
-------
This library provides a mean to compute Principal component analysis in a robust manner, using Grassmann averages. The method scales particularly well
for high dimensional problems. The full details can be found in the following article 

> **Grassmann averages for Scalable Robust PCA**, Søren Hauberg, Aasa Feragen and Michael J. Black
> Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2014

or on the web page of [the Robust PCA](http://ps.is.tuebingen.mpg.de/project/Robust_PCA). 

License
-------
All the code is licensed under the terms of the BSD-3 Clause license. A copy of this license can be found [here](http://opensource.org/licenses/BSD-3-Clause) and in the file LICENSE.txt distributed with the sources.


How to build and install
--------------

In order to build the project, you will need:
- cmake 2.8.11+
- Matlab (optional)
- Boost C++ library 1.54+

### Matlab extensions



Example code
---------------

@example example_robustpca.cpp
This is an example of how to use the non-trimmed version of the Robust PCA.