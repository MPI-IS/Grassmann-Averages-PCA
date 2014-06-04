"Grassmann averages for scalable robust PCA C++ library", Soren Hauberg, Michael J. Black, Raffi Enficiaud



1-Compilation
2-Installation
3-Test

1-Compilation
The library depends on boost (www.boost.org). The Matlab bindings are optional.

1.1 Boost installation


1.2 Compilation
cd robustpca_path
mkdir build
cmake ..
make 


1.1 Compilation options
Depending on the platform and the available instruction set (SSE, SSE2, SSE4, AVX, ...), the compilation options may be changed in order to produce a more computationally efficient binary. 

2-Installation


3-Test
