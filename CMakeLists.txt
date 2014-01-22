project(RobustPCA)


set(RobustPCA_SOURCES
     include/robust_pca.hpp
     src/robust_pca.cpp)
     
     
set(RobustPCA_Matlab
     extension/matlab/robust_pca.mex)
     

add_library(robust_pca ${RobustPCA_SOURCES})