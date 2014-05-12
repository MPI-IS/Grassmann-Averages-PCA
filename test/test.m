
addpath('D:\Code\gr_pca\')
mat = rand(10000,50);
%[toto, output] = trimmed_grassmann_pca(mat, 10, 5)
tic, [toto, output] = trimmed_grassmann_pca(mat, 10, 5); toc

addpath('D:\Code\robust_pca\build\Release\')
algorithm_config = {};
algorithm_config.max_dimensions = 5;
algorithm_config.nb_processing_threads = 7;
algorithm_config.max_chunk_size = 1000;
algorithm_config
tic, u = robustpca_m(mat', 10, algorithm_config); toc





addpath('D:\Code\gr_pca\')
mat = rand(100000,500);
%[toto, output] = trimmed_grassmann_pca(mat, 10, 5)
tic, [toto, output] = grassmann_pca(mat, 5); toc

addpath('D:\Code\robust_pca\build\Release\')
algorithm_config = {};
algorithm_config.max_dimensions = 5;
algorithm_config.nb_processing_threads = 7;
algorithm_config.max_chunk_size = 1000;
algorithm_config
tic, u = robustpca_m(mat', 0, algorithm_config); toc


