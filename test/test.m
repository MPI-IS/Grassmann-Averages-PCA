
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








D = 100; N = 25000;
tmp = rand(D); Sigma = tmp * tmp'; X = mvnrnd(zeros(D, 1), Sigma, N);
addpath('D:\Code\robust_pca\build\Debug\')
addpath('D:\Code\gr_pca\')
algorithm_config = {};
algorithm_config.max_dimensions = 5;
algorithm_config.nb_processing_threads = 7;
algorithm_config.max_chunk_size = 1000;
algorithm_config

%algorithm_config = 
%
%           max_dimensions: 5
%    nb_processing_threads: 7
%           max_chunk_size: 1000

tic, u0 = grassmann_pca(X, 5); toc
%Elapsed time is 1.139009 seconds.
tic, u = robustpca_m(X', 0, algorithm_config); toc
%Elapsed time is 16.603474 seconds.
algorithm_config.nb_processing_threads = 1;
tic, u = robustpca_m(X', 0, algorithm_config); toc
%Elapsed time is 65.046668 seconds.
