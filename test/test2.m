%-- 23/05/2014 08:56 --%
D = 100; N = 25000;
tmp = rand(D); Sigma = tmp * tmp'; X = mvnrnd(zeros(D, 1), Sigma, N);
addpath('D:\Code\gr_pca\')
tic, u0 = trimmed_grassmann_pca(X, 10, 5); toc
X = dlmread('D:\Code\robust_pca\test\mat_test.csv');
size(X)
tic, [u0, output] = trimmed_grassmann_pca(X, 10, 5); toc
size(output.initial_vectors)
% dlmwrite('d:\Code\robust_pca\test\mat_init_vectors.csv', output.initial_vectors, 'delimiter', '\t', 'precision', 8)
% dlmwrite('d:\Code\robust_pca\test\mat_test_desired_output.csv', u0, 'delimiter', '\t', 'precision', 8)




% new test
X = dlmread('D:\Code\robust_pca\test\mat_test.csv');
initial_vectors = dlmread('d:\Code\robust_pca\test\mat_test_init_vectors.csv');

addpath('D:\Code\robust_pca\build\Release\')
robustpca_m
% breakpoint can be set

algorithm_config = {};
algorithm_config.max_dimensions = 5;
algorithm_config.nb_processing_threads = 7;
algorithm_config.max_chunk_size = 1000;
algorithm_config.init = initial_vectors;
algorithm_config
tic, u = robustpca_m(X', 10, algorithm_config); toc




% trimmed grassmqn PCA
tic, u0 = trimmed_grassmann_pca(X, 10, 5, 'init', initial_vectors); toc
