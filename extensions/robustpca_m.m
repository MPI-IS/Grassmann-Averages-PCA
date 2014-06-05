% Grassmann Averages for robust PCA computation.
% ROBUSTPCA_M returns the eigenvectors of a dataset computed using the Grassmann average technique. This technique performs better
% scalability against the number of samples and the dimension of the dataset. 
%
% For more information about the technique, see the <a href="matlab:web('http://ps.is.tuebingen.mpg.de/project/Robust_PCA')">Robust PCA</a> page.
% The general forms of the function are the following:
%
%   eigenvectors = ROBUSTPCA_M(DATA)
%   eigenvectors = ROBUSTPCA_M(DATA, trimming_percentage)
%   eigenvectors = ROBUSTPCA_M(DATA, trimming_percentage, robustpca_configuration)
%
% DATA is always an N x D matrix, where N is the number of samples and D is the dimension of the problem. 
% 
% The first call computes the regular algorithm for all dimensions, the returned eigenvector is a 
% D x D matrix where each column is an eigenvector. 
%
% The second forms trims the data by the provided percentage during the subspace computation 
% (0 yields to the regular algorithm). 
% trimming_percentage should be in the range [0, 100].
%
% The third version allows to have more control on the algorithm by providing
% a structure. The optional fields on this structure are the following:
%
% Optional fields
% ---------------
%
% - "max_dimensions": indicates that the max_dimensions eigenvectors should be computed. 
%   Defaults to D. 0 < max_dimensions <= D.
%
% - "nb_processing_threads": indicates the number of threads used for processing.
%
% - "max_chunk_size": indicates the maximum amount of data processed by each of the threads. 
%   This allows finer tuning of the processing. 
%
% - "nb_pca_steps": indicates the number of PCA iterations prior to the robust pca algorithm in case
%   random initial points are generated. 
%
% - "initial_vectors": provides the set of vectors from which to start the computation of the eigenvectors. 
%   This is a D x max_dimensions matrix. If not provided, the initial vectors will be randomly generated
%   using a uniform distribution. After the random generation, nb_pca_steps "power iterations"  (Von Mises 
%   iterations) are performed.
%   
% Example of use
% --------------
%
% * ROBUSTPCA_M on D = 100 dimensions sample of N = 100000 random Gaussian variables:
%   D = 100;
%   N = 100000;
%   tmp = rand(D); 
%   sigma = tmp * tmp'; 
%   X = mvnrnd(zeros(D, 1), sigma, N);
%   rpca_eigenvectors = ROBUSTPCA_M(X);
% 
% * Trimming K = 5 percent of the distribution:
%   rpca_eigenvectors = ROBUSTPCA_M(X, 5)
%
% * Getting only the first 3 eigenvectors:
%   robustpca_config = {};
%   robustpca_config.max_dimensions = 3;
%   rpca_eigenvectors = ROBUSTPCA_M(X, 0, robustpca_config);
%
%
%  See also the <a href="matlab:web('http://ps.is.tuebingen.mpg.de/project/Robust_PCA')">Robust PCA</a> page
function eigen_vectors = robustpca_m(X, trimming_percentage, additionnal_options)
