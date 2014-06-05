% Grassmann Averages for robust PCA computation.
% ROBUSTPCA_M returns the eigenvectors of a data set X, computed using the Grassmann average technique. This technique performs better
% scalability against the number of points and the dimension of the problem. 
%
% For more information about the technique, see the <a href="matlab:web('http://ps.is.tuebingen.mpg.de/project/Robust_PCA')">Robust PCA</a> page.
% The general form of the function is the following:
%
%   eigenvectors = ROBUSTPCA_M(data)
%   eigenvectors = ROBUSTPCA_M(data, trimming_percentage)
%   eigenvectors = ROBUSTPCA_M(data, trimming_percentage, robustpca_configuration)
%
% * Example of use on a D = 100 dimensions sample of N = 100000 random Gaussian variables:
%   * D = 100;
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
%% BLABLABLA
%
%  See also the <a href="matlab:web('http://ps.is.tuebingen.mpg.de/project/Robust_PCA')">Robust PCA</a> page
function eigen_vectors = robustpca_m(X, trimming_percentage, additionnal_options)
