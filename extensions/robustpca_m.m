% Grassmann Averages for robust PCA computation.
% ROBUSTPCA_M returns the eigenvectors of a data set X, computed using the Grassmann average technique. This technique performs better
% scalability against the number of points and the dimension of the problem. 
%
% For more information about the technique, see the <a href="matlab:web('http://ps.is.tuebingen.mpg.de/project/Robust_PCA')">Robust PCA</a> page.
% 
%
%   eigen_vectors = ROBUSTPCA_M(X, 0, additionnal_options)
%
%  See also the <a href="matlab:web('http://ps.is.tuebingen.mpg.de/project/Robust_PCA')">Robust PCA</a> page
function eigen_vectors = robustpca_m(X, trimming, additionnal_options)
