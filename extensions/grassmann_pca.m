% function [vectors, output] = grassmann_pca (X, K)
% 

function [vectors, output] = grassmann_pca (X, K, varargin)
  
  %% Check input
  if (nargin == 0)
    error ('grassmann_pca: not enough input arguments');
  end % if
  
  if (nargin < 3)
    K = 1;
  end % if
  
  va = struct(varargin{:});  
  
  epsilon = 10*eps(ones(class(X)));
  
  %% Create output structure
  output.num_iters = NaN (K, 1);
  [N, D] = size (X);
  

  has_init = isfield(va, 'init');
  if (has_init)
    vectors = va.init;
    k0 = 1;
    %vectors
  elseif (isfield(va, 'extend'))
    vectors = rand (D, K, class (X)) - 0.5;
    ext = va.extend;
    X = X - (X * ext) * ext.';
    k0 = size(ext, 2)+1;
    vectors(:, 1:size(ext, 2)) = ext;
    vectors = gram_schmidt(vectors, k0);
  else
    vectors = gram_schmidt(rand (D, K, class (X)) - 0.5);
    k0 = 1;
    %vectors
  end % if  
  
  
  try
    t = tic ();
  
    for k = k0:K
       mu = vectors(:, k);
       
      %% Compute k'th principal component
      %% Initialize using a few EM iterations
      if (~has_init)
        for iter = 1:3
          dots = X * mu; % Nx1
          mu = dots.' * X; % 1xD
          %mu = trimmean(bsxfun(@times, X, dots), percent); % 1xD
          mu = mu(:) / norm(mu); % Dx1
        end % for
        dots = []; % clear up memory
      end % if

      
      %% Initialize using a few EM iterations
      % {
      for iter = 1:3
        dots = X * mu; % Nx1
        mu = dots.' * X; % 1xD
        mu = mu(:) / norm(mu); % Dx1
      end % for
      % }
      %print 'new mu starting the iterations is' mu
      
      
      %% Now the Grassmann average
      for iter = 1:N
        prev_mu = mu;

        %% Compute angles and flip
        dot_signs = sign (X * mu); % Nx1
        
        %% Compute weighted Grassmannian mean
        mu = (dot_signs).' * X; % 1xD
        mu = mu (:) / norm (mu); % Dx1

        %% Check for convergence
        if (max (abs (mu - prev_mu)) < epsilon)
          break;
        end % if
      end % for
  
      output.num_iters (k) = iter;
      vectors (:, k) = mu;
      
      %% Remove the k'th component from the data
      if (k < K)
        X = X - (X * mu) * mu.';
      end % if
    end % for

    output.time = toc (t);
  catch
    warning ('grassmann_pca: failed: %s', lasterr ());
    vectors = NaN (size (X, 2), 1);
    output.time = NaN;
  end % try_catch
end % function

