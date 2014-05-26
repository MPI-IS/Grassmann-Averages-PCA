% function [vectors, output] = grassmann_pca (X, K)
% 

function [vectors, output] = grassmann_pca (X, K)
  
  %% Check input
  if (nargin == 0)
    error ('grassmann_pca: not enough input arguments');
  elseif (nargin == 1)
    K = 1;
  end % if
  
  epsilon = 10*eps(ones(class(X)));
  
  %% Create output structure
  output.num_iters = NaN (K, 1);
  [N, D] = size (X);
  vectors = NaN (D, K, class (X));
  
  try
    t = tic ();
  
    for k = 1:K
      %% Compute k'th principal component
      mu = rand(D, 1); %G (1, :).'; % Dx1
      mu = mu / norm(mu(:));
      
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

