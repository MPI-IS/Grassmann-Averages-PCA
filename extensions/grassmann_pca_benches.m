% function [mean_mex, var_mex, mean_matlab, var_matlab] = grassmann_pca_benches (D, K, nb_trials, nb_threads)
% 
% This function serves as a bench for the grassman PCA extension: the
% matlab implementation against the mex file.


function [mean_mex, var_mex, mean_matlab, var_matlab] = grassmann_pca_benches (D, K, nb_trials, nb_threads, do_matlab)

  N = 100000;
  NB_steps = 10;
  %P = N/NB_steps;
  tmp = rand(D); 
  Sigma = tmp * tmp'; 
  X = mvnrnd(zeros(D, 1), Sigma, N);

  vectors = gram_schmidt(rand (D, K, class (X)) - 0.5);

  algorithm_config = {};
  algorithm_config.max_dimensions = K;
  algorithm_config.nb_processing_threads = nb_threads;
  %algorithm_config.max_chunk_size = 1000;
  algorithm_config.initial_vectors = vectors;

  
  mean_matlab = zeros(NB_steps, 1);
  var_matlab  = zeros(NB_steps, 1);

  if do_matlab
    display('Matlab version');
    time_matlab = zeros(nb_trials, 1);
    
    for j = 1:NB_steps
      display(sprintf('\tStep %d - #elements = %d', j, (N*j/NB_steps)))
      xt = X(1:(N*j/NB_steps), :);
      for i = 1:nb_trials
        tic
          out_soren = grassmann_pca(xt, K, 'init', vectors);
        u = toc;
        time_matlab(i) = u;

      end % for
      mean_matlab(j) = mean(time_matlab);
      var_matlab(j) = var(time_matlab);
    end % for
    display(mean_matlab')
    display(var_matlab')
  else
    display('Matlab results for comparison')
    out_soren = grassmann_pca(X, K, 'init', vectors);  
  end;
  
  time_mex = zeros(nb_trials, 1);
  mean_mex = zeros(NB_steps, 1);
  var_mex  = zeros(NB_steps, 1);

  display('C++ version');
  for j = 1:NB_steps
    display(sprintf('\tStep %d - #elements = %d', j, (N*j/NB_steps)))
    xt = (X(1:(N*j/NB_steps), :));
    for i = 1:nb_trials
      tic
        out_raffi = robustpca_m(xt, 0, algorithm_config);
      u = toc;
      time_mex(i) = u;

    end % for
    mean_mex(j) = mean(time_mex);
    var_mex(j) = var(time_mex);    
  end % for
  display(mean_mex')
  display(var_mex')
  if true
    display('max error=')
    display(max(abs(out_raffi - out_soren)))
  end
end % function