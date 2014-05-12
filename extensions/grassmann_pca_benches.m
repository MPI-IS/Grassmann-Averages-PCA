% function [mean, var, mean, var] = grassmann_pca_benches (X, K, nb_trials, nb_threads)
% 
% This function serves as a bench for the grassman PCA extension: the
% matlab implementation against the mex file.


function [mean_mex, var_mex, mean_matlab, var_matlab] = grassmann_pca_benches (X, K, nb_trials, nb_threads, do_matlab)

    algorithm_config = {};
    algorithm_config.max_dimensions = K;
    algorithm_config.nb_processing_threads = nb_threads;

    if do_matlab
        time_matlab = zeros(nb_trials, 1);
        for i = 1:nb_trials
            tic
                grassmann_pca(X, K);
            u = toc;
            time_matlab(i) = u;

        end % for
        mean_matlab = mean(time_matlab);
        var_matlab = var(time_matlab);
        
    end;
    
    xt = X';
    time_mex = zeros(nb_trials, 1);
    for i = 1:nb_trials
        tic
            robustpca_m(xt, 0, algorithm_config);
        u = toc;
        time_mex(i) = u;

    end % for
    mean_mex = mean(time_mex);
    var_mex = var(time_mex);    
    
    
end % function