
addpath('D:\Code\robust_pca\build\Release\')
addpath('D:\Code\gr_pca\')
addpath('D:\Code\robust_pca\extensions\')

redo_all = true;
avoid_recomputing_cpp = false;
trimmed = false;

N = 100000;
NBSteps= 10;
DIM = 100;

max_threads = 8;

timming_percent = 10;


mean_mex_output = zeros(max_threads+1, NBSteps);

if redo_all
  for nb_threads = 1:max_threads
    do_matlab = false;%nb_threads == 1;
    if trimmed
      [mean_mex, var_mex, mean_matlab, var_matlab] = grassmann_pca_trimmed_benches(DIM, 5, timming_percent, 10, nb_threads, do_matlab);
    else
      [mean_mex, var_mex, mean_matlab, var_matlab] = grassmann_pca_benches(DIM, 5, 10, nb_threads, do_matlab);
    end %
    if do_matlab
      mean_mex_output(1, :) = mean_matlab;
    end %if 
    if ~avoid_recomputing_cpp
      mean_mex_output(nb_threads+1, :) = mean_mex;
    else
      break
    end

  end % for

  if trimmed
    save(sprintf('benches_trimmed4_%d_nosignal_parallel_copy_aligned_lockfree_conservative_updates.mat', DIM));
  else
    save(sprintf('benches_%d_4_nosignal_parallel_copy_aligned_lockfree_conservative_updates.mat', DIM));
  end
else
  if trimmed
    S = load(sprintf('benches_trimmed4_%d_nosignal_parallel_copy_aligned_lockfree_conservative_updates.mat', DIM));
  else
    S = load(sprintf('benches_%d_4_nosignal_parallel_copy_aligned_lockfree_conservative_updates.mat', DIM));
  end
  mean_mex_output = S.mean_mex_output; 
end; %

hold off
XX = (1:NBSteps)*N/NBSteps;
XX = repmat(XX, max_threads+1, 1);
h = plot(XX(1,:), mean_mex_output(1, :), 'bd-');
hold on

%for i=1:(max_threads)
%  XX(i+1, :) = XX(i+1, :) + (i - max_threads/2)*N/(10*NBSteps);
%end
%plot(XX(2:end, :)', mean_mex_output(2:end, :)')
bar(XX(2:end, :)', mean_mex_output(2:end, :)')
colormap('Hot')

uistack(h,'top')
set(h(1), 'LineWidth', 2);
%set(h(1), 'Color', [0 0 0]);

title(sprintf('C++ vs MEX - dimension %d', DIM))
xlabel('data size')
ylabel('time (s)')

