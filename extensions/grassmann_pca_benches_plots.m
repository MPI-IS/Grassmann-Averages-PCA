addpath('D:\Code\robust_pca\extensions\')
addpath('D:\Code\robust_pca\build\Release\')

redo_all = true;

N = 100000;
NBSteps= 10;
DIM = 10;

max_threads = 8;


mean_mex_output = zeros(max_threads+1, NBSteps);

if redo_all
  for nb_threads = 1:max_threads
    do_matlab = nb_threads == 1;
    [mean_mex, var_mex, mean_matlab, var_matlab] = grassmann_pca_benches(DIM, 5, 10, nb_threads, do_matlab);
    mean_mex_output(nb_threads+1, :) = mean_mex;
    if do_matlab
      mean_mex_output(1, :) = mean_matlab;
    end %if 
  end % for

  save('benches.mat');
else
  S = load('benches.mat');
  mean_mex_output = S.mean_mex_output; 
end; %
  
XX = (1:NBSteps)*N/NBSteps;
XX = repmat(XX, max_threads+1, 1);
h = plot(XX(1,:), mean_mex_output(1, :));
hold on

bar(XX(2:end, :), mean_mex_output(2:end, :))
colormap('Hot')

set(h(1), 'LineWidth', 2);
set(h(1), 'Color', [0 0 0]);

title(sprintf('C++ vs MEX - dimension %d', DIM))
xlabel('data size')
ylabel('time (s)')