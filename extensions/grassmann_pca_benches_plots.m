addpath('D:\Code\robust_pca\extensions\')
addpath('D:\Code\robust_pca\build\Release\')

redo_all = true;

N = 100000;
NBSteps= 10;
DIM = 100;

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

  save(sprintf('benches_%d.mat', DIM));
else
  S = load(sprintf('benches_%d.mat', DIM));
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

