N = 100000;
NBSteps= 10;
DIM = 10;

max_threads = 8;


mean_mex_output = zeros(max_threads, NBSteps);

for nb_threads = 2:max_threads
  do_matlab = nb_threads == 1;
  [mean_mex, var_mex, mean_matlab, var_matlab] = grassmann_pca_benches(DIM, 5, 10, nb_threads, do_matlab);
  mean_mex_output(nb_threads, :) = mean_mex;
end % for

XX = (1:NBSteps)*N/NBSteps;
h = plot(XX, mean_matlab_10, XX, mean_mex_1_10);

set(h(1), 'LineWidth', 2);
set(h(1), 'Color', [0 0 0]);

title('C++ vs MEX')
xlabel('data size')
ylabel('time (s)')