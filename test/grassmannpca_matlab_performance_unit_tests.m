% Copyright 2014, Max Planck Society.
% Distributed under the BSD 3-Clause license.
% (See accompanying file LICENSE.txt or copy at
% http://opensource.org/licenses/BSD-3-Clause)


classdef grassmannpca_matlab_performance_unit_tests < matlab.unittest.TestCase
  
  properties
  end
  
  methods (Test)
        
    function testComputationTimeThreads(testCase)
      mat = rand(100000,500); % 1000 points, dimension 5

      algorithm_config = {};
      algorithm_config.max_dimensions = 5;
      algorithm_config.nb_processing_threads = feature('numCores'); % this is the number of physical cores
      
      tic
        u = GrassmannAveragesPCA(mat, 0, algorithm_config);
      t1 = toc

      algorithm_config.nb_processing_threads = 1;
      tic
        u = GrassmannAveragesPCA(mat, 0, algorithm_config);
      t2 = toc
      
      
      display(t1)
      display(t2)
      
      
      testCase.verifyLessThan(t1, t2);
      
    end
  end
end