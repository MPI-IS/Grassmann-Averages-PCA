% Copyright 2014, Max Planck Institute for Intelligent Systems.
% Distributed under the BSD 3-Clause license.
% (See accompanying file LICENSE.txt or copy at
% http://opensource.org/licenses/BSD-3-Clause)


classdef robustpca_matlab_unit_tests < matlab.unittest.TestCase
  % some simple unit test for the robust pca matlab bindings
  properties
  end
    
  methods (Test)
    function testSizesReturn(testCase)
      % very simple binding test
      mat = rand(3,3);
      ret = robustpca_m(mat');
      testCase.verifyEqual(size(ret), size(mat));
    end

    function testSizesReturnMaxDimension(testCase)
      % tests the max dimension stuff
      mat = rand(3,3);
      
      algorithm_config = {};
      algorithm_config.max_dimensions = 3;
      
      ret = robustpca_m(mat', 0, algorithm_config);
      testCase.verifyEqual(size(ret), size(mat));

      algorithm_config.max_dimensions = 2;
      
      ret = robustpca_m(mat', 0, algorithm_config);
      testCase.verifyEqual(size(ret, 2), 2);
    end
    
    
    
    function testComputationTimeMonothread(testCase)
      mat = rand(1000,5); % 1000 points, dimension 5

      algorithm_config = {};
      algorithm_config.max_dimensions = 5;
      algorithm_config.nb_processing_threads = 7;
      
      tic
        u = robustpca_m(mat', 0, algorithm_config);
      t1 = toc;

      % t1 = timeit(f, 1)
      
      algorithm_config.nb_processing_threads = 1;
      %t2 = timeit(f, 1)

      tic
        u = robustpca_m(mat', 0, algorithm_config);
      t2 = toc;
      
      
      display(t1)
      display(t2)
      
      
      testCase.verifyLessThan(t1, t2);
      
    end
    
    
    
    function testTrimming2Dimensions(testCase)
      mat = rand(1000, 100);

      algorithm_config = {};
      algorithm_config.max_dimensions = 2; % max dimension

      ret = robustpca_m(mat', 2, algorithm_config); % 2 percent
      
      
      testCase.verifyEqual(size(ret', 1), size(mat, 1));
      testCase.verifyEqual(size(ret', 2), algorithm_config.max_dimensions);
    end
  end
end
