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
      ret = robustpca_m(mat);
      testCase.verifyEqual(size(ret), size(mat));
    end

    function testSizesReturnMaxDimension(testCase)
      % tests the max dimension stuff
      mat = rand(3,3);
      
      algorithm_config = {};
      algorithm_config.max_dimensions = 3;
      
      ret = robustpca_m(mat, 0, algorithm_config);
      testCase.verifyEqual(size(ret), size(mat));

      algorithm_config.max_dimensions = 2;
      
      ret = robustpca_m(mat, 0, algorithm_config);
      testCase.verifyEqual(size(ret), 2);
    
    end
    
    
    
    function testComputationTimeMonothread(testCase)
      mat = rand(1000000,50000);

      algorithm_config = {};
      algorithm_config.max_dimensions = 200;
      algorithm_config.nb_processing_threads = 1;
      
      f = @() robustpca_m(mat, 0, algorithm_config);

      t1 = timeit(f, 1);
      
      algorithm_config.nb_processing_threads = 7;
      t2 = timeit(f, 1);
      
      details(t1);
      details(t2);
      
      testCase.verifyLess(t2, t1);
      
    end
    
    
    
    function testSizesReturnWithDimensions(testCase)
      mat = rand(3,3);
      ret = robustpca_m(mat, 2);
      
      sret = size(ret);
      smat = size(ret);
      
      testCase.verifyEqual(sret(1), smat(1));
      testCase.verifyEqual(sret(2), 2);
    end
  end
end
