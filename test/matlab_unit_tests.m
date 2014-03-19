% Copyright 2014, Max Planck Institute for Intelligent Systems.
% Distributed under the BSD 3-Clause license.
% (See accompanying file LICENSE.txt or copy at
% http://opensource.org/licenses/BSD-3-Clause)


classdef matlab_unit_tests < matlab.unittest.TestCase
  % some simple unit test for the robust pca matlab bindings
  properties
  end
    
  methods (Test)
    function testSizesReturn(testCase)
      mat = rand(3,3);
      ret = robust_pca_m(mat);
      testCase.verifyEqual(size(ret), size(mat));
    end
    
    function testSizesReturnWithDimensions(testCase)
      mat = rand(3,3);
      ret = robust_pca_m(mat, 2);
      
      sret = size(ret);
      smat = size(ret);
      
      testCase.verifyEqual(sret(1), smat(1));
      testCase.verifyEqual(sret(2), 2);
    end
  end
end
