


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
