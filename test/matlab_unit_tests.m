


classdef matlab_unit_tests < matlab.unittest.TestCase
  % some simple unit test
  properties
  end
    
  methods (Test)
    function testSizesReturn(testCase)
      mat = rand(3,3)
      %value = 1;
      %testCase.verifyEqual(value, 1, 'Something wrong');
      ret = robust_pca_m(mat)
      testCase.verifyEqual(size(ret), size(mat))
    end
  end
end
