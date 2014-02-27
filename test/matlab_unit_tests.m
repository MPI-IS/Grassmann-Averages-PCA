


classdef matlab_unit_tests < matlab.unittest.TestCase
  % some simple unit test
  properties
  end
    
  methods (Test)
    function test1(testCase)
      
      value = 1;
      testCase.verifyEqual(value, 1, 'Something wrong');
    end
  end
end
