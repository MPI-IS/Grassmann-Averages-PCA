% Copyright 2014, Max Planck Society.
% Distributed under the BSD 3-Clause license.
% (See accompanying file LICENSE.txt or copy at
% http://opensource.org/licenses/BSD-3-Clause)


classdef grassmannpca_matlab_unit_tests < matlab.unittest.TestCase
  % some simple unit test for the grassmann pca matlab bindings
  properties
  end
    
  methods (Test)
    function testSizesReturn(testCase)
      % very simple binding test
      mat = rand(3,3);
      ret = GrassmannAveragesPCA(mat);
      testCase.verifyEqual(size(ret), size(mat));


      % by default, uses all dimensions
      mat = rand(10,3);
      ret = GrassmannAveragesPCA(mat);
      testCase.verifyEqual(size(ret, 1), size(mat, 2));
      testCase.verifyEqual(size(ret, 2), size(mat, 2));


      % we ask for the first 2 basis vectors
      algorithm_config = {};
      algorithm_config.max_dimensions = 2;

      mat = rand(100,5);
      ret = GrassmannAveragesPCA(mat, 0, algorithm_config);
      testCase.verifyEqual(size(ret, 1), size(mat, 2)); % dimension
      testCase.verifyEqual(size(ret, 2), algorithm_config.max_dimensions); % number of basis vectors
      display(ret);

      % also the case for the trimmed version
      ret = GrassmannAveragesPCA(mat, 1, algorithm_config);
      testCase.verifyEqual(size(ret, 1), size(mat, 2)); % dimension
      testCase.verifyEqual(size(ret, 2), algorithm_config.max_dimensions); % number of basis vectors
      display(ret);
    end

    function testSizesReturnMaxDimension(testCase)
      % tests the max dimension stuff
      mat = rand(3,5); % 3 vectors, dimension 5
      
      algorithm_config = {};
      algorithm_config.max_dimensions = 3;
      
      ret = GrassmannAveragesPCA(mat, 0, algorithm_config);
      testCase.verifyEqual(size(ret, 2), algorithm_config.max_dimensions);
      testCase.verifyEqual(size(ret, 1), size(mat, 2));

      algorithm_config.max_dimensions = 2;
      
      ret = GrassmannAveragesPCA(mat, 0, algorithm_config);
      testCase.verifyEqual(size(ret, 1), size(mat, 2));
      testCase.verifyEqual(size(ret, 2), 2);
    end
    
    function testInitWithWrongDimensionIsCaught(testCase)
      % tests the dimension of the init vectors
      mat = rand(3,4); % 3 vectors dimension 4
      
      algorithm_config = {};
      algorithm_config.max_dimensions = 3;
      algorithm_config.initial_vectors = rand(3, 3); % 3 vectors dimension 3, should fail
      
      testCase.verifyError(@() GrassmannAveragesPCA(mat, 0, algorithm_config), 'GrassmannAveragesPCA:configuration');

      algorithm_config.initial_vectors = rand(2, 3); % 2 vectors dimension 3, should fail
      testCase.verifyError(@() GrassmannAveragesPCA(mat, 0, algorithm_config), 'GrassmannAveragesPCA:configuration');
    end    

    function testInitWithWrongNbElementsIsCaught(testCase)
      % tests the max dimension stuff
      mat = rand(3,4); % 3 vectors of dimension 4
      
      algorithm_config = {};
      algorithm_config.max_dimensions = 1;
      algorithm_config.initial_vectors = rand(2, 4); % two initial vectors, only one requested, correct dimension
      
      testCase.verifyError(@() GrassmannAveragesPCA(mat, 0, algorithm_config), 'GrassmannAveragesPCA:configuration');
    end    

    
    
    function testTrimming2Dimensions(testCase)
      mat = rand(1000, 100); % 1000 elements, dimension 100 each

      algorithm_config = {};
      algorithm_config.max_dimensions = 2; % max dimension
      algorithm_config.nb_processing_threads = 4;
      
      ret = GrassmannAveragesPCA(mat, 5, algorithm_config); % 5 percent
      
      testCase.verifyEqual(size(ret, 1), size(mat, 2));
      testCase.verifyEqual(size(ret, 2), algorithm_config.max_dimensions);
    end
    
    
    function testTrimming100Percent(testCase)
      mat = rand(1000, 1); % 1000 elements, dimension 100 each

      algorithm_config = {};
      algorithm_config.max_dimensions = 1; % max dimension
      algorithm_config.nb_processing_threads = 4;
      
      ret = GrassmannAveragesPCA(mat, 100, algorithm_config); % 100 percent
      
      testCase.verifyEqual(size(ret, 1), size(mat, 2));
      testCase.verifyEqual(size(ret, 2), algorithm_config.max_dimensions);
      
      testCase.verifyEqual(ret(1), 1);
      
      ret = GrassmannAveragesPCA([mat, mat], 100, algorithm_config); % 100 percent
      testCase.verifyEqual(ret(1,1), ret(2, 1));
      testCase.verifyEqual(abs(ret(1,1)), 1/sqrt(2), 'AbsTol', 1.E-3);
      
    end    
    
  end
end
