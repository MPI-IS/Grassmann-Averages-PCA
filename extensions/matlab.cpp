// Copyright 2014, Max Planck Institute for Intelligent Systems.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

//!@file
//! Mex wrapper file for robust PCA

// this first include is a small workaround for clang503/xcode5.1: algorithm should be included prior to mex.h
#include <algorithm>

#include "mex.h"

#include <boost/numeric/ublas/storage.hpp>

#include <include/robust_pca.hpp>
#include <include/robust_pca_trimming.hpp>
#include <include/private/boost_ublas_matlab_helper.hpp>
#include <include/private/boost_ublas_matrix_helper.hpp>







template <class output_type, class input_type>
output_type get_matlab_array_value_dispatch(mxArray const* array_, size_t index_row, size_t index_column)
{
  input_type* p_array = static_cast<input_type*>(mxGetData(array_));
  if(p_array == 0)
  {
    throw std::runtime_error("Unable to retrieve a pointer to the typed array");
  }

  // column major
  return p_array[index_row + mxGetN(array_) * index_column];
}

template <class output_type>
output_type get_matlab_array_value(mxArray const* array_, size_t index_row, size_t index_column)
{
  assert(index_row < mxGetM(array_));
  assert(index_column < mxGetN(array_));

  mxClassID classId = mxGetClassID(array_);
  switch(classId)
  {
  case mxDOUBLE_CLASS:
    return get_matlab_array_value_dispatch<output_type, double>(array_, index_row, index_column);
  case mxSINGLE_CLASS:
    return get_matlab_array_value_dispatch<output_type, float>(array_, index_row, index_column);
  default:
    throw std::runtime_error("Unable to dispatch to the correct type");
  }
}


struct s_algorithm_configuration
{
  size_t rows;
  size_t columns;

  size_t max_dimension;
  size_t max_iterations;
  size_t max_chunk_size;
  size_t nb_processors;
  double trimming_percentage;
  
  mxArray *initial_vectors;

  // we should also add the initial value of the eigen vectors
};


template <class input_array_type>
bool robust_pca_dispatch(
  mxArray const* X, 
  s_algorithm_configuration const& algorithm_configuration, 
  mxArray *outputMatrix)
{
  namespace ub = boost::numeric::ublas;

  using namespace robust_pca;
  using namespace robust_pca::ublas_adaptor;
  using namespace robust_pca::ublas_matlab_helper;


  typedef external_storage_adaptor<input_array_type> input_storage_t;
  typedef ub::matrix<input_array_type, ub::row_major, input_storage_t> input_matrix_t;

  typedef external_storage_adaptor<double> output_storage_t;
  typedef ub::matrix<double, ub::row_major, output_storage_t> output_matrix_t;

  const size_t &dimension = algorithm_configuration.columns;
  const size_t &nb_elements = algorithm_configuration.rows;
  const size_t &max_dimension = algorithm_configuration.max_dimension;
  const size_t &max_iterations = algorithm_configuration.max_iterations;


  // input data matrix, external storage.
  input_storage_t input_storage(nb_elements * dimension, mxGetPr(X));
  input_matrix_t input_data(nb_elements, dimension, input_storage);

  // output data matrix, also external storage for uBlas
  output_storage_t storageOutput(dimension * max_dimension, mxGetPr(outputMatrix));
  output_matrix_t output_eigen_vectors(max_dimension, dimension, storageOutput);




  // this is the form of the data extracted from the storage
  typedef ub::vector<double> data_t;
  typedef robust_pca_impl< data_t > robust_pca_t;

  typedef row_iter<const input_matrix_t> const_input_row_iter_t;
  typedef row_iter<output_matrix_t> output_row_iter_t;

  // main instance
  robust_pca_t instance;

  if(algorithm_configuration.nb_processors > 0)
  {
    if(!instance.set_nb_processors(algorithm_configuration.nb_processors))
    {
      mexWarnMsgTxt("Incorrect number of processors. Please consult the documentation.");
      return false;
    }
  }

  if(algorithm_configuration.max_chunk_size > 0)
  {
    if(!instance.set_max_chunk_size(algorithm_configuration.max_chunk_size))
    {
	    mexWarnMsgTxt("Incorrect chunk size. Please consult the documentation.");
      return false;
    }
  }
  
  // initialisation vector if given
  std::vector<data_t> init_vectors;
  if(algorithm_configuration.initial_vectors != 0)
  { 
    init_vectors.resize(max_dimension);
    input_storage_t input_init_vector_storage(max_dimension*dimension, mxGetPr(algorithm_configuration.initial_vectors ));
    input_matrix_t input_init_vector_data(max_dimension, dimension, input_init_vector_storage);
    size_t index = 0;
    for(const_input_row_iter_t it(input_init_vector_data, 0), ite(input_init_vector_data, max_dimension);
        it != ite;
        ++it)
    {
      init_vectors[index++] = *it;
    }
    
  }

  return instance.batch_process(
    max_iterations,
    max_dimension,
    const_input_row_iter_t(input_data, 0),
    const_input_row_iter_t(input_data, input_data.size1()),
    output_row_iter_t(output_eigen_vectors, 0),
    algorithm_configuration.initial_vectors ? &init_vectors: 0);


}

template <class input_array_type>
bool robust_pca_trimming_dispatch(
  mxArray const* X,
  s_algorithm_configuration const& algorithm_configuration, 
  mxArray *outputMatrix)
{
  namespace ub = boost::numeric::ublas;

  using namespace robust_pca;
  using namespace robust_pca::ublas_adaptor;
  using namespace robust_pca::ublas_matlab_helper;


  typedef external_storage_adaptor<input_array_type> input_storage_t;
  typedef ub::matrix<input_array_type, ub::row_major, input_storage_t> input_matrix_t;

  typedef external_storage_adaptor<double> output_storage_t;
  typedef ub::matrix<double, ub::row_major, output_storage_t> output_matrix_t;


  const size_t &dimension = algorithm_configuration.columns;
  const size_t &nb_elements = algorithm_configuration.rows;
  const size_t &max_dimension = algorithm_configuration.max_dimension;
  const size_t &max_iterations = algorithm_configuration.max_iterations;

  // input data matrix, external storage.
  input_storage_t input_storage(nb_elements*dimension, mxGetPr(X));
  input_matrix_t input_data(nb_elements, dimension, input_storage);

  // output data matrix, also external storage for uBlas
  output_storage_t storageOutput(dimension * max_dimension, mxGetPr(outputMatrix));
  output_matrix_t output_eigen_vectors(max_dimension, dimension, storageOutput);




  // this is the form of the data extracted from the storage
  typedef ub::vector<double> data_t;
  typedef robust_pca_with_trimming_impl< data_t > robust_pca_t;

  typedef row_iter<const input_matrix_t> const_input_row_iter_t;
  typedef row_iter<output_matrix_t> output_row_iter_t;


  // main instance
  robust_pca_t instance(algorithm_configuration.trimming_percentage / 100);


  if(algorithm_configuration.nb_processors > 0)
  {
    if(!instance.set_nb_processors(algorithm_configuration.nb_processors))
    {
      mexWarnMsgTxt("Incorrect number of processors. Please consult the documentation.");
      return false;
    }
  }

  if(algorithm_configuration.max_chunk_size > 0)
  {
    if(!instance.set_max_chunk_size(algorithm_configuration.max_chunk_size))
    {
	    mexWarnMsgTxt("Incorrect chunk size. Please consult the documentation.");
      return false;
    }
  }

  // initialisation vector if given
  std::vector<data_t> init_vectors;
  if(algorithm_configuration.initial_vectors != 0)
  { 
    init_vectors.resize(max_dimension);
    input_storage_t input_init_vector_storage(max_dimension*dimension, mxGetPr(algorithm_configuration.initial_vectors ));
    input_matrix_t input_init_vector_data(max_dimension, dimension, input_init_vector_storage);
    size_t index = 0;
    for(const_input_row_iter_t it(input_init_vector_data, 0), ite(input_init_vector_data, max_dimension);
        it != ite;
        ++it)
    {
      init_vectors[index++] = *it;
    }
    
  }


  return instance.batch_process(
    max_iterations,
    max_dimension,
    const_input_row_iter_t(input_data, 0),
    const_input_row_iter_t(input_data, input_data.size1()),
    output_row_iter_t(output_eigen_vectors, 0),
    algorithm_configuration.initial_vectors ? &init_vectors: 0);


}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  // arguments checking
  if (nrhs < 1 || nrhs > 4)
  {
	  mexErrMsgIdAndTxt("RobustPCA:configuration", "Incorrect number of arguments. Please consult the documentation.");
  }

  const mxArray* const X = prhs[0];
  assert(X);

  // checking the format of the data
  if (!mxIsDouble(X) && !mxIsSingle(X))
  {
	  mexErrMsgIdAndTxt("RobustPCA:configuration", "Unsupported input format (floating point required)");
  }

  if(mxIsComplex(X))
  {
    mexErrMsgIdAndTxt("RobustPCA:configuration", "Unsupported format (scalar data required)");
  }





  s_algorithm_configuration config;



  

  // Check the dimensions of inputs
  config.rows = mxGetM(X);
  config.columns = mxGetN(X);
  
  size_t dimension = config.columns;
  config.initial_vectors  = 0;


  // third argument is the optional trimming percentage
  bool b_trimming = false;
  config.trimming_percentage = -1;
  if(nrhs >= 2)
  {
    const mxArray* const trimmingArray = prhs[1];
    if(!mxIsNumeric(trimmingArray))
    {
      mexErrMsgIdAndTxt("RobustPCA:configuration", "Erroneous argument for the trimming percentage (non numeric argument)");
    }

    if(mxIsEmpty(trimmingArray))
    {
      mexErrMsgIdAndTxt("RobustPCA:configuration", "Erroneous argument for the trimming percentage (empty value)");
    }

    if(mxGetNumberOfElements(trimmingArray) > 1)
    {
      mexErrMsgIdAndTxt("RobustPCA:configuration", "Erroneous argument for the trimming percentage (non scalar)");
    }

    mxClassID classId = mxGetClassID(trimmingArray);
    if(classId == mxDOUBLE_CLASS || classId == mxSINGLE_CLASS)
    {
      //mexErrMsgTxt("Erroneous argument for the maximal dimension specification (floating point type)");
    }

    b_trimming = true;
    config.trimming_percentage = mxGetScalar(trimmingArray);

    if(config.trimming_percentage < 0 || config.trimming_percentage > 100)
    {
      mexErrMsgIdAndTxt("RobustPCA:configuration", "Erroneous argument for the trimming percentage (not within the range [0, 100])");
    }

    b_trimming = config.trimming_percentage > 0 && config.trimming_percentage < 100;
  }


  config.max_iterations = 1000;
  config.max_chunk_size = std::numeric_limits<size_t>::max();
  config.nb_processors = 1;
  config.max_dimension = dimension;

  if(nrhs == 3)
  {
    const mxArray* const algorithmConfiguration = prhs[2];

    if(!mxIsStruct(algorithmConfiguration))
    {
      mexErrMsgIdAndTxt("RobustPCA:configuration", "Erroneous argument for the algorithm configuration (not a structure)");
    }
    
    mxArray *nb_iteration_array = mxGetField(algorithmConfiguration, 0, "nb_iterations_max");
    if(nb_iteration_array != 0)
    {
      config.max_iterations = static_cast<int>(mxGetScalar(nb_iteration_array) + 0.5);
    }

    mxArray *max_chunk_size_array = mxGetField(algorithmConfiguration, 0, "max_chunk_size");
    if(max_chunk_size_array != 0)
    {
      config.max_chunk_size = static_cast<size_t>(mxGetScalar(max_chunk_size_array) + 0.5);
    }

    mxArray *nb_processing_threads_array = mxGetField(algorithmConfiguration, 0, "nb_processing_threads");
    if(nb_processing_threads_array != 0)
    {
      config.nb_processors = static_cast<size_t>(mxGetScalar(nb_processing_threads_array) + 0.5);
    }

    mxArray *nb_max_dimensions = mxGetField(algorithmConfiguration, 0, "max_dimensions");
    if(nb_max_dimensions != 0)
    {
      config.max_dimension = static_cast<size_t>(mxGetScalar(nb_max_dimensions) + 0.5);
    }
    
    config.initial_vectors = mxGetField(algorithmConfiguration, 0, "initial_vectors");
    if(config.initial_vectors != 0)
    {
      if (!mxIsDouble(config.initial_vectors) && !mxIsSingle(config.initial_vectors))
      {
        mexErrMsgIdAndTxt("RobustPCA:configuration", "Unsupported input format for initial directions (floating point required)");
      }
      if(mxIsComplex(config.initial_vectors))
      {
        mexErrMsgIdAndTxt("RobustPCA:configuration", "Unsupported format for initial directions (scalar data required)");
      }
      
      if(mxGetN(config.initial_vectors) != dimension)
      {
        mexErrMsgIdAndTxt("RobustPCA:configuration", "Error in the dimension of the initial values");
      }

      if(mxGetM(config.initial_vectors) != config.max_dimension)
      {
        mexErrMsgIdAndTxt("RobustPCA:configuration", "Error in the number of the initial values provided. Should be equal to \"max_dimensions\"");
      }
      
    }

    
    
  }


  // TODO put the dimension
  // for the library to work, we need to allocate some temporary storage
  // we also allocate the output if given
  plhs[0] = mxCreateDoubleMatrix(dimension, config.max_dimension, mxREAL);
  mxArray *outputMatrix = plhs[0];
  assert(outputMatrix);



  
  bool result = false;
  switch(mxGetClassID(X))
  {
  case mxDOUBLE_CLASS:
  {
    if(!b_trimming)
    {
      result = robust_pca_dispatch<double>(X, config, outputMatrix);
    }
    else
    {
      result = robust_pca_trimming_dispatch<double>(X, config, outputMatrix);
    }
    
    break;
  }
  default:
    break;
  }


  if(!result)
  {
    mexErrMsgTxt("Robust PCA: an error occurred in the call of the function.");
  }

}