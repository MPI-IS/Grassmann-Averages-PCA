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
#include <include/private/boost_ublas_external_storage.hpp>
#include <include/private/boost_ublas_row_iterator.hpp>







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
  size_t nb_pca_steps;
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

  using namespace grassmann_averages_pca;
  using namespace grassmann_averages_pca::details::ublas_helpers;


  typedef external_storage_adaptor<input_array_type> input_storage_t;
  typedef ub::matrix<input_array_type, ub::column_major, input_storage_t> input_matrix_t;

  typedef external_storage_adaptor<input_array_type> output_storage_t;
  typedef ub::matrix<input_array_type, ub::row_major, output_storage_t> output_matrix_t; // this is in fact column_major, it should be in accordance with the
                                                                               // dimension of the matrix output_eigen_vectors (we take the transpose of it)

  const size_t &dimension = algorithm_configuration.columns;
  const size_t &nb_elements = algorithm_configuration.rows;
  const size_t &max_dimension = algorithm_configuration.max_dimension;
  const size_t &max_iterations = algorithm_configuration.max_iterations;


  // input data matrix, external storage.
  input_storage_t input_storage(nb_elements * dimension, static_cast<input_array_type *>(mxGetData(X)));
  input_matrix_t input_data(nb_elements, dimension, input_storage);

  // output data matrix, also external storage for uBlas
  output_storage_t storageOutput(dimension * max_dimension, static_cast<input_array_type *>(mxGetData(outputMatrix)));
  output_matrix_t output_eigen_vectors(max_dimension, dimension, storageOutput);

  size_t nb_pca_steps = algorithm_configuration.nb_pca_steps;



  // this is the form of the data extracted from the storage
  typedef ub::vector<input_array_type> data_t;
  typedef robust_pca_impl< data_t > robust_pca_t;

  typedef row_iter<const input_matrix_t> const_input_row_iter_t;
  typedef row_iter<output_matrix_t> output_row_iter_t;

  // main instance
  robust_pca_t instance;

  if(algorithm_configuration.nb_processors > 0)
  {
    if(!instance.set_nb_processors(algorithm_configuration.nb_processors))
    {
      mexWarnMsgIdAndTxt("GrassmannAveragePCA:configuration", "Incorrect number of processors. Please consult the documentation.");
      return false;
    }
  }

  if(algorithm_configuration.max_chunk_size > 0)
  {
    if(!instance.set_max_chunk_size(algorithm_configuration.max_chunk_size))
    {
	    mexWarnMsgIdAndTxt("GrassmannAveragePCA:configuration", "Incorrect chunk size. Please consult the documentation.");
      return false;
    }
  }
  
  // initialisation vector if given
  std::vector<data_t> init_vectors;
  if(algorithm_configuration.initial_vectors != 0)
  { 
    init_vectors.resize(max_dimension);
    input_storage_t input_init_vector_storage(max_dimension*dimension, static_cast<input_array_type*>(mxGetData(algorithm_configuration.initial_vectors)));
    input_matrix_t input_init_vector_data(dimension, max_dimension, input_init_vector_storage);
    for(size_t index = 0;
        index < max_dimension;
        index++)
    {
      init_vectors[index] = ub::column(input_init_vector_data, index);
    }
    
    // if the initial vectors are set, we avoid the computation of the regular PCA.
    nb_pca_steps = 0;
  }

  if(!instance.set_nb_steps_pca(nb_pca_steps))
  {
    mexWarnMsgIdAndTxt("GrassmannAveragePCA:configuration", "Incorrect number of regular PCA steps (%d). Please consult the documentation.", nb_pca_steps);
    return false;
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

  using namespace grassmann_averages_pca;
  using namespace grassmann_averages_pca::details::ublas_helpers;


  typedef external_storage_adaptor<input_array_type> input_storage_t;
  typedef ub::matrix<input_array_type, ub::column_major, input_storage_t> input_matrix_t;

  typedef external_storage_adaptor<input_array_type> output_storage_t;
  typedef ub::matrix<input_array_type, ub::row_major, output_storage_t> output_matrix_t; // this is in fact column_major, it should be in accordance with the
                                                                               // dimension of the matrix output_eigen_vectors (we take the transpose of it)


  const size_t &dimension = algorithm_configuration.columns;
  const size_t &nb_elements = algorithm_configuration.rows;
  const size_t &max_dimension = algorithm_configuration.max_dimension;
  const size_t &max_iterations = algorithm_configuration.max_iterations;

  // input data matrix, external storage.
  input_storage_t input_storage(nb_elements*dimension, static_cast<input_array_type*>(mxGetData(X)));
  input_matrix_t input_data(nb_elements, dimension, input_storage);

  // output data matrix, also external storage for uBlas
  output_storage_t storageOutput(dimension * max_dimension, static_cast<input_array_type *>(mxGetData(outputMatrix)));
  output_matrix_t output_eigen_vectors(max_dimension, dimension, storageOutput);


  size_t nb_pca_steps = algorithm_configuration.nb_pca_steps;



  // this is the form of the data extracted from the storage
  typedef ub::vector<input_array_type> data_t;
  typedef robust_pca_with_trimming_impl< data_t > robust_pca_t;

  typedef row_iter<const input_matrix_t> const_input_row_iter_t;
  typedef row_iter<output_matrix_t> output_row_iter_t;


  // main instance
  robust_pca_t instance(algorithm_configuration.trimming_percentage / 100);


  if(algorithm_configuration.nb_processors > 0)
  {
    if(!instance.set_nb_processors(algorithm_configuration.nb_processors))
    {
      mexWarnMsgIdAndTxt("GrassmannAveragePCA:configuration", "Incorrect number of processors. Please consult the documentation.");
      return false;
    }
  }

  if(algorithm_configuration.max_chunk_size > 0)
  {
    if(!instance.set_max_chunk_size(algorithm_configuration.max_chunk_size))
    {
	    mexWarnMsgIdAndTxt("GrassmannAveragePCA:configuration", "Incorrect chunk size. Please consult the documentation.");
      return false;
    }
  }

  // initialisation vector if given
  std::vector<data_t> init_vectors;
  if(algorithm_configuration.initial_vectors != 0)
  { 
    init_vectors.resize(max_dimension);
    input_storage_t input_init_vector_storage(max_dimension*dimension, static_cast<input_array_type*>(mxGetData(algorithm_configuration.initial_vectors)));
    input_matrix_t input_init_vector_data(dimension, max_dimension, input_init_vector_storage);
    for(size_t index = 0;
        index < max_dimension;
        index++)
    {
      init_vectors[index] = ub::column(input_init_vector_data, index);
    }

    // if the initial vectors are set, we avoid the computation of the regular PCA.
    nb_pca_steps = 0;
  }


  if(!instance.set_nb_steps_pca(nb_pca_steps))
  {
    mexWarnMsgIdAndTxt("GrassmannAveragePCA:configuration", "Incorrect number of regular PCA steps (%d). Please consult the documentation.", nb_pca_steps);
    return false;
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
	  mexErrMsgIdAndTxt("GrassmannAveragePCA:configuration", "Incorrect number of arguments. Please consult the documentation.");
  }

  const mxArray* const X = prhs[0];
  assert(X);

  // checking the format of the data
  if (!mxIsDouble(X) && !mxIsSingle(X))
  {
	  mexErrMsgIdAndTxt("GrassmannAveragePCA:configuration", "Unsupported input format (floating point required)");
  }

  if(mxIsComplex(X))
  {
    mexErrMsgIdAndTxt("GrassmannAveragePCA:configuration", "Unsupported format (scalar data required)");
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
      mexErrMsgIdAndTxt("GrassmannAveragePCA:configuration", "Erroneous argument for the trimming percentage (non numeric argument)");
    }

    if(mxIsEmpty(trimmingArray))
    {
      mexErrMsgIdAndTxt("GrassmannAveragePCA:configuration", "Erroneous argument for the trimming percentage (empty value)");
    }

    if(mxGetNumberOfElements(trimmingArray) > 1)
    {
      mexErrMsgIdAndTxt("GrassmannAveragePCA:configuration", "Erroneous argument for the trimming percentage (non scalar)");
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
      mexErrMsgIdAndTxt("GrassmannAveragePCA:configuration", "Erroneous argument for the trimming percentage (not within the range [0, 100])");
    }

    b_trimming = config.trimming_percentage > 0 && config.trimming_percentage < 100;
  }


  config.max_iterations = 1000;
  config.max_chunk_size = std::numeric_limits<size_t>::max();
  config.nb_processors = 1;
  config.max_dimension = dimension;
  config.nb_pca_steps = 3;

  if(nrhs == 3)
  {
    const mxArray* const algorithmConfiguration = prhs[2];

    if(!mxIsStruct(algorithmConfiguration))
    {
      mexErrMsgIdAndTxt("GrassmannAveragePCA:configuration", "Erroneous argument for the algorithm configuration (not a structure)");
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
        mexErrMsgIdAndTxt("GrassmannAveragePCA:configuration", "Unsupported input format for initial directions (floating point required)");
      }
      if(mxIsComplex(config.initial_vectors))
      {
        mexErrMsgIdAndTxt("GrassmannAveragePCA:configuration", "Unsupported format for initial directions (scalar data required)");
      }
      
      if(mxGetM(config.initial_vectors) != dimension)
      {
        mexErrMsgIdAndTxt("GrassmannAveragePCA:configuration", "Error in the dimension of the initial values");
      }

      if(mxGetN(config.initial_vectors) != config.max_dimension)
      {
        mexErrMsgIdAndTxt("GrassmannAveragePCA:configuration", "Error in the number of the initial values provided. Should be equal to \"max_dimensions\"");
      }
      
    }

    mxArray *nb_pca_steps = mxGetField(algorithmConfiguration, 0, "nb_pca_steps");
    if(nb_pca_steps != 0)
    {
      config.nb_pca_steps = static_cast<size_t>(mxGetScalar(nb_pca_steps) + 0.5);
    }

    
    
  }



  plhs[0] = mxCreateNumericMatrix(dimension, config.max_dimension, mxGetClassID(X), mxREAL);
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
  case mxSINGLE_CLASS:
  {
    if(!b_trimming)
    {
      result = robust_pca_dispatch<float>(X, config, outputMatrix);
    }
    else
    {
      result = robust_pca_trimming_dispatch<float>(X, config, outputMatrix);
    }
    
    break;
  }
  default:
    break;
  }


  if(!result)
  {
    mexErrMsgIdAndTxt("GrassmannAveragePCA:configuration", "An error occurred in the call of the function.");
  }

}