// Copyright 2014, Max Planck Institute for Intelligent Systems.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

//!@file
//! Mex wrapper file for robust PCA


#include "mex.h"

#include <boost/numeric/ublas/storage.hpp>

#include <include/robust_pca.hpp>
#include <include/private/boost_ublas_matlab_helper.hpp>
#include <include/private/boost_ublas_matrix_helper.hpp>


//! Implements allocate and free.
template <class T>
struct matlab_matrix_allocation;


template <>
struct matlab_matrix_allocation<double>
{
  static mxArray * allocate(size_t s)
  {
    return mxCreateDoubleMatrix(s, 1, mxREAL);
  }

  static void free(mxArray *mat)
  {
    mxDestroyArray(mat);
  }
};





//! An helper class that implements the storage concept of boost::ublas
//! and performs the allocation on matlab side. It is also possible
//! to provide a matlab array.
template <class T>
class matlab_matrix_storage
{ 
  typedef T& reference;

  //! Default constructible
  matlab_matrix_storage() : pData(0)
  {
  }

  //! Size constructible
  matlab_matrix_storage(size_t s) : pData(0)
  {
    mxArray *v = matlab_matrix_allocation<T>::allocate(s);

  }



  //! Random access container
  reference operator[](size_t i)
  {
    return pData[i];
  }


private:
  T* pData;
};


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



template <class input_array_type>
bool robust_pca_dispatch(
  mxArray const* X, 
  size_t rows, 
  size_t columns, 
  size_t max_dimension, 
  int max_iterations, 
  mxArray *outputMatrix)
{
  namespace ub = boost::numeric::ublas;

  using namespace robust_pca;
  using namespace robust_pca::ublas_adaptor;
  using namespace robust_pca::ublas_matlab_helper;


  typedef external_storage_adaptor<input_array_type> input_storage_t;
  typedef ub::matrix<input_array_type, ub::column_major, input_storage_t> input_matrix_t;

  typedef external_storage_adaptor<double> output_storage_t;
  typedef ub::matrix<double, ub::column_major, output_storage_t> output_matrix_t;


  const size_t dimension = columns;

  // input data matrix, external storage.
  input_storage_t input_storage(rows*columns, mxGetPr(X));
  input_matrix_t input_data(rows, columns, input_storage);

  // output data matrix, also external storage for uBlas
  output_storage_t storageOutput(dimension * max_dimension, mxGetPr(outputMatrix));
  output_matrix_t output_eigen_vectors(max_dimension, dimension, storageOutput);




  // this is the form of the data extracted from the storage
  typedef ub::vector<double> data_t;
  typedef robust_pca_impl< data_t > robust_pca_t;

  typedef row_iter<const input_matrix_t> const_input_row_iter_t;
  typedef row_iter<output_matrix_t> output_row_iter_t;

  // should be matlab style
  std::vector<data_t> temporary_data(rows);

  robust_pca_t instance;
  return instance.batch_process(
    max_iterations,
    max_dimension,
    const_input_row_iter_t(input_data, 0),
    const_input_row_iter_t(input_data, input_data.size1()),
    temporary_data.begin(),
    output_row_iter_t(output_eigen_vectors, 0));


}

template <class input_array_type>
bool robust_pca_trimming_dispatch(
  mxArray const* X,
  size_t rows,
  size_t columns,
  size_t max_dimension,
  double trimming_percentage,
  int max_iterations,
  mxArray *outputMatrix)
{
  namespace ub = boost::numeric::ublas;

  using namespace robust_pca;
  using namespace robust_pca::ublas_adaptor;
  using namespace robust_pca::ublas_matlab_helper;


  typedef external_storage_adaptor<input_array_type> input_storage_t;
  typedef ub::matrix<input_array_type, ub::column_major, input_storage_t> input_matrix_t;

  typedef external_storage_adaptor<double> output_storage_t;
  typedef ub::matrix<double, ub::column_major, output_storage_t> output_matrix_t;


  const size_t dimension = columns;

  // input data matrix, external storage.
  input_storage_t input_storage(rows*columns, mxGetPr(X));
  input_matrix_t input_data(rows, columns, input_storage);

  // output data matrix, also external storage for uBlas
  output_storage_t storageOutput(dimension * max_dimension, mxGetPr(outputMatrix));
  output_matrix_t output_eigen_vectors(max_dimension, dimension, storageOutput);




  // this is the form of the data extracted from the storage
  typedef ub::vector<double> data_t;
  typedef robust_pca_trimming_dispatch< data_t > robust_pca_t;

  typedef row_iter<const input_matrix_t> const_input_row_iter_t;
  typedef row_iter<output_matrix_t> output_row_iter_t;

  // should be matlab style
  std::vector<data_t> temporary_data(rows);

  robust_pca_t instance(trimming_percentage / 200, 1 - trimming_percentage / 200);
  return instance.batch_process(
    max_iterations,
    max_dimension,
    const_input_row_iter_t(input_data, 0),
    const_input_row_iter_t(input_data, input_data.size1()),
    temporary_data.begin(),
    output_row_iter_t(output_eigen_vectors, 0));


}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  // arguments checking
  if (nrhs < 1 || nrhs > 4)
  {
	  mexErrMsgTxt("Incorrect number of arguments. Please consult the documentation.");
  }

  const mxArray* const X = prhs[0];
  assert(X);

  // checking the format of the data
  if (!mxIsDouble(X) && !mxIsSingle(X))
  {
	  mexErrMsgTxt("Unsupported input format (floating point required)");
  }

  if(mxIsComplex(X))
  {
    mexErrMsgTxt("Unsupported format (scalar data required)");
  }






  

  // Check the dimensions of inputs
  const size_t rows = mxGetM(X);
  const size_t columns = mxGetN(X);
  size_t dimension = columns;
  size_t max_dimension = dimension;

  // second argument is the maximum numbers of dimensions
  if (nrhs >= 2)
  {
    const mxArray* const maxDimArray = prhs[1];
    if(!mxIsNumeric(maxDimArray))
    {
      mexErrMsgTxt("Erroneous argument for the maximal dimension specification (non numeric argument)");
    }

    if(mxIsEmpty(maxDimArray))
    {
      mexErrMsgTxt("Erroneous argument for the maximal dimension specification (empty value)");
    }

    if(mxGetNumberOfElements(maxDimArray) > 1)
    {
      mexErrMsgTxt("Erroneous argument for the maximal dimension specification (non scalar)");
    }

    mxClassID classId = mxGetClassID(maxDimArray);
    if(classId == mxDOUBLE_CLASS || classId == mxSINGLE_CLASS)
    {
      //mexErrMsgTxt("Erroneous argument for the maximal dimension specification (floating point type)");
    }
    max_dimension = static_cast<size_t>(mxGetScalar(maxDimArray) + 0.5);


    if(max_dimension >= dimension)
    {
      mexErrMsgTxt("Erroneous argument for the maximal dimension specification (exceeds the dimension of the data)");
    }
  }


  // third argument is the optional trimming percentage
  bool b_trimming = false;
  double trimming_percentage = -1;
  if(nrhs == 3)
  {
    const mxArray* const trimmingArray = prhs[2];
    if(!mxIsNumeric(trimmingArray))
    {
      mexErrMsgTxt("Erroneous argument for the trimming percentage (non numeric argument)");
    }

    if(mxIsEmpty(trimmingArray))
    {
      mexErrMsgTxt("Erroneous argument for the trimming percentage (empty value)");
    }

    if(mxGetNumberOfElements(trimmingArray) > 1)
    {
      mexErrMsgTxt("Erroneous argument for the trimming percentage (non scalar)");
    }

    mxClassID classId = mxGetClassID(trimmingArray);
    if(classId == mxDOUBLE_CLASS || classId == mxSINGLE_CLASS)
    {
      //mexErrMsgTxt("Erroneous argument for the maximal dimension specification (floating point type)");
    }

    b_trimming = true;
    trimming_percentage = mxGetScalar(trimmingArray);



    if(trimming_percentage < 0 || trimming_percentage > 100)
    {
      mexErrMsgTxt("Erroneous argument for the trimming percentage (not within the range [0, 100])");
    }
  }


  int nb_iterations_max = 1000;

  if(nrhs == 4)
  {
    const mxArray* const algorithmConfiguration = prhs[3];

    if(!mxIsStruct(algorithmConfiguration))
    {
      mexErrMsgTxt("Erroneous argument for the algorithm configuration (not a structure)");
    }
    mxArray *nb_iteration_array = mxGetField(algorithmConfiguration, 0, "nb_iterations_max");
    if(nb_iteration_array != 0)
    {
      nb_iterations_max = static_cast<int>(mxGetScalar(nb_iteration_array) + 0.5);
    }
  }


  // TODO put the dimension
  // for the library to work, we need to allocate some temporary storage
  // we also allocate the output if given
  plhs[0] = mxCreateDoubleMatrix(dimension, max_dimension, mxREAL);
  mxArray *outputMatrix = plhs[0];
  assert(outputMatrix);



  
  bool result = false;
  switch(mxGetClassID(X))
  {
  case mxDOUBLE_CLASS:
  {
    if(!b_trimming)
    {
      result = robust_pca_dispatch<double>(X, rows, columns, max_dimension, nb_iterations_max, outputMatrix);
    }
    else
    {
      result = robust_pca_trimming_dispatch<double>(X, rows, columns, max_dimension, trimming_percentage, nb_iterations_max, outputMatrix);
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