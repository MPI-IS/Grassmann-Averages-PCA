// Mex wrapper file for robust PCA


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


#if 0


//! A simple vector used with two pointers. 
//! The vector is read-only
template <class T>
struct pointer_vector_t
{
  typedef typename boost::add_reference<T>::type reference;
  //! @param length_ length of the vector (number of elements)
  pointer_vector_t(T * initial_position_, size_t const length_) : 
    initial_position(initial_position_)
    length(length_)
  {}


  reference operator[](size_t const index)
  {
    return 
  }


private:
  T* initial_position;
  size_t length;
};
#endif









void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  namespace ub = boost::numeric::ublas;
  using namespace robust_pca;
  using namespace robust_pca::ublas_adaptor;
  using namespace robust_pca::ublas_matlab_helper;


  // arguments checking
  if (nrhs != 1) 
  {
	  mexErrMsgTxt("One input argument required.");
  }
  else if (nlhs > 1) 
  {
	  mexErrMsgTxt("Too many output arguments.");
  }

  const mxArray* const X = prhs[0];
  assert(X);

  // checking the format of the data
  if (!mxIsDouble(X) && !mxIsSingle(X))
  {
	  mexErrMsgTxt("Unsupported format");
  }

  if(mxIsComplex(X))
  {
    mexErrMsgTxt("Unsupported format (should be scalar)");
  }



  

  // Check the dimensions of inputs
  const size_t rows = mxGetM(X);
  const size_t columns = mxGetN(X);
  const size_t &dimension = columns; // an alias


  if(mxIsDouble(X))
  {
    typedef readonly_array_adaptor<double> storage_t;
    typedef ub::matrix<double, ub::column_major, storage_t> matrix_t;

    storage_t storage(rows*columns, mxGetPr(X));

    matrix_t mat_data(rows, columns, storage);


    // this is the form of the data extracted from the storage
    typedef ub::vector<double> data_t;
    typedef robust_pca_impl< data_t > robust_pca_t;

    // TODO put the dimension
    // for the library to work, we need to allocate some temporary storage
    // we also allocate the output if given
    plhs[0] = mxCreateDoubleMatrix(dimension, dimension, mxREAL);
    mxArray *outputMatrix = plhs[0];
    assert(outputMatrix);


    storage_t storageOutput(dimension * dimension, mxGetPr(outputMatrix));

    matrix_t mat_eigen_vectors(dimension, dimension , storageOutput);


    typedef row_iter<matrix_t> const_row_iter_t;

    


    // should be matlab style
    std::vector<data_t> temporary_data(rows);
    std::vector<data_t> eigenvectors;

    robust_pca_t instance;
    instance.batch_process(
      100,
      dimension,
      const_row_iter_t(mat_data, 0),
      const_row_iter_t(mat_data, mat_data.size1()),
      temporary_data.begin(),
      eigenvectors);

  }





  //typedef robust_pca::robust_pca_impl< boost::numeric::ublas::vector<double> > robust_pca_t;

  //robust_pca_t instance;



}