// Mex wrapper file for robust PCA


#include "mex.h"

#include <boost/iterator/iterator_adaptor.hpp>

//! An helper class for iterating over rows of a matrix created in matlab
//! To be adapted
template <class matrix_t>
class row_iter_matlab_matrix
  : public boost::iterator_facade<
        row_iter_matlab_matrix<matrix_t>              // Derived
      , boost::numeric::ublas::matrix_row<matrix_t>   // Value
      , boost::random_access_traversal_tag            // CategoryOrTraversal
      , boost::numeric::ublas::matrix_row<matrix_t>   // reference
    >
{
private:
  typedef row_iter_matlab_matrix<matrix_t> this_type;

  struct enabler {};

  size_t index;
  matrix_t *matrix;

  typedef boost::numeric::ublas::matrix_row<matrix_t> return_t;
  static const size_t max_size;

public:


  row_iter_matlab_matrix() : index(max_size), matrix(0)
  {}

  row_iter_matlab_matrix(matrix_t &mat, size_t index_) : index(index_), matrix(&mat)
  {}


  template <class other_matrix_t>
  row_iter_matlab_matrix(
      row_iter<other_matrix_t> const& other
    , typename boost::enable_if<
          boost::is_convertible<typename other_matrix_t::iterator1, typename matrix_t::iterator1>
        , enabler
      >::type = enabler())
    : index(other.index), matrix(other.matrix) 
  {}

private:
  friend class boost::iterator_core_access;
  
  void increment()
  { 
    assert(index < max_size);
    index++; 
  }

  bool equal(this_type const& other) const
  {
    assert(matrix == other.matrix);
    return this->index == other.index;
  }

  return_t dereference() const
  {
    assert(matrix);
    return return_t(*matrix, index);
  }

  typename this_type::difference_type distance_to(this_type const& r) const
  {
    assert((matrix != 0) && (r.matrix == matrix));
    return typename this_type::difference_type(r.index) - typename this_type::difference_type(index); // sign promotion
  }

};
template <class matrix_t>
const size_t row_iter_matlab_matrix<matrix_t>::max_size = std::numeric_limits<size_t>::max();


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

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





  // for the library to work, we need to allocate some temporary storage
  // we also allocate the output if given



  plhs[0] = mxCreateDoubleMatrix(dimension , dimension , mxREAL); 




  typedef robust_pca::robust_pca_impl< boost::numeric::ublas::vector<double> > robust_pca_t;

  robust_pca_t instance;



}