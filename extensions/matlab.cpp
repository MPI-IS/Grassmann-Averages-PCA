// Mex wrapper file for robust PCA


#include "mex.h"

#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <include/robust_pca.hpp>


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

template<class T>
class readonly_array_adaptor : public boost::numeric::ublas::storage_array<readonly_array_adaptor<T> >
{

  typedef readonly_array_adaptor<T> self_type;

public:
  typedef size_t size_type;
  typedef std::ptrdiff_t difference_type;
  typedef T value_type;
  typedef const T &const_reference;
  typedef const T *const_pointer;


public:

  // Construction and destruction
  readonly_array_adaptor() :
    size_(0), data_(0)
  {}

  readonly_array_adaptor(size_type size, const_pointer data) :
    size_(size), data_(data) {
  }

  ~readonly_array_adaptor() {}

  readonly_array_adaptor(const self_type& rhs) : size_(rhs.size_), data_(rhs.data_)
  {}

  // Resizing
  void resize(size_type size)
  {
    size_ = size;
  }

  void resize(size_type size, const_pointer data)
  {
    size_ = size;
    data_ = data;
  }

  // Random Access Container
  size_type max_size() const
  {
    return std::numeric_limits<size_type>::max();
  }

  bool empty() const
  {
    return size_ == 0;
  }

  size_type size() const
  {
    return size_;
  }

  // Element access
  const_reference operator [] (size_type i) const
  {
    assert(i < size_);
    return data_[i];
  }

  // Iterators simply are pointers.
  typedef const_pointer const_iterator;

  const_iterator begin() const
  {
    return data_;
  }
  const_iterator end() const
  {
    return data_ + size_;
  }

  // this typedef is used by vector and matrix classes
  typedef const_pointer iterator;

  // Reverse iterators
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
  typedef std::reverse_iterator<iterator> reverse_iterator;

  const_reverse_iterator rbegin() const
  {
    return const_reverse_iterator(end());
  }

  const_reverse_iterator rend() const
  {
    return const_reverse_iterator(begin());
  }

private:
  size_type size_;
  const_pointer data_;
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

#if 0

//! An helper class for iterating over rows of a matrix created in matlab
//! To be adapted
template <class matrix_value_t>
class row_iter_pointer_matrix_column_major
  : public boost::iterator_facade<
        row_iter_pointer_matrix_column_major<matrix_value_t>              // Derived
      , boost::numeric::ublas::matrix_row<matrix_value_t>   // Value
      , boost::random_access_traversal_tag            // CategoryOrTraversal
      , boost::numeric::ublas::matrix_row<matrix_value_t>   // reference
    >
{
private:
  typedef row_iter_pointer_matrix_column_major<matrix_value_t> this_type;

  struct enabler {};

  size_t index;
  size_t nb_rows;
  size_t nb_columns;

  
  matrix_value_t *p_matrix;



  typedef boost::numeric::ublas::matrix_row<matrix_value_t> return_t;
  static const size_t max_size;

public:


  row_iter_pointer_matrix_column_major() : 
    index(max_size), nb_rows(0), nb_columns(0)
    p_matrix(0)
  {}


  //! @param index the current row position
  row_iter_pointer_matrix_column_major(
    matrix_value_t *mat, 
    size_t nb_rows_, 
    size_t nb_columns_, 
    size_t index_) : 
    index(index_), nb_rows(nb_rows_), nb_columns(nb_columns_) 
    p_matrix(mat)
  {}


  template <class other_matrix_value_t>
  row_iter_pointer_matrix_column_major(
    row_iter_pointer_matrix_column_major<other_matrix_value_t> const& other
    , typename boost::enable_if<
        boost::is_convertible<other_matrix_value_t, matrix_value_t>, enabler>::type = enabler())
    : index(other.index), nb_rows(other.nb_rows), nb_columns(other.nb_columns), 
      p_matrix(other.p_matrix)
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
    assert(p_matrix == other.p_matrix);
    return this->index == other.index;
  }

  return_t dereference() const
  {
    assert(p_matrix);
    return return_t(*p_matrix, index);
  }

  typename this_type::difference_type distance_to(this_type const& r) const
  {
    assert((p_matrix != 0) && (r.p_matrix == p_matrix));
    return typename this_type::difference_type(r.index) - typename this_type::difference_type(index); // sign promotion
  }

};
template <class matrix_value_t>
const size_t row_iter_pointer_matrix_column_major<matrix_value_t>::max_size = std::numeric_limits<size_t>::max();

#endif




// copy of the class from the unit tests, need to clean a bit
template <class matrix_t>
class row_iter
  : public boost::iterator_facade<
      row_iter<matrix_t>                            // Derived
    , boost::numeric::ublas::matrix_row<matrix_t>   // Value
    , boost::random_access_traversal_tag            // CategoryOrTraversal
    , boost::numeric::ublas::matrix_row<matrix_t>   // reference
  >
{
private:
  typedef row_iter<matrix_t> this_type;

  struct enabler {};

  size_t index;
  matrix_t *matrix;

  typedef boost::numeric::ublas::matrix_row<matrix_t> return_t;
  static const size_t max_size;

public:


  row_iter() : index(max_size), matrix(0)
  {}

  row_iter(matrix_t &mat, size_t index_) : index(index_), matrix(&mat)
  {}


  template <class other_matrix_t>
  row_iter(
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
const size_t row_iter<matrix_t>::max_size = std::numeric_limits<size_t>::max();






void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  namespace ub = boost::numeric::ublas;
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
    typedef boost::numeric::ublas::vector<double> data_t;
    typedef robust_pca::robust_pca_impl< data_t > robust_pca_t;

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