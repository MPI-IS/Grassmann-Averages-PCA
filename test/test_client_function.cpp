#include <boost/test/unit_test.hpp>

#include <include/robust_pca.hpp>

// data stored into a matrix
#include <boost/numeric/ublas/matrix.hpp>

// generating data randomly
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include <boost/iterator/transform_iterator.hpp>


// does not work, the operator() receives the it->operator* and it becomes a bit messy
template <class matrix_t>
struct get_row_proxy
{
  typedef boost::numeric::ublas::matrix_row<matrix_t> matrix_row_t;

  typedef matrix_row_t result_type;
  

  matrix_t const& mat;

  get_row_proxy(matrix_t const& mat_) : mat(mat_) {}

  template <class row_iterator_t>
  matrix_row_t operator()(row_iterator_t const& it) const
  {
    return matrix_row_t(mat, it.index());
  }
};


// this does not work because the iterator on matrice is not fully indexed (iterator version of the class)
// or because int/size_t does not have an iterator semantic.
#if 0
template <class matrix_t>
class row_iter
  : public boost::iterator_adaptor<
        row_iter<matrix_t>                            // Derived
      , size_t                                        // Base
      , boost::numeric::ublas::matrix_row<matrix_t>   // Value
      , boost::random_access_traversal_tag            // CategoryOrTraversal
      , boost::numeric::ublas::matrix_row<matrix_t>   // reference
    >
{
private:
  struct enabler {};
  matrix_t *matrix;
  typedef row_iter<matrix_t> this_type;
  typedef boost::numeric::ublas::matrix_row<matrix_t> return_t;
  static const size_t max_size;

public:


  row_iter() : row_iter::iterator_adaptor_(), matrix(0)
  {}

  row_iter(matrix_t &mat, typename this_type::base_type it) : row_iter::iterator_adaptor_(it), matrix(&mat)
  {}


  template <class other_matrix_t>
  row_iter(
      row_iter<other_matrix_t> const& other
    , typename boost::enable_if<
          boost::is_convertible<typename other_matrix_t::iterator1, typename matrix_t::iterator1>
        , enabler
      >::type = enabler())
    : row_iter::iterator_adaptor_(other.base()), matrix(other.matrix) 
  {}

private:
  friend class boost::iterator_core_access;
  typename iterator_adaptor::reference dereference() const
  {
    assert(matrix);
    return typename iterator_adaptor::reference(*matrix, this->base_reference());
  }
};

template <class matrix_t>
const size_t row_iter<matrix_t>::max_size = std::numeric_limits<size_t>::max();
#endif

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







boost::random::mt19937 rng;

BOOST_AUTO_TEST_CASE(instance_test)
{
  namespace bubla = boost::numeric::ublas;
  using namespace robust_pca;
  using namespace boost::numeric::ublas;
  using namespace boost::random;

  uniform_real_distribution<double> dist(-1000, 1000);
  

  typedef matrix<double> matrix_t;
  typedef robust_pca_impl<bubla::vector<double> > robust_pca_t;
  

  robust_pca_t instance;

  const int nb_elements = 1000;
  const int dimensions  = 5;

  // creating some data, 1000 lines of a vector of length 5
  matrix_t mat_data(nb_elements, dimensions, 0);

  for(int i = 0; i < nb_elements; i++)
  {
    for(int j = 0; j < dimensions; j++)
    {
      mat_data(i, j) = dist(rng);
    }
  }


  std::vector<double> norms(nb_elements);


  typedef get_row_proxy<matrix_t> row_proxy_tranform_t;


  //BOOST_CHECK(instance.batch_process(
  //              boost::make_transform_iterator(mat_data.begin1(), row_proxy_tranform_t(mat_data)), 
  //              boost::make_transform_iterator(mat_data.end1(), row_proxy_tranform_t(mat_data)), 
  //              norms.begin()));

  typedef row_iter<const matrix_t> const_raw_iter_t;

  BOOST_CHECK(!instance.batch_process(
    const_raw_iter_t(mat_data, 2),
    const_raw_iter_t(mat_data, 0),
    norms.begin()));

  BOOST_CHECK(!instance.batch_process(
    const_raw_iter_t(mat_data, 0),
    const_raw_iter_t(mat_data, 0),
    norms.begin()));


  
  BOOST_CHECK(instance.batch_process(
    const_raw_iter_t(mat_data, 0),
    const_raw_iter_t(mat_data, mat_data.size1()),
    norms.begin()));

}



boost::unit_test::test_suite* init_unit_test_suite( int argc, char* argv[] )
{
  return 0;
}