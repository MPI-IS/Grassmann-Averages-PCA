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


template <class matrix_t>
class row_iter
  : public boost::iterator_adaptor<
        row_iter<matrix_t>              // Derived
      , row_iter<matrix_t>              // Base
      , boost::numeric::ublas::matrix_row<matrix_t>              // Value
      , boost::random_access_traversal_tag // CategoryOrTraversal
    >
{

 public:
    row_iter()
      : row_iter::iterator_adaptor_(0) {}

    explicit row_iter(Value* p)
      : node_iter::iterator_adaptor_(p) {}

    template <class OtherValue>
    node_iter(
        node_iter<OtherValue> const& other
      , typename boost::enable_if<
            boost::is_convertible<OtherValue*,Value*>
          , enabler
        >::type = enabler()
    )
      : node_iter::iterator_adaptor_(other.base()) {}

 private:
    friend class boost::iterator_core_access;
    void increment() { this->base_reference() = this->base()->next(); }
};



boost::random::mt19937 rng;

BOOST_AUTO_TEST_CASE(instance_test)
{
  using namespace robust_pca;
  using namespace boost::numeric::ublas;
  using namespace boost::random;

  uniform_real_distribution<double> dist(-1000, 1000);
  

  typedef matrix<double> matrix_t;
  typedef robust_pca_impl<double> robust_pca_t;
  

  robust_pca_t instance;

  const int nb_elements = 1000;
  const int dimensions  = 5;

  // creating some data, 1000 lines of a vector of length 5
  matrix_t mat_data(nb_elements, dimensions, 0);

  matrix_t::iterator1 it(mat_data.begin1());
  for(int i = 0; i < nb_elements; i++)
  {
    for(int j = 0; j < dimensions; j++, ++it)
    {
      *it = dist(rng);
    }
  }


  std::vector<double> norms(nb_elements);


  typedef get_row_proxy<matrix_t> row_proxy_tranform_t;


  BOOST_CHECK(instance.batch_process(
                boost::make_transform_iterator(mat_data.begin1(), row_proxy_tranform_t(mat_data)), 
                boost::make_transform_iterator(mat_data.end1(), row_proxy_tranform_t(mat_data)), 
                norms.begin()));
}



boost::unit_test::test_suite* init_unit_test_suite( int argc, char* argv[] )
{
  return 0;
}