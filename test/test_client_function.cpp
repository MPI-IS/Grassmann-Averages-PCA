#include <boost/test/unit_test.hpp>

#include <include/robust_pca.hpp>
#include <boost/numeric/ublas/matrix.hpp>


BOOST_AUTO_TEST_CASE(instance_test)
{
  using namespace robust_pca;

  typedef boost::numeric::ublas::matrix<double> matrix_t;
  typedef robust_pca_impl<double> robust_pca_t;
  

  robust_pca_t instance;


  // creating some data, 1000 lines of a vector of length 5
  matrix_t mat_data(1000, 5, 0);

  BOOST_CHECK(instance.batch_process(mat_data.begin1(), mat_data.end1()));
}



boost::unit_test::test_suite* init_unit_test_suite( int argc, char* argv[] )
{

}