#include <boost/test/unit_test.hpp>

#include <include/robust_pca.hpp>

// data stored into a matrix
#include <boost/numeric/ublas/matrix.hpp>

// generating data randomly
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>


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

  BOOST_CHECK(instance.batch_process(mat_data.begin1(), mat_data.end1(), norms.begin()));
}



boost::unit_test::test_suite* init_unit_test_suite( int argc, char* argv[] )
{
  return 0;
}