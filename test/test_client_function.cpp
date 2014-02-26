#include <boost/test/unit_test.hpp>

#include <include/robust_pca.hpp>

// data stored into a matrix
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

// generating data randomly
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include <boost/iterator/transform_iterator.hpp>

#include <fstream>

// the random number generator
boost::random::mt19937 rng;


//! An helper class for iterating over rows of a matrix
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




// Fixture for the tests
struct fixture_simple_matrix_creation
{
  typedef boost::numeric::ublas::matrix<double> matrix_t;


  static const int nb_elements = 1000;
  static const int dimensions  = 5;
  matrix_t mat_data;
  //std::vector<double> norms;

  boost::random::uniform_real_distribution<double> dist;


  fixture_simple_matrix_creation() : dist(-1000, 1000) 
  {
    // creating some data, 1000 lines of a vector of length 5
    mat_data.resize(nb_elements, dimensions);
    //norms.resize(nb_elements);

    // default seed for reproductible sequences
    rng.seed();

    //std::cout << "current seed : " << ;

    //const std::string filename = "./toto.txt";
    //std::ofstream ff(filename.c_str());

    //BOOST_REQUIRE(ff.is_open());

    for(int i = 0; i < nb_elements; i++)
    {
      for(int j = 0; j < dimensions; j++)
      {
        mat_data(i, j) = dist(rng);
        //ff << mat_data(i, j) << " ";
      }
      //ff << std::endl;
    }

    
    

  }
};


BOOST_FIXTURE_TEST_SUITE(basic_checks, fixture_simple_matrix_creation)

BOOST_AUTO_TEST_CASE(returns_false_for_inapropriate_inputs)
{
  typedef robust_pca::robust_pca_impl<boost::numeric::ublas::vector<double> > robust_pca_t;

  robust_pca_t instance;


  typedef row_iter<const matrix_t> const_row_iter_t;
  typedef boost::numeric::ublas::vector<double> data_t;

  std::vector<data_t> temporary_data(nb_elements);
  std::vector<data_t> eigen_vectors(dimensions);
  const int max_iterations = 1000;

  BOOST_CHECK(!instance.batch_process(
    max_iterations,
    dimensions,
    const_row_iter_t(mat_data, 2),
    const_row_iter_t(mat_data, 0),
    temporary_data.begin(),
    //norms.begin(),
    eigen_vectors));

  BOOST_CHECK(!instance.batch_process(
    max_iterations,
    dimensions,
    const_row_iter_t(mat_data, 2),
    const_row_iter_t(mat_data, 0),
    temporary_data.begin(),
    //norms.begin(),
    eigen_vectors));
}


BOOST_AUTO_TEST_CASE(smoke_and_orthogonality_tests)
{
  namespace ub = boost::numeric::ublas;
  typedef robust_pca::robust_pca_impl< ub::vector<double> > robust_pca_t;  
  robust_pca_t instance;
  typedef row_iter<const matrix_t> const_row_iter_t;
  
  typedef boost::numeric::ublas::vector<double> data_t;


  std::vector<data_t> temporary_data(nb_elements);
  std::vector<data_t> eigen_vectors(dimensions);
  const int max_iterations = 1000;

  const double initial_point[] = {0.2097, 0.3959, 0.5626, 0.2334, 0.6545};
  BOOST_REQUIRE_EQUAL(dimensions, sizeof(initial_point)/sizeof(initial_point[0])); // just in case

  ub::vector<double> vec_initial_point(dimensions);
  for(int i = 0; i < dimensions; i++)
  {
    vec_initial_point(i) = initial_point[i];
  }

  BOOST_CHECK(instance.batch_process(
    max_iterations,
    dimensions,
    const_row_iter_t(mat_data, 0),
    const_row_iter_t(mat_data, mat_data.size1()),
    temporary_data.begin(),
    //norms.begin(),
    eigen_vectors,
    &vec_initial_point));



  // testing the output sizes
  BOOST_REQUIRE_EQUAL(eigen_vectors.size(), dimensions);
  for(int i = 0; i < dimensions; i++)
  {
    BOOST_CHECKPOINT("testing eigenvector size for vector " << i);
    BOOST_REQUIRE_EQUAL(eigen_vectors[i].size(), dimensions);
  }


  BOOST_MESSAGE("Generated eigen vectors are:");

  for(int i = 0; i < dimensions; i++)
  {
    BOOST_MESSAGE("vector " << i << " :" << eigen_vectors[i]);
  }

  // testing orthogonality of all eigenvectors
  for(int i = 0; i < dimensions-1; i++)
  {
    for(int j = i + 1; j < dimensions; j++)
    {
      BOOST_CHECK_LE(ub::inner_prod(eigen_vectors[i], eigen_vectors[j]), 1E-6);
    }
  }

  for(int i = 0; i < dimensions; i++)
  {
    BOOST_CHECK_CLOSE(ub::inner_prod(eigen_vectors[i], eigen_vectors[i]), 1, 1E-6);
  }


  // testing against the matlab output, the script being given by Sorent Hauberg, and the init between
  // each dimension iteration being given by the vector "initial_point" above. Each column represents 
  // an eigenvector.
  static const double matlab_data[] = {
   -0.0219,    0.0905,   -0.0057,    0.0914,    0.9914,
   -0.0512,   -0.0875,    0.9873,   -0.1202,    0.0236,
   -0.0624,    0.9889,    0.0908,    0.0336,   -0.0942,
    0.0280,   -0.0508,    0.1193,    0.9875,   -0.0851,
    0.9961,    0.0609,    0.0529,   -0.0298,    0.0195
  };


  for(int i = 0; i < dimensions; i++)
  {
    ub::vector<double> current_matlab_vector(dimensions);
    for(int j = 0; j < dimensions; j++)
    {
      current_matlab_vector(j) = matlab_data[i + j*dimensions];
    }
    BOOST_CHECKPOINT("iteration " << i);
    BOOST_CHECK_LE(ub::norm_2(eigen_vectors[i] - current_matlab_vector), 1E-3);
  }
  



}

#if 0
BOOST_AUTO_TEST_CASE(checking_against_matlab)
{
  typedef robust_pca::robust_pca_impl< boost::numeric::ublas::vector<double> > robust_pca_t;  
  robust_pca_t instance;
  typedef row_iter<const matrix_t> const_row_iter_t;
  
  typedef boost::numeric::ublas::vector<double> data_t;

  std::vector<data_t> temporary_data(nb_elements);
  std::vector<data_t> eigen_vectors(dimensions);
  const int max_iterations = 1000;


  BOOST_CHECK(instance.batch_process(
    max_iterations,
    const_row_iter_t(mat_data, 0),
    const_row_iter_t(mat_data, mat_data.size1()),
    temporary_data.begin(),
    norms.begin(),
    eigen_vectors));


  BOOST_MESSAGE(
    "Generated eigen vectors are:");

  for(int i = 0; i < dimensions; i++)
  {
    BOOST_MESSAGE("vector " << i << " :" << eigen_vectors[i]);
  }
}
#endif

BOOST_AUTO_TEST_SUITE_END();

boost::unit_test::test_suite* init_unit_test_suite( int argc, char* argv[] )
{
  return 0;
}