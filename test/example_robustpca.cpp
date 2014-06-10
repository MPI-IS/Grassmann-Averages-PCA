
#include <boost/random/mersenne_twister.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random/uniform_real_distribution.hpp>



// generating a matrix of random numbers, each row will contain
// a data vector.
void generate_matrix(const int nb_elements, const int dimension matrix_t& mat)
{
  static boost::random::mt19937 rng;
  boost::random::uniform_real_distribution<double> dist(-1000, 1000);
  
  mat_data.resize(nb_elements, dimensions);
  for(int i = 0; i < nb_elements; i++)
  {
    for(int j = 0; j < dimensions; j++)
    {
      mat_data(i, j) = dist(rng);
    }
  }
}

void example_robust_pca_impl()
{
  using namespace robust_pca;
  using namespace robust_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;

  
  typedef ub::vector<double> data_t;                        // type of the vectors 
  typedef ub::matrix<double> matrix_t;                      // type of the structure holding the data

  typedef robust_pca_impl< data_t > robust_pca_t;           // the type of the structure for the computation of the robust pca
                                                             // data_t tells which kind of structure will be used for internal computations
                                                             // and will be received from the data iterators.
  
  typedef row_iter<const matrix_t> const_row_iter_t;        // iterator on the data, deferencing one
                                                             // of these iterator should yield a structure
                                                             // convertible to data_t

  
  const int dimensions = 5;                   // each vector is of dimension 5
  const int max_dimension_to_compute = 3;     // we want only the first 3 eigen-vectors
  const int max_iterations = 1000;            // at most 1000 iterations before giving up for the current dimension


  // configure the first points
  const double initial_point[] = {0.2097, 0.3959, 0.5626, 0.2334, 0.6545}; // some dummy initial point
  std::vector<data_t> eigen_vectors(dimensions);
  ub::vector<double> vec_initial_point(initial_point, initial_point + dimensions);
  std::vector< ub::vector<double> > v_init_points(max_dimension_to_compute, vec_initial_point);

  
  // the instance
  robust_pca_t instance;
  
  // some configuration
  instance.set_nb_processors(4);              // using 4 processors
  instance.set_nb_processors(4);              // using 4 processors
  
  bool return_calue = instance.batch_process(
    max_iterations,
    max_dimension_to_compute,
    const_row_iter_t(mat_data, 0),
    const_row_iter_t(mat_data, mat_data.size1()),
    eigen_vectors.begin(),
    &v_init_points));


}