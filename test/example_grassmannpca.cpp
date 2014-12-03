
#include <boost/random/mersenne_twister.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random/uniform_real_distribution.hpp>



#include <include/grassmann_pca.hpp>
#include <include/private/boost_ublas_row_iterator.hpp>



// generating a matrix of random numbers, each row will contain
// a data vector.
template <class matrix_t>
void generate_matrix(const int nb_elements, const int dimension, matrix_t& mat_data)
{
  static boost::random::mt19937 rng;
  boost::random::uniform_real_distribution<double> dist(-1000, 1000);
  
  mat_data.resize(nb_elements, dimension);
  for(int i = 0; i < nb_elements; i++)
  {
    for(int j = 0; j < dimension; j++)
    {
      mat_data(i, j) = dist(rng);
    }
  }
}

void example_grassmann_pca_impl()
{
  using namespace grassmann_averages_pca;
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;

  
  typedef ub::vector<double> data_t;                         // type of the vectors 
  typedef ub::matrix<double> matrix_t;                       // type of the structure holding the data

  typedef grassmann_pca< data_t > grassmann_pca_t;           // the type of the structure for the computation of the grassmann averages pca
                                                             // data_t tells which kind of structure will be used for internal computations
                                                             // and will be received from the data iterators.
  
  typedef row_iter<const matrix_t> const_row_iter_t;         // iterator on the data, deferencing one
                                                             // of these iterator should yield a structure
                                                             // convertible to data_t. This iterator iterates 
                                                             // over the rows of a given matrix.

  
  const int dimensions = 5;                   // each vector is of dimension 5
  const int max_dimension_to_compute = 3;     // we want only the first 3 basis-vectors
  const int max_iterations = 1000;            // at most 1000 iterations before giving up for the current dimension
  
  // generating the data points
  matrix_t mat_data;
  generate_matrix(10000, dimensions, mat_data);


  // configure the first points
  const double initial_point[] = {0.2097, 0.3959, 0.5626, 0.2334, 0.6545}; // some dummy initial point
  data_t vec_initial_point(dimensions);
  std::copy(initial_point, initial_point + dimensions, vec_initial_point.begin());
  std::vector< data_t > v_init_points(max_dimension_to_compute, vec_initial_point);

  // allocate the output
  std::vector<data_t> basis_vectors(max_dimension_to_compute);
  
  // the instance
  grassmann_pca_t instance;
  
  // some configuration
  if(!instance.set_nb_processors(4))              // using 4 processors
  {
    std::cerr << "Error while setting the number of threads" << std::endl;
    return;
  }
  
  if(!instance.set_max_chunk_size(1000))          // using max chunk size of 1000
  {
    std::cerr << "Error while setting the chunk size" << std::endl;
    return;
  }
  
  if(!instance.set_nb_steps_pca(3))               // using 3 PCA steps
  {
    std::cerr << "Error while setting the number of PCA steps" << std::endl;
    return;
  }
  
  
  bool return_calue = instance.batch_process(
    max_iterations,
    max_dimension_to_compute,
    const_row_iter_t(mat_data, 0),
    const_row_iter_t(mat_data, mat_data.size1()),
    basis_vectors.begin(),
    &v_init_points);

  if(!return_calue)
  {
    std::cerr << "Error during the computation of the GrassmannAveragesPCA" << std::endl;  
  }
}
