// Copyright 2014, Max Planck Institute for Intelligent Systems.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

/*!@file
 * This file includes the tests for the trimmed version of the robust pca
 */

#include <boost/test/unit_test.hpp>
#include <test/test_main.hpp>

#include <include/robust_pca.hpp>
#include <include/private/boost_ublas_matrix_helper.hpp>

// data stored into a matrix
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

// generating data randomly
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <fstream>



BOOST_FIXTURE_TEST_SUITE(robust_pca_trimmed, fixture_simple_matrix_creation)



BOOST_AUTO_TEST_CASE(returns_false_for_inapropriate_inputs)
{
  using namespace robust_pca;
  using namespace robust_pca::ublas_adaptor;
  namespace ub = boost::numeric::ublas;

  typedef robust_pca_with_trimming_impl< ub::vector<double> > robust_pca_t;

  robust_pca_t instance;


  typedef row_iter<const matrix_t> const_row_iter_t;
  typedef ub::vector<double> data_t;

  std::vector<data_t> temporary_data(nb_elements);
  std::vector<data_t> eigen_vectors(dimensions);
  const int max_iterations = 1000;

  BOOST_CHECK(!instance.batch_process(
    max_iterations,
    dimensions,
    const_row_iter_t(mat_data, 2),
    const_row_iter_t(mat_data, 0),
    temporary_data.begin(),
    eigen_vectors.begin()));

  BOOST_CHECK(!instance.batch_process(
    max_iterations,
    dimensions,
    const_row_iter_t(mat_data, 2),
    const_row_iter_t(mat_data, 0),
    temporary_data.begin(),
    eigen_vectors.begin()));
}


BOOST_AUTO_TEST_CASE(smoke_and_orthogonality_tests)
{
  using namespace robust_pca;
  using namespace robust_pca::ublas_adaptor;
  namespace ub = boost::numeric::ublas;

  typedef robust_pca_with_trimming_impl< ub::vector<double> > robust_pca_t;
  robust_pca_t instance;
  typedef row_iter<const matrix_t> const_row_iter_t;

  typedef ub::vector<double> data_t;


  std::vector<data_t> temporary_data(nb_elements);
  std::vector<data_t> eigen_vectors(dimensions);
  const int max_iterations = 1000;

  //const double initial_point[] = { 0.2097, 0.3959, 0.5626, 0.2334, 0.6545 };
  
  // this is the initialisation of the sequence of random vectors for each dimension 
  // and some gram_schmidt orthogonalisation was also applied on it.
  const double initial_point[] = {
    0.1843, -0.5685, -0.4177,  0.5001, -0.4672,
    0.9318, -0.0753,  0.0646, -0.3479,  0.0291,
   -0.0042,  0.5226, -0.0230, -0.1980, -0.8289,
    0.1709,  0.5421, -0.7280,  0.2310,  0.3059,
   -0.2617, -0.3227, -0.5393, -0.7324, -0.0122,
  };
  //BOOST_REQUIRE_EQUAL(dimensions, sizeof(initial_point) / sizeof(initial_point[0])); // just in case


  std::vector< ub::vector<double> > vec_initial_point(dimensions);
  for(int i = 0; i < dimensions; i++)
  {
    vec_initial_point[i].resize(dimensions);
    for(int j = 0; j < dimensions; j++)
    {
      vec_initial_point[i](j) = initial_point[j];
    }
  }

  BOOST_CHECK(instance.batch_process(
    max_iterations,
    dimensions,
    const_row_iter_t(mat_data, 0),
    const_row_iter_t(mat_data, mat_data.size1()),
    temporary_data.begin(),
    eigen_vectors.begin(),
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
  for(int i = 0; i < dimensions - 1; i++)
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
   -0.0355,  0.1244, -0.0150,  0.9881,  0.0816,
    0.9664, -0.0598, -0.2245,  0.0470, -0.0992,
   -0.0607,  0.1019,  0.1580,  0.0682, -0.9779,
    0.0755,  0.9848, -0.0033, -0.1287,  0.0884,
   -0.2353,  0.0254, -0.9615, -0.0148, -0.1391,
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





BOOST_AUTO_TEST_SUITE_END();
