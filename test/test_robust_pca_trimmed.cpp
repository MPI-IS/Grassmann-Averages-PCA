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
  robust_pca_t instance(0., 1);
  typedef row_iter<const matrix_t> const_row_iter_t;

  typedef ub::vector<double> data_t;


  std::vector<data_t> temporary_data(nb_elements);
  std::vector<data_t> eigen_vectors(dimensions);
  const int max_iterations = 1000;

  //const double initial_point[] = { 0.2097, 0.3959, 0.5626, 0.2334, 0.6545 };
  
  // this is the initialisation of the sequence of random vectors for each dimension 
  // and some gram_schmidt orthogonalisation was also applied on it.
  const double initial_point[] = {
    0.2658, -0.4880, 0.4029, 0.4855, 0.5414,
    0.8306, 0.3194, 0.2228, -0.3925, 0.0663,
    -0.2066, -0.4473, 0.6346, -0.4964, -0.3288,
    -0.3310, 0.0890, -0.0642, -0.5345, 0.7699,
    -0.2953, 0.6722, 0.6174, 0.2794, 0.0408,
  };
  //BOOST_REQUIRE_EQUAL(dimensions, sizeof(initial_point) / sizeof(initial_point[0])); // just in case


  std::vector< ub::vector<double> > vec_initial_point(dimensions);
  for(int i = 0; i < dimensions; i++)
  {
    vec_initial_point[i].resize(dimensions);
    for(int j = 0; j < dimensions; j++)
    {
      vec_initial_point[i](j) = initial_point[i*dimensions + j];
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


  // testing unitarity
  for(int i = 0; i < dimensions; i++)
  {
    BOOST_CHECK_CLOSE(ub::inner_prod(eigen_vectors[i], eigen_vectors[i]), 1, 1E-6);
  }


  // testing against the matlab output, the script being given by Sorent Hauberg, and the init between
  // each dimension iteration being given by the vector "initial_point" above. Each column represents 
  // an eigenvector.
  static const double matlab_data[] = {
   -0.0506,  0.9918, -0.0085,  0.1167, -0.0055,
   -0.1120, -0.0348, -0.4900,  0.2504,  0.8267,
   -0.0555,  0.0622,  0.6842, -0.4771,  0.5452,
    0.0689,  0.0965, -0.5397, -0.8317, -0.0546,
    0.9885,  0.0436,  0.0200,  0.0656,  0.1278,
  };


  for(int i = 0; i < dimensions; i++)
  {
    ub::vector<double> current_matlab_vector(dimensions);
    for(int j = 0; j < dimensions; j++)
    {
      current_matlab_vector(j) = matlab_data[i + j*dimensions];
    }
    BOOST_CHECKPOINT("iteration " << i);
    BOOST_CHECK_LE(ub::norm_2(eigen_vectors[i] - current_matlab_vector), 2.6E-2); 
    // there is a slight difference between the two implementation when enabling the P^2 algorithm
    // for producing the quantiles.
  }




}





BOOST_AUTO_TEST_SUITE_END();
