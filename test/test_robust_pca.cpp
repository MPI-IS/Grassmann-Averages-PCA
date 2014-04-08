// Copyright 2014, Max Planck Institute for Intelligent Systems.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

/*!@file
 * This file includes the tests for the robust pca
 */


#include <boost/test/unit_test.hpp>

#include <test/test_main.hpp>

#include <include/robust_pca.hpp>
#include <include/private/boost_ublas_matrix_helper.hpp>

// data stored into a matrix
#include <boost/numeric/ublas/io.hpp>




#include <fstream>






BOOST_FIXTURE_TEST_SUITE(robust_pca, fixture_simple_matrix_creation)

BOOST_AUTO_TEST_CASE(returns_false_for_inapropriate_inputs)
{
  using namespace robust_pca;
  using namespace robust_pca::ublas_adaptor;
  namespace ub = boost::numeric::ublas;

  typedef robust_pca_impl< ub::vector<double> > robust_pca_t;

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

  typedef robust_pca_impl< ub::vector<double> > robust_pca_t;  
  robust_pca_t instance;
  typedef row_iter<const matrix_t> const_row_iter_t;
  
  typedef ub::vector<double> data_t;


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

  std::vector< ub::vector<double> > v_init_points(dimensions, vec_initial_point);

  BOOST_CHECK(instance.batch_process(
    max_iterations,
    dimensions,
    const_row_iter_t(mat_data, 0),
    const_row_iter_t(mat_data, mat_data.size1()),
    temporary_data.begin(),
    eigen_vectors.begin(),
    &v_init_points));



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
  using namespace robust_pca;
  using namespace robust_pca::ublas_adaptor;
  namespace ub = boost::numeric::ublas;


  typedef robust_pca_impl< ub::vector<double> > robust_pca_t;  
  robust_pca_t instance;
  typedef row_iter<const matrix_t> const_row_iter_t;
  
  typedef boost::numeric::ublas::vector<double> data_t;

  std::vector<data_t> temporary_data(nb_elements);
  std::vector<data_t> eigen_vectors(dimensions);
  const int max_iterations = 1000;


  BOOST_CHECK(instance.batch_process(
    max_iterations,
    dimensions,
    const_row_iter_t(mat_data, 0),
    const_row_iter_t(mat_data, mat_data.size1()),
    temporary_data.begin(),
    eigen_vectors.begin()));


  BOOST_MESSAGE("Generated eigen vectors are:");

  for(int i = 0; i < dimensions; i++)
  {
    BOOST_MESSAGE("vector " << i << " :" << eigen_vectors[i]);
  }
}
#endif

BOOST_AUTO_TEST_SUITE_END();

