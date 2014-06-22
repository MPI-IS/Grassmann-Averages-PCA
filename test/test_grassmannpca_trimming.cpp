// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

/*!@file
 * This file includes the tests for the trimmed version of the grassmann pca
 */

#include <boost/test/unit_test.hpp>
#include <test/test_main.hpp>

#include <include/robust_pca_trimming.hpp>
#include <include/private/boost_ublas_row_iterator.hpp>

// data stored into a matrix
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

// generating data randomly
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>


// boost chrono
#include <boost/chrono/include.hpp>

#include <fstream>



BOOST_FIXTURE_TEST_SUITE(grassmann_pca_trimmed, fixture_simple_matrix_creation)



BOOST_AUTO_TEST_CASE(returns_false_for_inapropriate_inputs)
{
  using namespace grassmann_averages_pca;
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;

  typedef grassmann_pca_with_trimming< ub::vector<double> > grassmann_pca_t;

  grassmann_pca_t instance;


  typedef row_iter<const matrix_t> const_row_iter_t;
  typedef ub::vector<double> data_t;

  std::vector<data_t> basis_vectors(dimensions);
  const int max_iterations = 1000;

  BOOST_CHECK(!instance.batch_process(
    max_iterations,
    dimensions,
    const_row_iter_t(mat_data, 2),
    const_row_iter_t(mat_data, 0),
    basis_vectors.begin()));

  BOOST_CHECK(!instance.batch_process(
    max_iterations,
    dimensions,
    const_row_iter_t(mat_data, 2),
    const_row_iter_t(mat_data, 0),
    basis_vectors.begin()));
}


BOOST_AUTO_TEST_CASE(smoke_and_orthogonality_tests)
{
  using namespace grassmann_averages_pca;
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;
  typedef boost::chrono::steady_clock clock_type;
  

  typedef grassmann_pca_with_trimming< ub::vector<double> > grassmann_pca_t;
  grassmann_pca_t instance(0.);
  typedef row_iter<const matrix_t> const_row_iter_t;

  typedef ub::vector<double> data_t;

  BOOST_CHECK(instance.set_nb_processors(1));
  BOOST_CHECK(instance.set_nb_steps_pca(0));


  std::vector<data_t> basis_vectors(dimensions);
  const int max_iterations = 1000;

  clock_type::duration elapsed;
  if(DATA_DIMENSION == 5)
  { 
  
    // this is the initialisation of the sequence of random vectors for each dimension 
    // and some gram_schmidt orthogonalisation was also applied on it.
    const double initial_point[] = {
       0.2658,   -0.4880,    0.4029,    0.4855,    0.5414,
       0.8306,    0.3194,    0.2228,   -0.3925,    0.0663,
      -0.2066,   -0.4473,    0.6346,   -0.4964,   -0.3288,
      -0.3310,    0.0890,   -0.0642,   -0.5345,    0.7699,
      -0.2953,    0.6722,    0.6174,    0.2794,    0.0408
    };
    //BOOST_REQUIRE_EQUAL(dimensions, sizeof(initial_point) / sizeof(initial_point[0])); // just in case


    std::vector< ub::vector<double> > vec_initial_point(dimensions);
    for(int i = 0; i < dimensions; i++)
    {
      vec_initial_point[i].resize(dimensions);
      for(int j = 0; j < dimensions; j++)
      {
        vec_initial_point[i](j) = initial_point[j*dimensions + i];
      }
    }

    clock_type::time_point start = clock_type::now();
    BOOST_CHECK(instance.batch_process(
      max_iterations,
      dimensions,
      const_row_iter_t(mat_data, 0),
      const_row_iter_t(mat_data, mat_data.size1()),
      basis_vectors.begin(),
      &vec_initial_point));
    elapsed = clock_type::now() - start;
  }
  else
  {
    clock_type::time_point start = clock_type::now();
    BOOST_CHECK(instance.batch_process(
      max_iterations,
      dimensions,
      const_row_iter_t(mat_data, 0),
      const_row_iter_t(mat_data, mat_data.size1()),
      basis_vectors.begin()));
    elapsed = clock_type::now() - start;
  }
  
  
  std::cout << "processing " << nb_elements << " elements "
    << "in " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed) << std::endl;
  


  // testing the output sizes
  BOOST_REQUIRE_EQUAL(basis_vectors.size(), dimensions);
  for(int i = 0; i < dimensions; i++)
  {
    BOOST_CHECKPOINT("testing basis vector size for vector " << i);
    BOOST_REQUIRE_EQUAL(basis_vectors[i].size(), dimensions);
  }


  if(DATA_DIMENSION <= 5)
  {
    BOOST_MESSAGE("Generated basis vectors are:");

    for(int i = 0; i < dimensions; i++)
    {
      BOOST_MESSAGE("vector " << i << " :" << basis_vectors[i]);
    }
  }

  // testing orthogonality of all basis vectors
  for(int i = 0; i < dimensions - 1; i++)
  {
    for(int j = i + 1; j < dimensions; j++)
    {
      BOOST_CHECK_LE(ub::inner_prod(basis_vectors[i], basis_vectors[j]), 1E-6);
    }
  }


  // testing unitarity
  for(int i = 0; i < dimensions; i++)
  {
    BOOST_CHECK_CLOSE(ub::inner_prod(basis_vectors[i], basis_vectors[i]), 1, 1E-6);
  }

  if(DATA_DIMENSION == 5)
  {
    // testing against the matlab output, the script being given by Sorent Hauberg, and the init between
    // each dimension iteration being given by the vector "initial_point" above. Each column represents 
    // an basis vector.
    static const double matlab_data[] = {
      -0.0365,   -0.0529,    0.0556,    0.9964,   -0.0058,
       0.9983,   -0.0457,    0.0088,    0.0337,    0.0073,
      -0.0076,   -0.0120,    0.9982,   -0.0565,    0.0165,
      -0.0084,   -0.0242,   -0.0166,    0.0051,    0.9995,
       0.0436,    0.9972,    0.0150,    0.0538,    0.0245,
    };



    for(int i = 0; i < dimensions; i++)
    {
      ub::vector<double> current_matlab_vector(dimensions);
      for(int j = 0; j < dimensions; j++)
      {
        current_matlab_vector(j) = matlab_data[i + j*dimensions];
      }
      BOOST_CHECKPOINT("iteration " << i);
      BOOST_CHECK_LE(ub::norm_2(basis_vectors[i] - current_matlab_vector), 2.6E-2); 
      // there is a slight difference between the two implementation when enabling the P^2 algorithm
      // for producing the quantiles.
    }

  }


}



BOOST_AUTO_TEST_CASE(smoke_and_orthogonality_tests_several_workers)
{

  using namespace grassmann_averages_pca;
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;
  typedef boost::chrono::steady_clock clock_type;
  

  typedef grassmann_pca_with_trimming< ub::vector<double> > grassmann_pca_t;
  grassmann_pca_t instance(0.);
  typedef row_iter<const matrix_t> const_row_iter_t;

  typedef ub::vector<double> data_t;


  std::vector<data_t> basis_vectors(DATA_DIMENSION == 5 ? dimensions : 5);
  const int max_iterations = 1000;

  BOOST_CHECK(instance.set_nb_processors(7)); // each chunk is floor(1000/7) = 142. Last chunk is 148. Just to test sthg different from 10.
  BOOST_CHECK(instance.set_nb_steps_pca(0));

  clock_type::duration elapsed;
  if(DATA_DIMENSION == 5)
  { 
  
    // this is the initialisation of the sequence of random vectors for each dimension 
    // and some gram_schmidt orthogonalisation was also applied on it.
    const double initial_point[] = {
       0.2658,   -0.4880,    0.4029,    0.4855,    0.5414,
       0.8306,    0.3194,    0.2228,   -0.3925,    0.0663,
      -0.2066,   -0.4473,    0.6346,   -0.4964,   -0.3288,
      -0.3310,    0.0890,   -0.0642,   -0.5345,    0.7699,
      -0.2953,    0.6722,    0.6174,    0.2794,    0.0408
    };
    //BOOST_REQUIRE_EQUAL(dimensions, sizeof(initial_point) / sizeof(initial_point[0])); // just in case


    std::vector< ub::vector<double> > vec_initial_point(dimensions);
    for(int i = 0; i < dimensions; i++)
    {
      vec_initial_point[i].resize(dimensions);
      for(int j = 0; j < dimensions; j++)
      {
        vec_initial_point[i](j) = initial_point[j*dimensions + i];
      }
    }

    clock_type::time_point start = clock_type::now();
    BOOST_CHECK(instance.batch_process(
      max_iterations,
      dimensions,
      const_row_iter_t(mat_data, 0),
      const_row_iter_t(mat_data, mat_data.size1()),
      basis_vectors.begin(),
      &vec_initial_point));
    elapsed = clock_type::now() - start;
  }
  else
  {
    clock_type::time_point start = clock_type::now(); 
    BOOST_CHECK(instance.batch_process(
      max_iterations,
      5,
      const_row_iter_t(mat_data, 0),
      const_row_iter_t(mat_data, mat_data.size1()),
      basis_vectors.begin()));
    elapsed = clock_type::now() - start;
  }
  
  std::cout << "processing " << nb_elements << " elements "
    << "in " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed) << std::endl;


  // testing the output sizes
  BOOST_REQUIRE_EQUAL(basis_vectors.size(), DATA_DIMENSION == 5 ? dimensions : 5);
  for(int i = 0; i < basis_vectors.size(); i++)
  {
    BOOST_CHECKPOINT("testing basis vector size for vector " << i);
    BOOST_REQUIRE_EQUAL(basis_vectors[i].size(), dimensions);
  }


  if(DATA_DIMENSION <= 5)
  {
    BOOST_MESSAGE("Generated basis vectors are:");

    for(int i = 0; i < dimensions; i++)
    {
      BOOST_MESSAGE("vector " << i << " :" << basis_vectors[i]);
    }
  }

  // testing orthogonality of all basis vectors
  for(int i = 0; i < basis_vectors.size() - 1; i++)
  {
    for(int j = i + 1; j < basis_vectors.size(); j++)
    {
      BOOST_CHECK_LE(ub::inner_prod(basis_vectors[i], basis_vectors[j]), 1E-6);
    }
  }


  // testing unitarity
  for(int i = 0; i < basis_vectors.size(); i++)
  {
    BOOST_CHECK_CLOSE(ub::inner_prod(basis_vectors[i], basis_vectors[i]), 1, 1E-6);
  }

  if(DATA_DIMENSION == 5)
  {
    // testing against the matlab output, the script being given by Sorent Hauberg, and the init between
    // each dimension iteration being given by the vector "initial_point" above. Each column represents 
    // an basis vector.
    static const double matlab_data[] = {
      -0.0365,   -0.0529,    0.0556,    0.9964,   -0.0058,
       0.9983,   -0.0457,    0.0088,    0.0337,    0.0073,
      -0.0076,   -0.0120,    0.9982,   -0.0565,    0.0165,
      -0.0084,   -0.0242,   -0.0166,    0.0051,    0.9995,
       0.0436,    0.9972,    0.0150,    0.0538,    0.0245,
    };


    for(int i = 0; i < dimensions; i++)
    {
      ub::vector<double> current_matlab_vector(dimensions);
      for(int j = 0; j < dimensions; j++)
      {
        current_matlab_vector(j) = matlab_data[i + j*dimensions];
      }
      BOOST_CHECKPOINT("iteration " << i);
      BOOST_CHECK_LE(ub::norm_2(basis_vectors[i] - current_matlab_vector), 2.6E-2); 
      // there is a slight difference between the two implementation when enabling the P^2 algorithm
      // for producing the quantiles.
    }

  }


}



BOOST_AUTO_TEST_SUITE_END();


