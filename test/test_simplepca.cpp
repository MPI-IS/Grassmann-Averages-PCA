// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

/*!@file
 * This file includes the tests for the grassmann pca
 */


#include <boost/test/unit_test.hpp>

#include <test/test_main.hpp>

#include <include/private/em_pca.hpp>
#include <include/private/boost_ublas_row_iterator.hpp>

// data stored into a matrix
#include <boost/numeric/ublas/io.hpp>

// boost chrono
#include <boost/chrono/include.hpp>



#include <fstream>






BOOST_FIXTURE_TEST_SUITE(simple_pca_test_suite, fixture_simple_matrix_creation)

BOOST_AUTO_TEST_CASE(returns_false_for_inapropriate_inputs)
{
  using namespace grassmann_averages_pca;
  namespace ub = boost::numeric::ublas;

  typedef em_pca< ub::vector<double> > em_pca_t;

  em_pca_t instance;


  typedef details::ublas_helpers::row_iter<const matrix_t> const_row_iter_t;
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


BOOST_AUTO_TEST_CASE(check_centering_not_called)
{
  // in this implementation, the centering should not be called
  // the observer throws an exception that is intercepted by the test here.
  using namespace grassmann_averages_pca;
  namespace ub = boost::numeric::ublas;
  typedef ub::vector<double> data_t;

  typedef test_mean_observer<data_t> observer_t;

  typedef em_pca< ub::vector<double>, observer_t > em_pca_t;

  em_pca_t instance;

  observer_t observer;

  BOOST_CHECK(instance.set_observer(&observer));
  BOOST_CHECK(instance.set_centering(false));


  typedef details::ublas_helpers::row_iter<const matrix_t> const_row_iter_t;

  std::vector<data_t> basis_vectors(dimensions);
  const int max_iterations = 10;

  BOOST_CHECK_NO_THROW(instance.batch_process(
                       max_iterations,
                       dimensions,
                       const_row_iter_t(mat_data, 0),
                       const_row_iter_t(mat_data, 10),
                       basis_vectors.begin()));

}



BOOST_AUTO_TEST_CASE(check_centering_of_data)
{
  // checks that the multithreaded centering is performing well
  // the observer throws an exception that is intercepted by the test here.
  using namespace grassmann_averages_pca;
  namespace ub = boost::numeric::ublas;
  typedef ub::vector<double> data_t;

  typedef test_mean_observer<data_t> observer_t;

  typedef em_pca< ub::vector<double>, observer_t > em_pca_t;

  em_pca_t instance;

  observer_t observer;

  BOOST_CHECK(instance.set_observer(&observer));
  BOOST_CHECK(instance.set_centering(true));


  typedef details::ublas_helpers::row_iter<const matrix_t> const_row_iter_t;

  std::vector<data_t> basis_vectors(dimensions);
  const int max_iterations = 1000;

  BOOST_CHECK_THROW(instance.batch_process(
                      max_iterations,
                      dimensions,
                      const_row_iter_t(mat_data, 0),
                      const_row_iter_t(mat_data, 10),
                      basis_vectors.begin()), 
                    observer_t::s_signal_exception);

  for(int j = 0; j < dimensions; j++)
  {
    double acc = 0;
    for(int i = 0; i < 10; i++)
    {
      acc += mat_data(i, j);
    }

    BOOST_CHECK_CLOSE(acc / 10, observer.mean(j), 1E-3);
  }


}



BOOST_AUTO_TEST_CASE(smoke_and_orthogonality_tests)
{
  using namespace grassmann_averages_pca;
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;
  typedef boost::chrono::steady_clock clock_type;


  typedef em_pca< ub::vector<double> > em_pca_t;  
  em_pca_t instance;
  typedef row_iter<const matrix_t> const_row_iter_t;
  
  typedef ub::vector<double> data_t;


  std::vector<data_t> basis_vectors(dimensions);
  const int max_iterations = 1000;


  clock_type::duration elapsed;
  if(DATA_DIMENSION == 5)
  {
    const double initial_point[] = {0.2097, 0.3959, 0.5626, 0.2334, 0.6545};
    BOOST_REQUIRE_EQUAL(dimensions, sizeof(initial_point)/sizeof(initial_point[0])); // just in case

    ub::vector<double> vec_initial_point(dimensions);
    for(int i = 0; i < dimensions; i++)
    {
      vec_initial_point(i) = initial_point[i];
    }

    std::vector< ub::vector<double> > v_init_points(dimensions, vec_initial_point);
    clock_type::time_point start = clock_type::now();
    BOOST_CHECK(instance.batch_process(
      max_iterations,
      dimensions,
      const_row_iter_t(mat_data, 0),
      const_row_iter_t(mat_data, mat_data.size1()),
      basis_vectors.begin(),
      &v_init_points));
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
    BOOST_TEST_CHECKPOINT("testing basis vector size for vector " << i);
    BOOST_REQUIRE_EQUAL(basis_vectors[i].size(), dimensions);
  }


  if(DATA_DIMENSION <= 5)
  {
    BOOST_TEST_MESSAGE("Generated basis vectors are:");

    for(int i = 0; i < dimensions; i++)
    {
      BOOST_TEST_MESSAGE("vector " << i << " :" << basis_vectors[i]);
    }
  }

  // testing orthogonality of all basis vectors
  for(int i = 0; i < dimensions-1; i++)
  {
    for(int j = i + 1; j < dimensions; j++)
    {
      BOOST_CHECK_LE(ub::inner_prod(basis_vectors[i], basis_vectors[j]), 1E-6);
    }
  }

  // testing unitarity of all basis vectors
  for(int i = 0; i < dimensions; i++)
  {
    BOOST_CHECK_CLOSE(ub::inner_prod(basis_vectors[i], basis_vectors[i]), 1, 1E-6);
  }

}


BOOST_AUTO_TEST_CASE(smoke_and_orthogonality_tests_several_workers)
{
  // this test uses several worker threads, but should provide exactly the same values at the previous test. 
  // its body is almost the same.
  using namespace grassmann_averages_pca;
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;
  typedef boost::chrono::steady_clock clock_type;


  typedef em_pca< ub::vector<double> > em_pca_t;  
  em_pca_t instance;
  typedef row_iter<const matrix_t> const_row_iter_t;
  
  typedef ub::vector<double> data_t;


  std::vector<data_t> basis_vectors(DATA_DIMENSION == 5 ? dimensions : 5);
  const int max_iterations = 1000;

  BOOST_CHECK(instance.set_centering(true));

  BOOST_CHECK(instance.set_nb_processors(7)); // each chunk is floor(1000/7) = 142. Last chunk is 148. Just to test sthg different from 10.

  
  clock_type::duration elapsed;
  if(DATA_DIMENSION == 5)
  {
    const double initial_point[] = {0.2097, 0.3959, 0.5626, 0.2334, 0.6545};
    BOOST_REQUIRE_EQUAL(dimensions, sizeof(initial_point)/sizeof(initial_point[0])); // just in case

    ub::vector<double> vec_initial_point(dimensions);
    for(int i = 0; i < dimensions; i++)
    {
      vec_initial_point(i) = initial_point[i];
    }

    std::vector< ub::vector<double> > v_init_points(dimensions, vec_initial_point);

    // setting the number of workers

    // main call
    clock_type::time_point start = clock_type::now();
    BOOST_CHECK(instance.batch_process(
      max_iterations,
      dimensions,
      const_row_iter_t(mat_data, 0),
      const_row_iter_t(mat_data, mat_data.size1()),
      basis_vectors.begin(),
      &v_init_points));
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
    BOOST_TEST_CHECKPOINT("testing basis vector size for vector " << i);
    BOOST_REQUIRE_EQUAL(basis_vectors[i].size(), dimensions);
  }

  if(DATA_DIMENSION <= 5)
  {
    BOOST_TEST_MESSAGE("Generated basis vectors are:");

    for(int i = 0; i < dimensions; i++)
    {
      BOOST_TEST_MESSAGE("vector " << i << " :" << basis_vectors[i]);
    }
  }

  // testing orthogonality of all basis vectors
  for(int i = 0; i < basis_vectors.size()-1; i++)
  {
    for(int j = i + 1; j < basis_vectors.size(); j++)
    {
      BOOST_CHECK_LE(ub::inner_prod(basis_vectors[i], basis_vectors[j]), 1E-6);
    }
  }

  for(int i = 0; i < basis_vectors.size(); i++)
  {
    BOOST_CHECK_CLOSE(ub::inner_prod(basis_vectors[i], basis_vectors[i]), 1, 1E-6);
  }


}



BOOST_AUTO_TEST_SUITE_END();

