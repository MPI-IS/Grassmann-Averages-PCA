// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

/*!@file
 * This file includes the tests for the grassmann pca
 */


#include <boost/test/unit_test.hpp>

#include <test/test_main.hpp>

#include <include/grassmann_pca.hpp>
#include <include/private/boost_ublas_row_iterator.hpp>

// data stored into a matrix
#include <boost/numeric/ublas/io.hpp>

// boost chrono
#include <boost/chrono/include.hpp>
#include <fstream>


BOOST_FIXTURE_TEST_SUITE(grassmann_pca_test_suite, fixture_simple_matrix_creation)

BOOST_AUTO_TEST_CASE(returns_false_for_inapropriate_inputs)
{
  using namespace grassmann_averages_pca;
  namespace ub = boost::numeric::ublas;

  typedef grassmann_pca< ub::vector<double> > grassmann_pca_t;

  grassmann_pca_t instance;


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

  typedef grassmann_pca< ub::vector<double>, observer_t > grassmann_pca_t;

  grassmann_pca_t instance;

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

  typedef grassmann_pca< ub::vector<double>, observer_t > grassmann_pca_t;

  grassmann_pca_t instance;

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


  typedef grassmann_pca< ub::vector<double> > grassmann_pca_t;  
  grassmann_pca_t instance;
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



  if(DATA_DIMENSION == 5)
  {
    // testing against the matlab output, the script being given by Sorent Hauberg, and the init between
    // each dimension iteration being given by the vector "initial_point" above. Each column represents 
    // an basis vector.
    static const double matlab_data[] = {
      -0.0318,    0.0564,   -0.0290,   -0.0136,    0.9974,
       0.0242,    0.0061,    0.9993,   -0.0013,    0.0295,
      -0.0118,    0.9982,   -0.0041,    0.0153,   -0.0567,
      -0.0307,   -0.0149,    0.0017,    0.9993,    0.0136,
       0.9987,    0.0130,   -0.0252,    0.0305,    0.0308,
    };


    for(int i = 0; i < dimensions; i++)
    {
      ub::vector<double> current_matlab_vector(dimensions);
      for(int j = 0; j < dimensions; j++)
      {
        current_matlab_vector(j) = matlab_data[i + j*dimensions];
      }
      BOOST_TEST_CHECKPOINT("iteration " << i);
      BOOST_CHECK_LE(ub::norm_2(basis_vectors[i] - current_matlab_vector), 1E-3);
      //std::cout << "computed = " << basis_vectors[i] << std::endl;
      //std::cout << "matlab = " << current_matlab_vector << std::endl;
    }
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


  typedef grassmann_pca< ub::vector<double> > grassmann_pca_t;  
  grassmann_pca_t instance;
  typedef row_iter<const matrix_t> const_row_iter_t;
  
  typedef ub::vector<double> data_t;


  std::vector<data_t> basis_vectors(DATA_DIMENSION == 5 ? dimensions : 5);
  const int max_iterations = 1000;


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

  if(DATA_DIMENSION == 5)
  {
    // testing against the matlab output, the script being given by Sorent Hauberg, and the init between
    // each dimension iteration being given by the vector "initial_point" above. Each column represents 
    // an basis vector.
    static const double matlab_data[] = {
      -0.0318,    0.0564,   -0.0290,   -0.0136,    0.9974,
       0.0242,    0.0061,    0.9993,   -0.0013,    0.0295,
      -0.0118,    0.9982,   -0.0041,    0.0153,   -0.0567,
      -0.0307,   -0.0149,    0.0017,    0.9993,    0.0136,
       0.9987,    0.0130,   -0.0252,    0.0305,    0.0308,
    };


    for(int i = 0; i < dimensions; i++)
    {
      ub::vector<double> current_matlab_vector(dimensions);
      for(int j = 0; j < dimensions; j++)
      {
        current_matlab_vector(j) = matlab_data[i + j*dimensions];
      }
      BOOST_TEST_CHECKPOINT("iteration " << i);
      BOOST_CHECK_LE(ub::norm_2(basis_vectors[i] - current_matlab_vector), 1E-3);
    }
  }  



}



#if 0
BOOST_AUTO_TEST_CASE(checking_against_matlab)
{
  using namespace grassmann_averages_pca;
  using namespace grassmann_averages_pca::ublas_adaptor;
  namespace ub = boost::numeric::ublas;


  typedef grassmann_pca< ub::vector<double> > grassmann_pca_t;  
  grassmann_pca_t instance;
  typedef row_iter<const matrix_t> const_row_iter_t;
  
  typedef boost::numeric::ublas::vector<double> data_t;

  std::vector<data_t> temporary_data(nb_elements);
  std::vector<data_t> basis_vectors(dimensions);
  const int max_iterations = 1000;


  BOOST_CHECK(instance.batch_process(
    max_iterations,
    dimensions,
    const_row_iter_t(mat_data, 0),
    const_row_iter_t(mat_data, mat_data.size1()),
    temporary_data.begin(),
    basis_vectors.begin()));


  BOOST_TEST_MESSAGE("Generated basis vectors are:");

  for(int i = 0; i < dimensions; i++)
  {
    BOOST_TEST_MESSAGE("vector " << i << " :" << basis_vectors[i]);
  }
}
#endif

BOOST_AUTO_TEST_SUITE_END();
