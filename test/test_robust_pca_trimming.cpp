// Copyright 2014, Max Planck Institute for Intelligent Systems.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

/*!@file
 * This file includes the tests for the trimmed version of the robust pca
 */

#include <boost/test/unit_test.hpp>
#include <test/test_main.hpp>

#include <include/robust_pca_trimming.hpp>
#include <include/private/boost_ublas_matrix_helper.hpp>

// data stored into a matrix
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

// generating data randomly
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>


// boost chrono
#include <boost/chrono/include.hpp>

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
  typedef boost::chrono::steady_clock clock_type;
  

  typedef robust_pca_with_trimming_impl< ub::vector<double> > robust_pca_t;
  robust_pca_t instance(0.);
  typedef row_iter<const matrix_t> const_row_iter_t;

  typedef ub::vector<double> data_t;

  BOOST_CHECK(instance.set_nb_processors(1));


  std::vector<data_t> temporary_data(nb_elements);
  std::vector<data_t> eigen_vectors(dimensions);
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
      temporary_data.begin(),
      eigen_vectors.begin(),
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
      temporary_data.begin(),
      eigen_vectors.begin()));
    elapsed = clock_type::now() - start;
  }
  
  
  std::cout << "processing " << nb_elements << " elements "
    << "in " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed) << "microseconds" << std::endl;
  


  // testing the output sizes
  BOOST_REQUIRE_EQUAL(eigen_vectors.size(), dimensions);
  for(int i = 0; i < dimensions; i++)
  {
    BOOST_CHECKPOINT("testing eigenvector size for vector " << i);
    BOOST_REQUIRE_EQUAL(eigen_vectors[i].size(), dimensions);
  }


  if(DATA_DIMENSION <= 5)
  {
    BOOST_MESSAGE("Generated eigen vectors are:");

    for(int i = 0; i < dimensions; i++)
    {
      BOOST_MESSAGE("vector " << i << " :" << eigen_vectors[i]);
    }
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

  if(DATA_DIMENSION == 5)
  {
    // testing against the matlab output, the script being given by Sorent Hauberg, and the init between
    // each dimension iteration being given by the vector "initial_point" above. Each column represents 
    // an eigenvector.
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
      BOOST_CHECK_LE(ub::norm_2(eigen_vectors[i] - current_matlab_vector), 2.6E-2); 
      // there is a slight difference between the two implementation when enabling the P^2 algorithm
      // for producing the quantiles.
    }

  }


}



BOOST_AUTO_TEST_CASE(smoke_and_orthogonality_tests_several_workers)
{

  using namespace robust_pca;
  using namespace robust_pca::ublas_adaptor;
  namespace ub = boost::numeric::ublas;
  typedef boost::chrono::steady_clock clock_type;
  

  typedef robust_pca_with_trimming_impl< ub::vector<double> > robust_pca_t;
  robust_pca_t instance(0.);
  typedef row_iter<const matrix_t> const_row_iter_t;

  typedef ub::vector<double> data_t;


  std::vector<data_t> temporary_data(nb_elements);
  std::vector<data_t> eigen_vectors(dimensions);
  const int max_iterations = 1000;

  BOOST_CHECK(instance.set_nb_processors(7)); // each chunk is floor(1000/7) = 142. Last chunk is 148. Just to test sthg different from 10.

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
      temporary_data.begin(),
      eigen_vectors.begin(),
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
      temporary_data.begin(),
      eigen_vectors.begin()));
    elapsed = clock_type::now() - start;
  }
  
  std::cout << "processing " << nb_elements << " elements "
    << "in " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed) << "microseconds" << std::endl;


  // testing the output sizes
  BOOST_REQUIRE_EQUAL(eigen_vectors.size(), dimensions);
  for(int i = 0; i < dimensions; i++)
  {
    BOOST_CHECKPOINT("testing eigenvector size for vector " << i);
    BOOST_REQUIRE_EQUAL(eigen_vectors[i].size(), dimensions);
  }


  if(DATA_DIMENSION <= 5)
  {
    BOOST_MESSAGE("Generated eigen vectors are:");

    for(int i = 0; i < dimensions; i++)
    {
      BOOST_MESSAGE("vector " << i << " :" << eigen_vectors[i]);
    }
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

  if(DATA_DIMENSION == 5)
  {
    // testing against the matlab output, the script being given by Sorent Hauberg, and the init between
    // each dimension iteration being given by the vector "initial_point" above. Each column represents 
    // an eigenvector.
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
      BOOST_CHECK_LE(ub::norm_2(eigen_vectors[i] - current_matlab_vector), 2.6E-2); 
      // there is a slight difference between the two implementation when enabling the P^2 algorithm
      // for producing the quantiles.
    }

  }


}



BOOST_AUTO_TEST_SUITE_END();



template <class T>
std::vector<T>* readLines(std::istream& str)
{
  std::vector<T>* result = new std::vector<T>();
  std::string line;
  std::getline(str,line);

  std::stringstream lineStream(line);
  std::string cell;

  
  while(std::getline(lineStream, cell, '\t'))
  {
    std::stringstream current_stream(cell);
    T current_element;
    current_stream >> current_element;
    assert(!current_stream.fail());
    result->push_back(current_element);
  }
  return result;
}



BOOST_AUTO_TEST_CASE(convergence_rate_tests_several_workers)
{

  using namespace robust_pca;
  using namespace robust_pca::ublas_adaptor;
  namespace ub = boost::numeric::ublas;
  typedef boost::chrono::steady_clock clock_type;
  typedef boost::numeric::ublas::matrix<double> matrix_t;


  typedef robust_pca_with_trimming_impl< ub::vector<double> > robust_pca_t;
  robust_pca_t instance(.1);
  typedef row_iter<const matrix_t> const_row_iter_t;

  typedef ub::vector<double> data_t;

  std::string const filename_to_read = "D:/Code/mat_test.csv";
  std::ifstream ff(filename_to_read.c_str());

  std::vector< std::vector<double>* > read_vectors;

  std::cout << "Reading data" << std::endl;
  matrix_t mat_data;
  while(!ff.eof())
  {
    std::vector<double>* v = readLines<double>(ff);
    if(v->size())
    {
      read_vectors.push_back(v);
      if((read_vectors.size() % 1000) == 0)
      {
        std::cout << ".";
        std::cout.flush();
      }
    }
  }
  std::cout << std::endl << "copying" << std::endl;

  mat_data.resize(read_vectors.size(), read_vectors[0]->size());
  for(size_t i = 0; i < read_vectors.size(); i++)
  {
    std::copy(read_vectors[i]->begin(), read_vectors[i]->end(), ub::row(mat_data, i).begin());
    delete read_vectors[i];
  }
  read_vectors.clear();

  const size_t nb_elements = mat_data.size1();
  const size_t dimensions = mat_data.size2();
  std::cout << "Data ok : dimensions = " << dimensions << " #elements = " << nb_elements << std::endl;


  const size_t max_dimensions = 5;

  std::vector<data_t> temporary_data(nb_elements);
  std::vector<data_t> eigen_vectors(max_dimensions);
  const int max_iterations = 1000;

  BOOST_CHECK(instance.set_nb_processors(7)); // each chunk is floor(1000/7) = 142. Last chunk is 148. Just to test sthg different from 10.

  clock_type::duration elapsed;

  clock_type::time_point start = clock_type::now(); 
  BOOST_CHECK(instance.batch_process(
    max_iterations,
    max_dimensions,
    const_row_iter_t(mat_data, 0),
    const_row_iter_t(mat_data, mat_data.size1()),
    temporary_data.begin(),
    eigen_vectors.begin()));
  elapsed = clock_type::now() - start;
  
  std::cout << "processing " << nb_elements << " elements "
    << "in " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed) << "microseconds" << std::endl;


  // testing the output sizes
  BOOST_REQUIRE_EQUAL(eigen_vectors.size(), max_dimensions);
  for(int i = 0; i < max_dimensions; i++)
  {
    BOOST_CHECKPOINT("testing eigenvector size for vector " << i);
    BOOST_REQUIRE_EQUAL(eigen_vectors[i].size(), dimensions);
  }

  // testing orthogonality of all eigenvectors
  for(int i = 0; i < max_dimensions - 1; i++)
  {
    for(int j = i + 1; j < max_dimensions; j++)
    {
      BOOST_CHECK_LE(ub::inner_prod(eigen_vectors[i], eigen_vectors[j]), 1E-6);
    }
  }


  // testing unitarity
  for(int i = 0; i < max_dimensions; i++)
  {
    BOOST_CHECK_CLOSE(ub::inner_prod(eigen_vectors[i], eigen_vectors[i]), 1, 1E-6);
  }


}


