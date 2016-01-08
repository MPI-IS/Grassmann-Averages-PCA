// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

/*!@file
 * This file includes an acceptance tests for the trimmed version of the grassmann pca
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
// data stored into a matrix
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


// boost chrono
#include <boost/chrono/include.hpp>

#include <fstream>
#include <sstream>

#include <include/grassmann_pca_with_trimming.hpp>
#include <include/private/boost_ublas_row_iterator.hpp>



std::string filename_data;
std::string filename_basis_vectors;
std::string filename_expected_result;


struct ConfigurationFiles
{

  ConfigurationFiles()
  {
    int argc  = boost::unit_test::framework::master_test_suite().argc;
    char**argv= boost::unit_test::framework::master_test_suite().argv;

    // from the command line, we expect the file that should be read.
    for(int i = 0; i < argc - 1; i++)
    {
      if(std::string(argv[i]) == "--data")
      {
        if(!filename_data.empty())
        {
          std::cerr << "Test initialisation error: data given several times" << std::endl;
        }
        filename_data = argv[i+1];
      }
      if(std::string(argv[i]) == "--basis_vectors")
      {
        if(!filename_basis_vectors.empty())
        {
          std::cerr << "Test initialisation error: basis_vectors given several times" << std::endl;
        }
        filename_basis_vectors = argv[i+1];
      }
      if(std::string(argv[i]) == "--expected_result")
      {
        if(!filename_expected_result.empty())
        {
          std::cerr << "Test initialisation error: expected_result given several times" << std::endl;
        }
        filename_expected_result = argv[i+1];
      }
    }
  }
};

BOOST_GLOBAL_FIXTURE( ConfigurationFiles );







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

template <class matrix_t>
bool read_matrix(const std::string &filename, matrix_t &mat_data)
{
  std::vector< std::vector<double>* > read_vectors;
  std::ifstream ff(filename.c_str());
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

  mat_data.resize(read_vectors.size(), read_vectors[0]->size());
  for(size_t i = 0; i < read_vectors.size(); i++)
  {
    std::copy(read_vectors[i]->begin(), read_vectors[i]->end(), boost::numeric::ublas::row(mat_data, i).begin());
    delete read_vectors[i];
  }
  read_vectors.clear();

  return true;
}

BOOST_AUTO_TEST_CASE(convergence_rate_tests_several_workers)
{

  using namespace grassmann_averages_pca;
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;
  typedef boost::chrono::steady_clock clock_type;
  typedef boost::numeric::ublas::matrix<double> matrix_t;


  BOOST_REQUIRE(!filename_data.empty());

  typedef ub::vector<double> data_t;
  typedef grassmann_pca_with_trimming< data_t > grassmann_pca_t;
  grassmann_pca_t instance(.1);
  typedef row_iter<const matrix_t> const_row_iter_t;

  

  // reading the data, one vector per line
  std::cout << "Reading data" << std::endl;
  matrix_t mat_data;
  BOOST_REQUIRE(read_matrix(filename_data, mat_data));



  std::cout << std::endl << "Reading init basis vectors" << std::endl;
  // reading the init vectors, one vector per column
  std::vector<data_t> v_basisvectors_init;
  matrix_t mat_vectors;
  BOOST_REQUIRE(read_matrix(filename_basis_vectors, mat_vectors));
  for(int i = 0; i < mat_vectors.size2(); i++)
  {
    v_basisvectors_init.push_back(ub::column(mat_vectors, i));
  }

  std::cout << std::endl << "Reading expected result" << std::endl;
  std::vector<data_t> v_known_result;
  matrix_t mat_result;
  BOOST_REQUIRE(read_matrix(filename_expected_result, mat_result));
  for(int i = 0; i < mat_vectors.size2(); i++)
  {
    v_known_result.push_back(ub::column(mat_result, i));
  }


  const size_t nb_elements = mat_data.size1();
  const size_t dimensions = mat_data.size2();
  std::cout << std::endl << "Data: dimensions = " << dimensions << " #elements = " << nb_elements << std::endl;


  const size_t max_dimensions = mat_vectors.size2();

  std::vector<data_t> basis_vectors(max_dimensions);
  const int max_iterations = 1000;


  BOOST_CHECK(instance.set_nb_processors(7));
  BOOST_CHECK(instance.set_nb_steps_pca(0));

  clock_type::duration elapsed;

  clock_type::time_point start = clock_type::now(); 
  BOOST_CHECK(instance.batch_process(
    max_iterations,
    max_dimensions,
    const_row_iter_t(mat_data, 0),
    const_row_iter_t(mat_data, mat_data.size1()),
    basis_vectors.begin(),
    &v_basisvectors_init));
  elapsed = clock_type::now() - start;
  
  std::cout << "processing " << nb_elements << " elements "
    << "in " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed) << std::endl;


  // testing the output sizes
  BOOST_REQUIRE_EQUAL(basis_vectors.size(), max_dimensions);
  for(int i = 0; i < max_dimensions; i++)
  {
    BOOST_TEST_CHECKPOINT("testing basis vector size for vector " << i);
    BOOST_REQUIRE_EQUAL(basis_vectors[i].size(), dimensions);
  }

  // testing orthogonality of all basis vectors
  for(int i = 0; i < max_dimensions - 1; i++)
  {
    for(int j = i + 1; j < max_dimensions; j++)
    {
      BOOST_CHECK_LE(ub::inner_prod(basis_vectors[i], basis_vectors[j]), 1E-6);
    }
  }


  // testing unitarity
  for(int i = 0; i < max_dimensions; i++)
  {
    BOOST_CHECK_CLOSE(ub::inner_prod(basis_vectors[i], basis_vectors[i]), 1, 1E-6);
  }

  // testing for coherence with matlab
  for(int i = 0; i < max_dimensions; i++)
  {
    details::norm_infinity norm;
    BOOST_CHECK_LE(norm(basis_vectors[i] - v_known_result[i]), 1E-2);
  }

}
