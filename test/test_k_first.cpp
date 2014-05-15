// Copyright 2014, Max Planck Institute for Intelligent Systems.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

/*!@file
 * This file contains some tests for the heaps implementation
 */

#include <boost/test/unit_test.hpp>
#include <test/test_main.hpp>

// data stored into a matrix
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>


// boost heaps
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/heap/pairing_heap.hpp>

// boost chrono
#include <boost/chrono/include.hpp>

#include <include/robust_pca_trimming.hpp>
#include <include/private/boost_ublas_matrix_helper.hpp>


BOOST_AUTO_TEST_CASE(test_heap)
{
  // a simple test for looking at the heaps functionalities
  typedef boost::heap::fibonacci_heap<double, boost::heap::compare< std::less<double> > > low_heap_t;
  typedef boost::heap::fibonacci_heap<double, boost::heap::compare< std::greater<double> > > high_heap_t;
  
  low_heap_t lowh;
  high_heap_t highh;
  
  const int K = 7;
  
  boost::random::uniform_real_distribution<double> dist;
  for(int i = 0; i < 100; i++)
  {
    double current = dist(rng);
    if(lowh.size() < K)
    {
      lowh.push(current);
    }
    else if(lowh.value_comp()(current, lowh.top()))
    {
      lowh.push(current);
      lowh.pop();
    }
    
    if(highh.size() < K)
    {
      highh.push(current);
    }
    else if(highh.value_comp()(current, highh.top()))
    {
      highh.push(current);
      highh.pop();
    }    
  }
  
  BOOST_CHECK_EQUAL(lowh.size(), K);
  BOOST_CHECK_EQUAL(highh.size(), K);
  
  BOOST_CHECK_LE(lowh.top(), highh.top());
  
  BOOST_MESSAGE("low heap top = " << lowh.top());
  for(low_heap_t::ordered_iterator it(lowh.ordered_begin()), ite(lowh.ordered_end());
      it != ite;
      ++it)
  {
    BOOST_MESSAGE(" " << *it);
  }  

  BOOST_MESSAGE("high heap top = " << highh.top());
  for(high_heap_t::ordered_iterator it(highh.ordered_begin()), ite(highh.ordered_end());
      it != ite;
      ++it)
  {
    BOOST_MESSAGE(" " << *it);
  }  
  
}



BOOST_AUTO_TEST_CASE(test_heap_vector)
{
  typedef boost::heap::pairing_heap<double, boost::heap::compare< std::less<double> > > low_heap_t;
  
  low_heap_t lowh;
  
  const int K = 7;
  
  boost::random::uniform_real_distribution<double> dist;
  for(int i = 0; i < 100; i++)
  {
    double current = dist(rng);
    if(lowh.size() < K)
    {
      lowh.push(current);
    }
    else if(lowh.value_comp()(current, lowh.top()))
    {
      lowh.push(current);
      lowh.pop();
    }
  }
  
  BOOST_CHECK_EQUAL(lowh.size(), K);
#if 1
  std::vector<low_heap_t> v(10, lowh);
  for(int i = 0; i < 10; i++)
  {
    BOOST_CHECK_EQUAL(v[i].size(), K);
  }
#endif
}

BOOST_AUTO_TEST_CASE(test_double_heap_performances)
{
  namespace ub = boost::numeric::ublas;
  using namespace robust_pca::ublas_adaptor;
  //using namespace robust_pca::ublas_matlab_helper;

  typedef ub::matrix<double> matrix_t;
  typedef boost::chrono::steady_clock clock_type;

  const int nb_elements = 100000;
  const int dimensions = 500;
  const int nb_elements_to_keep = 50;
  const int nb_chunks = 100;

  matrix_t mat_data;

  boost::random::uniform_real_distribution<double> dist;

  typedef row_iter<matrix_t> row_iter_t;

  row_iter_t const it1(mat_data, 0);
  row_iter_t const ite(mat_data, mat_data.size1());



  mat_data.resize(nb_elements, dimensions);

  for(int i = 0; i < nb_elements; i++)
  {
    for(int j = 0; j < dimensions; j++)
    {
      mat_data(i, j) = dist(rng);
    }
  }


  typedef ub::vector<double> data_t;
  typedef data_t::value_type scalar_t;

  robust_pca::details::s_double_heap_vector<data_t> double_heap_merged;
  

  clock_type::duration elapsed;
  clock_type::duration elapsed_push;
  clock_type::duration elapsed_push_ignore;
  clock_type::duration elapsed_merge;
  clock_type::time_point start = clock_type::now();
  
  double_heap_merged.set_dimension(dimensions);
  
  for(int i = 0; i < nb_chunks; i++)
  {
    robust_pca::details::s_double_heap_vector<data_t> double_heap;
    double_heap.set_dimension(dimensions);
    row_iter_t it_data(it1 + i*nb_elements/nb_chunks);
    for(size_t s = 0; s < nb_elements/nb_chunks; ++it_data, s++)
    {
      row_iter_t::reference current_data = *it_data;
    
      if(s < nb_elements_to_keep)
      { 
        clock_type::time_point start_push = clock_type::now();
        double_heap.push(current_data, s % 2 == 0);
        elapsed_push += clock_type::now() - start_push;
      }
      else
      {
        clock_type::time_point start_push = clock_type::now();
        double_heap.push_or_ignore(current_data, s % 2 == 0);
        elapsed_push_ignore += clock_type::now() - start_push;
      }
    }

    {
      clock_type::time_point start_merge = clock_type::now();
      double_heap_merged.merge(double_heap);
      elapsed_merge += clock_type::now() - start_merge;
    }
  }

  elapsed = clock_type::now() - start;


  std::cout << "processing " << nb_elements << " elements "
    << "in " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed) << std::endl;
  std::cout << "\tpush " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed_push) << std::endl;
  std::cout << "\tpush-ignore " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed_push_ignore) << std::endl;
  std::cout << "\tmerge " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed_merge) << std::endl;


  std::vector<double> lower_bounds(dimensions, std::numeric_limits<double>::max());
  std::vector<double> upper_bounds(dimensions, -std::numeric_limits<double>::max());
  clock_type::time_point start_std = clock_type::now();
  for(int i = 0; i < dimensions; i++)
  {
    std::vector<double> copy(ub::matrix_column<matrix_t>(mat_data, i).begin(), ub::matrix_column<matrix_t>(mat_data, i).end());
    std::partial_sort(copy.begin(), copy.begin()+nb_elements_to_keep, copy.end());
    lower_bounds[i] = copy[nb_elements_to_keep];
    copy.erase(copy.begin(), copy.begin()+nb_elements_to_keep);

    std::partial_sort(copy.begin(), copy.begin()+nb_elements_to_keep, copy.end(), std::greater<double>());
    upper_bounds[i] = copy[nb_elements_to_keep];
  }

  clock_type::duration elapsed_std = clock_type::now() - start_std;
  std::cout << "\tstd " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed_std) << std::endl;


  clock_type::time_point start_std_t = clock_type::now();
  matrix_t transposed = ub::trans(mat_data);
  for(int i = 0; i < dimensions; i++)
  {
    std::vector<double> copy(ub::matrix_row<matrix_t>(mat_data, i).begin(), ub::matrix_row<matrix_t>(mat_data, i).end());
    std::partial_sort(copy.begin(), copy.begin()+nb_elements_to_keep, copy.end());
    lower_bounds[i] = copy[nb_elements_to_keep];

    copy.erase(copy.begin(), copy.begin()+nb_elements_to_keep);

    std::partial_sort(copy.begin(), copy.begin()+nb_elements_to_keep, copy.end(), std::greater<double>());
    upper_bounds[i] = copy[nb_elements_to_keep];
  }

  clock_type::duration elapsed_std_t = clock_type::now() - start_std_t;
  std::cout << "\tstd transposed " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed_std_t) << std::endl;


}