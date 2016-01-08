// Copyright 2014, Max Planck Society.
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

#include <include/grassmann_pca_with_trimming.hpp>
#include <include/private/boost_ublas_row_iterator.hpp>


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
  
  BOOST_TEST_MESSAGE("low heap top = " << lowh.top());
  for(low_heap_t::ordered_iterator it(lowh.ordered_begin()), ite(lowh.ordered_end());
      it != ite;
      ++it)
  {
    BOOST_TEST_MESSAGE(" " << *it);
  }  

  BOOST_TEST_MESSAGE("high heap top = " << highh.top());
  for(high_heap_t::ordered_iterator it(highh.ordered_begin()), ite(highh.ordered_end());
      it != ite;
      ++it)
  {
    BOOST_TEST_MESSAGE(" " << *it);
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

  std::vector<low_heap_t> v(10, lowh);
  for(int i = 0; i < 10; i++)
  {
    BOOST_CHECK_EQUAL(v[i].size(), K);
  }
}
