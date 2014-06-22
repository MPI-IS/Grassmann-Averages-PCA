// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)


/*!@file
 * This file contains tests for the row iterators. These iterators take some matrix like structure and provide 
 * iterators on the rows. These matrix-like structure can be matrices or some external array (eg. arrays that
 * are provided by matlab).
 */

#include <boost/test/unit_test.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <include/private/boost_ublas_row_iterator.hpp>
#include <include/private/boost_ublas_external_storage.hpp>

#include <boost/scoped_array.hpp>


BOOST_TEST_DONT_PRINT_LOG_VALUE(grassmann_averages_pca::details::ublas_helpers::row_iter<boost::numeric::ublas::matrix<double> >);
BOOST_TEST_DONT_PRINT_LOG_VALUE(grassmann_averages_pca::details::ublas_helpers::row_iter<const boost::numeric::ublas::matrix<double> >);

BOOST_AUTO_TEST_CASE(test_row_proxy_read)
{
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;

  typedef ub::matrix<double> matrix_t;
  typedef row_iter<const matrix_t> const_row_iter_t;


  matrix_t mat(3, 4);
  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 4; j++)
    {
      mat(i, j) = i * 10 + j;
    }
  }

  {
    const_row_iter_t it(mat, 0), ite(mat, 3);
    for(int i = 0; i < 3; i++, ++it)
    {
      BOOST_CHECK_NE(it, ite);
      const_row_iter_t::value_type v(*it);

      int j = 0;
      for(const_row_iter_t::value_type::iterator it1(v.begin()), it1e(v.end()); it1 != it1e; ++it1, j++)
      {
        BOOST_CHECK_EQUAL(*it1, i * 10 + j);
      }
    }
  }


}

BOOST_AUTO_TEST_CASE(test_row_proxy_advance)
{
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;

  typedef ub::matrix<double> matrix_t;
  typedef row_iter<matrix_t> row_iter_t;


  matrix_t mat(3, 4);

  row_iter_t it(mat, 0), ite(mat, 3);
  std::advance(it, 3);
  BOOST_CHECK(it == ite);
}




BOOST_AUTO_TEST_CASE(test_row_proxy_write)
{
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;

  typedef ub::matrix<double> matrix_t;
  typedef row_iter<matrix_t> row_iter_t;


  matrix_t mat(3, 4);

  {
    row_iter_t it(mat, 0), ite(mat, 3);
    for(int i = 0; i < 3; i++, ++it)
    {
      BOOST_CHECK_NE(it, ite);
      row_iter_t::value_type v(*it);

      int j = 0;
      for(row_iter_t::value_type::iterator it1(v.begin()), it1e(v.end()); it1 != it1e; ++it1, j++)
      {
        *it1 = i * 10 + j;
      }
    }
  }

  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 4; j++)
    {
      BOOST_CHECK_EQUAL(mat(i, j), i * 10 + j);
    }
  }
}


BOOST_AUTO_TEST_CASE(test_row_proxy_read_specific_storage)
{
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;

  typedef ub::matrix<double, ub::row_major, ub::bounded_array<double, 1000> > matrix_t;
  typedef row_iter<const matrix_t> const_row_iter_t;


  matrix_t mat(3, 4);
  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 4; j++)
    {
      mat(i, j) = i * 10 + j;
    }
  }

  {
    const_row_iter_t it(mat, 0), ite(mat, 3);
    for(int i = 0; i < 3; i++, ++it)
    {
      BOOST_CHECK(it != ite); // not using BOOST_CHECK_NE because the stream output should be disabled.
      const_row_iter_t::value_type v(*it);

      int j = 0;
      for(const_row_iter_t::value_type::iterator it1(v.begin()), it1e(v.end()); it1 != it1e; ++it1, j++)
      {
        BOOST_CHECK_EQUAL(*it1, i * 10 + j);
      }
    }
  }


}

// testing the row write access capability on for the row proxy, element wise
BOOST_AUTO_TEST_CASE(test_row_proxy_write_specific_storage)
{
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;

  typedef ub::matrix<double, ub::row_major, ub::bounded_array<double, 1000> > matrix_t;
  typedef row_iter<matrix_t> row_iter_t;

  matrix_t mat(3, 4);

  {
    row_iter_t it(mat, 0), ite(mat, 3);
    for(int i = 0; i < 3; i++, ++it)
    {
      BOOST_CHECK(it != ite);
      row_iter_t::value_type v(*it);

      int j = 0;
      for(row_iter_t::value_type::iterator it1(v.begin()), it1e(v.end()); it1 != it1e; ++it1, j++)
      {
        *it1 = i * 10 + j;
      }
    }
  }

  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 4; j++)
    {
      BOOST_CHECK_EQUAL(mat(i, j), i * 10 + j);
    }
  }

}


// testing the row write access capability for the row proxy, vector access
BOOST_AUTO_TEST_CASE(test_row_proxy_write_specific_storage_row_access)
{
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;

  typedef ub::matrix<double, ub::row_major, ub::bounded_array<double, 1000> > matrix_t;
  typedef row_iter<matrix_t> row_iter_t;

  matrix_t mat(3, 4);

  {
    row_iter_t it(mat, 0), ite(mat, 3);
    for(int i = 0; i < 3; i++, ++it)
    {
      BOOST_CHECK(it != ite);
      *it = ub::zero_vector<double>(4);
    }
  }

  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 4; j++)
    {
      BOOST_CHECK_EQUAL(mat(i, j), 0);
    }
  }

}



// testing the row write access capability for the row proxy, vector access, external storage
BOOST_AUTO_TEST_CASE(test_row_proxy_write_specific_external_storage)
{
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;

  typedef ub::matrix<double, ub::row_major, external_storage_adaptor<double> > matrix_t;
  typedef row_iter<matrix_t> row_iter_t;

  boost::scoped_array<double> arr(new double[1000]);
  external_storage_adaptor<double> storage(1000, arr.get());

  matrix_t mat(3, 4, storage);

  {
    row_iter_t it(mat, 0), ite(mat, 3);
    for(int i = 0; i < 3; i++, ++it)
    {
      BOOST_CHECK(it != ite);
      row_iter_t::value_type v(*it);

      int j = 0;
      for(row_iter_t::value_type::iterator it1(v.begin()), it1e(v.end()); it1 != it1e; ++it1, j++)
      {
        *it1 = i * 10 + j;
      }
    }
  }

  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 4; j++)
    {
      BOOST_CHECK_EQUAL(mat(i, j), i * 10 + j);
    }
  }


  // this is exactly the kind of access we want:
  *row_iter_t(mat, 0) = ub::zero_vector<double>(4);

  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 4; j++)
    {
      BOOST_CHECK_EQUAL(mat(i, j), i == 0 ? 0 : i * 10 + j);
    }
  }

}


// testing the row write access capability for the row proxy, vector access, external storage
BOOST_AUTO_TEST_CASE(test_row_proxy_random_access_iterator)
{
  using namespace grassmann_averages_pca::details::ublas_helpers;
  namespace ub = boost::numeric::ublas;
  
  typedef ub::matrix<double, ub::row_major > matrix_t;
  typedef row_iter<matrix_t> row_iter_t;

  matrix_t mat(20, 30);
  row_iter_t it(mat, 0), ite(mat, 3);
  BOOST_CHECK_EQUAL(ite-it, 3);

  BOOST_CHECK(
    (
      boost::is_same<
        std::iterator_traits<row_iter_t>::iterator_category,
        std::random_access_iterator_tag
      >::value
    ));
}

