
#include <boost/test/unit_test.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <include/private/boost_ublas_matrix_helper.hpp>
#include <include/private/boost_ublas_matlab_helper.hpp>

#include <boost/scoped_array.hpp>


BOOST_TEST_DONT_PRINT_LOG_VALUE(robust_pca::ublas_adaptor::row_iter<boost::numeric::ublas::matrix<double> >);
BOOST_TEST_DONT_PRINT_LOG_VALUE(robust_pca::ublas_adaptor::row_iter<const boost::numeric::ublas::matrix<double> >);

BOOST_AUTO_TEST_CASE(test_row_proxy_read)
{
  using namespace robust_pca::ublas_adaptor;
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


BOOST_AUTO_TEST_CASE(test_row_proxy_write)
{
  using namespace robust_pca::ublas_adaptor;
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
  using namespace robust_pca::ublas_adaptor;
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

BOOST_AUTO_TEST_CASE(test_row_proxy_write_specific_storage)
{
  using namespace robust_pca::ublas_adaptor;
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



BOOST_AUTO_TEST_CASE(test_row_proxy_write_specific_storage_row_access)
{
  using namespace robust_pca::ublas_adaptor;
  namespace ub = boost::numeric::ublas;

  typedef ub::matrix<double, ub::row_major, ub::bounded_array<double, 1000> > matrix_t;
  typedef row_iter<matrix_t> row_iter_t;

  matrix_t mat(3, 4);

  {
    row_iter_t it(mat, 0), ite(mat, 3);
    for(int i = 0; i < 3; i++, ++it)
    {
      BOOST_CHECK(it != ite);
      
      int j = 0;
      *it = ub::zero_vector<double>(3);
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


namespace boost {
  namespace numeric {
    namespace ublas {


      /// For the creation of temporary vectors in the assignment of proxies
      template <class T, class L>
      struct vector_temporary_traits<matrix<T, L, robust_pca::ublas_matlab_helper::external_storage_adaptor<T> > >
      {
        typedef vector<T> type;
      };

    }
  }
}

BOOST_AUTO_TEST_CASE(test_row_proxy_write_specific_external_storage)
{
  using namespace robust_pca::ublas_adaptor;
  using namespace robust_pca::ublas_matlab_helper;
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
  *row_iter_t(mat, 0) = ub::zero_vector<double>(3);


}

