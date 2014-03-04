
#include <boost/test/unit_test.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <include/private/boost_ublas_matrix_helper.hpp>

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
