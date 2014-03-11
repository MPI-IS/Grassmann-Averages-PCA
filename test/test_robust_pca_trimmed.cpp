
/*!@file
 * This file includes the tests for the trimmed version of the robust pca
 */

#include <boost/test/unit_test.hpp>
#include <test/test_main.hpp>

#include <include/robust_pca.hpp>
#include <include/private/boost_ublas_matrix_helper.hpp>

// data stored into a matrix
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

// generating data randomly
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <fstream>



BOOST_FIXTURE_TEST_SUITE(basic_checks, fixture_simple_matrix_creation)




BOOST_AUTO_TEST_SUITE_END();
