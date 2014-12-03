// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <test/test_main.hpp>


// the random number generator
boost::random::mt19937 rng;

const int fixture_simple_matrix_creation::nb_elements = 10000;
const int fixture_simple_matrix_creation::dimensions = DATA_DIMENSION;
