// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

/*!@file
 * This file contains the tests for the proximities of mus in order to avoid
 * some of the computations of the inner products
 */

#include <boost/test/unit_test.hpp>
#include <test/test_main.hpp>


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <set>

template <class matrix_t>
void print_matrix(const matrix_t& mat)
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  for(size_t i = 0; i < mat.size1(); i++)
  {
    for(size_t j = 0; j < mat.size2(); j++)
    {
      std::cout << std::fixed << mat(i, j) << " ";
    }
    std::cout << std::endl;
  }

}


BOOST_AUTO_TEST_CASE(test_add_mu_all_similar_gives_null_matrix)
{
  // mainly testing some numerical stability of the computations, as acos is surprisingly sensitive
  namespace ub = boost::numeric::ublas;
  namespace ga = grassmann_averages_pca;

  typedef ub::vector<double> data_t;
  typedef ga::details::mus_distance_optimisation<data_t> mu_angle_helper_t;
  mu_angle_helper_t instance;

  const int nb_elements = 3;
  const int dimensions = 5;

  for(int i = 0; i < nb_elements; i++)
  {
    data_t current(dimensions);
    for(int j = 0; j < dimensions; j++)
    {
      current(j) = (i + 1);
    }
    
    current /= ga::details::norm2()(current);

    instance.add_mu(current, i, 10);
  }

  // the max of the abs should be close to 0 since all vectors are the same
  double norm_matrix = ub::norm_1(instance.get_angle_matrix());
  BOOST_CHECK_SMALL(norm_matrix, 1E-5);

  
  //print_matrix(instance.get_angle_matrix());
  
}


BOOST_AUTO_TEST_CASE(test_add_mu_update_first_available)
{
  namespace ub = boost::numeric::ublas;
  namespace ga = grassmann_averages_pca;

  typedef ub::vector<double> data_t;
  typedef ga::details::mus_distance_optimisation<data_t> mu_angle_helper_t;
  mu_angle_helper_t instance;

  const int nb_elements = 3;
  const int dimensions = 5;

  int i = 0;
  for(; i < nb_elements; i++)
  {
    data_t current(dimensions);
    for(int j = 0; j < dimensions; j++)
    {
      current(j) = (i + j + 1) % 17;
    }
    
    current /= ga::details::norm2()(current);

    instance.add_mu(current, i, 10);
  }
  //print_matrix(instance.get_angle_matrix());

  // should not prune
  instance.update_count(2, -9);

  {
    data_t current(dimensions);
    for(int j = 0; j < dimensions; j++)
    {
      current(j) = (i + j + 1) % 17;
    }
    
    current /= ga::details::norm2()(current);

    instance.add_mu(current, i, 10);
  }
  i++;

  BOOST_CHECK_EQUAL(instance.get_mus().size(), nb_elements + 1);
  //print_matrix(instance.get_angle_matrix());
  const double first_element_column_added = instance.get_angle_matrix()(0, instance.get_angle_matrix().size1() - 1);
  const double first_element_column2 = instance.get_angle_matrix()(0, 2);



  // should prune
  instance.update_count(2, -1);
  {
    data_t current(dimensions);
    for(int j = 0; j < dimensions; j++)
    {
      current(j) = (i + j + 1) % 17;
    }
    
    current /= ga::details::norm2()(current);

    instance.add_mu(current, i, 10);
  }
  i++;

  BOOST_CHECK_EQUAL(instance.get_mus().size(), nb_elements + 1);
  BOOST_CHECK_EQUAL(first_element_column_added, instance.get_angle_matrix()(0, instance.get_angle_matrix().size1() - 1));
  BOOST_CHECK_NE(first_element_column2, instance.get_angle_matrix()(0, 2));



  
  //print_matrix(instance.get_angle_matrix());
  
}


BOOST_AUTO_TEST_CASE(test_prune)
{
  namespace ub = boost::numeric::ublas;
  namespace ga = grassmann_averages_pca;

  typedef ub::vector<double> data_t;
  typedef ga::details::mus_distance_optimisation<data_t> mu_angle_helper_t;
  mu_angle_helper_t instance;
  mu_angle_helper_t instance2;

  const int nb_elements = 5;
  const int dimensions = 5;

  int table_count[nb_elements] = {0,1,13,10,10};

  
  for(int i = 0; i < nb_elements; i++)
  {
    data_t current(dimensions);
    for(int j = 0; j < dimensions; j++)
    {
      current(j) = (i + j + 1) % 17;
    }
    
    current /= ga::details::norm2()(current);

    instance.add_mu(current, i, 10);
    if(table_count[i])
    {
      instance2.add_mu(current, i, 10);
    }
  }
  //print_matrix(instance.get_angle_matrix());

  for(int i = 0; i <  nb_elements; i++)
  {
    instance.update_count(i, table_count[i] - 10);
  }

  double angle_between_indices_3_and_4 = instance.get_angle_from_indices(3, 4);
  //std::cout << "angle_between_indices_3_and_4 = " << angle_between_indices_3_and_4 << std::endl;
  //print_matrix(instance.get_angle_matrix());

  // only one gets pruned
  BOOST_CHECK_EQUAL(instance.prune(), 1);
  BOOST_CHECK_EQUAL(instance.get_angle_matrix().size1(), nb_elements - 1);
  BOOST_REQUIRE_EQUAL(instance.get_angle_matrix().size1(), instance2.get_angle_matrix().size1());

  double norm_difference = ub::norm_1(instance.get_angle_matrix() - instance2.get_angle_matrix());
  BOOST_CHECK_SMALL(norm_difference, 1E-8);


  // checks the mapping of the indexes are correct after a prune
  instance.update_count(1, -1); // should be 0 now
  BOOST_CHECK_EQUAL(instance.prune(), 1);

  //std::cout << "angle_between_indices_3_and_4 = " << instance.get_angle_from_indices(3, 4) << std::endl;

  // this vectors did not move, so the angle should remain equal
  BOOST_CHECK_EQUAL(angle_between_indices_3_and_4, instance.get_angle_from_indices(3, 4));

  //print_matrix(instance.get_angle_matrix());


  // now should add the new element at position 4
  instance.update_count(2, -13);
  {
    data_t current(dimensions);
    for(int j = 0; j < dimensions; j++)
    {
      current(j) = (nb_elements + j + 1) % 17;
    }
    
    current /= ga::details::norm2()(current);

    instance.add_mu(current, nb_elements, 10);
  }
  
  // this vectors did not move, so the angle should remain equal
  BOOST_CHECK_EQUAL(angle_between_indices_3_and_4, instance.get_angle_from_indices(3, 4));


}

