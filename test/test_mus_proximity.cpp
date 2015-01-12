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

#include <boost/numeric/ublas/symmetric.hpp>
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


/*!brief Handles the changes between two consecutive estimations of mus.
 *
 */
template <class data_t>
struct mus_proximity
{

private:
  //! Structure intended to store the mus and the number 
  //! of vectors referencing them. 
  struct mu_reference_count_data
  {
    size_t count;               //!< the number of vectors referencing the mu
    size_t iteration_index;     //!< the index of the iteration on which this mu appeared
    data_t mu;                  //!< the data
    
    mu_reference_count_data(const data_t& mu_) : mu(mu_) {}
  };

  typedef typename data_t::value_type scalar_t;

  // we need to keep the mus for computing the inner product
  typedef std::list<mu_reference_count_data> mu_container_t;
  mu_container_t mus;
  
  // This will map the index of the mu from the algorithm to the index of the list @c mus
  typedef std::map<size_t, size_t> map_index_to_list_index_t;
  map_index_to_list_index_t mapping_indices;
  
  // The matrix containing the angles between the vectors stored in the mu container. 
  // This is a symetric matrix and the indices follow the one of the list @c mus
  typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower> mu_angles_t;
  mu_angles_t mu_angles;
 

public:

  mus_proximity() : mus(), mapping_indices(), mu_angles()
  {}
  
  //! Computes the angles between two vectors. 
  //!
  //! This function is meant to be used by general vectors (non normalized). If the vectors are 
  //! normalised already, consider using angle_normalized_vectors. 
  double angle(const data_t& left, const data_t& right) const
  {
    using namespace std; // bringing acos

    typedef typename data_t::value_type scalar_t;
    double out(0);
    double norm_l(0), norm_r(0);
    
    scalar_t const *const left_p = &left.data()[0];
    scalar_t const *const right_p= &right.data()[0];
    for(size_t i = 0, j = left.size(); i < j; i++)
    {
      out += left_p[i] * right_p[i];
      norm_l += left_p[i]*left_p[i];
      norm_r += right_p[i]*right_p[i];
    }
    
    if(norm_l < 1E-10 || norm_r < 1E-10)
      return acos(0); // pi/2

    out /= sqrt(norm_l) * sqrt(norm_r);
    if(out > 1) 
    {
      out = 1;
    }
    else if(out < -1)
    {
      out = -1;
    }
    return acos(out);
  }
  
  //! Computes the angles between two vectors. 
  //!
  //! This version supposes the vectors are already normalized
  //! and the corresponding computations are simplified compared to @ref angle
  double angle_normalized_vectors(const data_t& left, const data_t& right) const
  {
    using namespace std; // bringing acos

    typedef typename data_t::value_type scalar_t;
    double out(0);
    
    scalar_t const *const left_p = &left.data()[0];
    scalar_t const *const right_p= &right.data()[0];
    for(size_t i = 0, j = left.size(); i < j; i++)
    {
      out += left_p[i] * right_p[i];
    }

    if(out > 1) 
    {
      out = 1;
    }
    else if(out < -1)
    {
      out = -1;
    }

    return acos(out);  
  }
  
  //! Removes any unnecessary element from the list of mus
  //! in order to lower the computation of the symmetric matrix to its bare minimum
  //! @return the number of pruned elements
  size_t prune()
  {
    std::set<size_t> index_to_remove;

    // identifies the rows/columns to be pruned
    {
      size_t index(0);
      for(typename mu_container_t::iterator it(mus.begin()), ite(mus.end()); it != ite; ++it, index++)
      {
        if(!it->count)
        {
          index_to_remove.insert(index);
          it = mus.erase(it);
        }
      }
    }
    

    // reconstructs the matrix
    mu_angles_t new_matrix(mu_angles.size1() - index_to_remove.size());
    
    for(size_t index = 0, index_to_copy_to = 0; index < mu_angles.size1(); index++)
    {
      if(index_to_remove.count(index))
        continue;
      
      for(size_t index2 = index, index_to_copy_to2 = index_to_copy_to; index2 < mu_angles.size1(); index2++)
      {
        if(index_to_remove.count(index2))
          continue;
        
        new_matrix(index_to_copy_to, index_to_copy_to2) = mu_angles(index, index2);
        index_to_copy_to2++;
      }
      
      
      index_to_copy_to++;
      
    }
    
    new_matrix.swap(mu_angles);


    // reconstruct the mapping index
    mapping_indices.clear();
    {
      size_t index(0);
      for(typename mu_container_t::iterator it(mus.begin()), ite(mus.end()); it != ite; ++it, index++)
      {
        mapping_indices[it->iteration_index] = index;
      }
    }

    return index_to_remove.size();
  }


  //! Returns the number of vectors that are currently held by the instance
  size_t get_nb_mus() const
  {
    return mu_angles.size1();
  }

  //! Returns the angles matrix. 
  //!
  //! This should not be used for accessing the angles directly as the indices of the elements
  //! do not follow the indices of the algorithm. The mapping provided by mu_reference_count_data::iteration_index
  //! should be used jointly to the access to this matrix.
  mu_angles_t const& get_angle_matrix() const
  {
    return mu_angles;
  }

  //! Returns the list of elements that are stored into this container. 
  mu_container_t const& get_mus() const
  {
    return mus;
  }


  //! Gets the value of the angle between two iteration index vectors
  double get_angle_from_indices(size_t index1, size_t index2) const
  {
    assert(mapping_indices.count(index1) > 0);
    assert(mapping_indices.count(index2) > 0);
    return mu_angles(mapping_indices.at(index1), mapping_indices.at(index2));
  }
  
  //! Updates the count of a particular mu
  //!
  //! @param[in] iteration_index index manipulated by the algorithm, indicating the time of arrival of the mu that
  //!            is targetted by the update.
  //! @param[in] delta_count positive count means added reference, while negative count means
  //!            removed reference. 
  //! 
  //! @pre 
  //! This value of delta_count + current_count should be positive.
  void update_count(size_t iteration_index, int delta_count)
  {
    assert(mapping_indices.count(iteration_index) > 0);
    size_t index_in_list = mapping_indices[iteration_index];
    
    typename mu_container_t::iterator it(mus.begin());
    std::advance(it, index_in_list);
    assert(it->count + delta_count >= 0);
    it->count += delta_count;
  }


  //! Adds a mu to the managed list of mus
  //!
  //! The mu is stored in the first available place. 
  //!
  //! @param[in] new_mu new value of mu to be stored
  //! @param[in] iteration_index iteration on the algorithm side at which new_mu arrived. This is the value that should
  //!            be consecutevely used for accessing the angles. 
  //! @param[in] nb_references number of references to this mu
  //! 
  //! @note 
  //! It is possible to have a nb_references equal to 0 since this is what happens for a newly computed mu. 
  mu_reference_count_data const& add_mu(const data_t& new_mu, size_t iteration_index, size_t nb_references)
  {
    
    // finding the first available place to put this new mu:
    size_t index_in_list(0);
    typename mu_container_t::iterator it(mus.begin());
    typename mu_container_t::iterator ite(mus.end());
    
    for(; it != ite; ++it, index_in_list++)
    {
      if(!it->count)
      {
        break;
      }
    }
    
    if(it != ite)
    {
      it->mu = new_mu;
      it->count= nb_references;
      it->iteration_index = iteration_index;
    }
    else
    {
      // there is nothing to prune here anyway, so we add a new column/row to the 
      // matrix
      
      mus.push_back(mu_reference_count_data(new_mu));
      it = --mus.end();
      it->iteration_index   = iteration_index;
      index_in_list    = mus.size()-1;
      it->count = nb_references;      
      
      // adding a row / column with preservation of the content
      mu_angles.resize(mu_angles.size1() + 1, true);
      mu_angles(mu_angles.size1()-1, mu_angles.size1()-1) = 0;
            
    }
    
    mapping_indices[iteration_index] = index_in_list;
    
    
    size_t current_index(0);
    for(typename mu_container_t::const_iterator it2(mus.begin()); it2 != ite; ++it2, current_index++)
    {
      // computing the inner product and storing into the appropriate place of
      // the "distance" matrix
    
      if(it2 == it)
        continue;
       
      // not updating against the one that will be pruned
      if(it2->count == 0)
        continue;
    
      mu_angles(current_index, index_in_list) = angle(it->mu, it2->mu);
    }
    
    return *it;
  }

};




BOOST_AUTO_TEST_CASE(test_add_mu_all_similar_gives_null_matrix)
{
  // mainly testing some numerical stability of the computations, as acos is surprisingly sensitive
  namespace ub = boost::numeric::ublas;
  namespace ga = grassmann_averages_pca;

  typedef ub::vector<double> data_t;
  typedef mus_proximity<data_t> mu_angle_helper_t;
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
  typedef mus_proximity<data_t> mu_angle_helper_t;
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
  typedef mus_proximity<data_t> mu_angle_helper_t;
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

