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
  struct mu_reference
  {
    size_t count;         //!< the number of vectors referencing the mu
    size_t index_matrix;  //!< the index of the mu in the matrix
    size_t iteration_index;    //!< the index of the mu in the list @todo check if relevant
    data_t mu;            //!< the data
    
    mu_reference(const data_t& mu_) : mu(mu_) {}
  };

  typedef typename data_t::value_type scalar_t;

  // we need to keep the mus in order to compute the inner product
  typedef std::list<mu_reference> mu_container_t;
  mu_container_t mus;
  
  // This will map the index of the algorithm to the index of the list above
  typedef std::map<size_t, size_t> map_index_to_list_index_t;
  map_index_to_list_index_t mapping_indices;
  
  typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower> mu_angles_t;
  mu_angles_t mu_angles;
  




  mus_proximity() : mus(), mu_angles()
  {}
  
  // may be optimized
  double angle(const data_t& left, const data_t& right) const
  {
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
    std::cout << " -- before acos " << out / (sqrt(norm_l) * sqrt(norm_r)) << std::endl;
    std::cout << " -- bigger than one? " << (out / (sqrt(norm_l) * sqrt(norm_r)) > 1 ? "yes" : "no") << std::endl;
    return acos(out / (sqrt(norm_l) * sqrt(norm_r)));
  }
  
  // This version supposes the vectors are already normalized
  // and the corresponding computations are simplified compared to @ref angle
  double angle_normalized_vectors(const data_t& left, const data_t& right) const
  {
    typedef typename data_t::value_type scalar_t;
    double out(0);
    
    scalar_t const *const left_p = &left.data()[0];
    scalar_t const *const right_p= &right.data()[0];
    for(size_t i = 0, j = left.size(); i < j; i++)
    {
      out += left_p[i] * right_p[i];
    }
    std::cout << " -- before acos " << out << std::endl;
    std::cout << " -- bigger than one? " << (out > 1 ? "yes" : "no") << std::endl;
    return acos(out);  
  }
  
  //! Removes any unnecessary element from the list of mus
  //! in order to lower the computation of the symmetric matrix to its bare minimum
  //! @return the number of pruned elements
  size_t prune()
  {
    std::set<size_t> index_to_remove;

    //
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

    return index_to_remove.size();
  }


  //! Returns the number of vectors that are currently held by the instance
  size_t get_nb_mus() const
  {
    return mu_angles.size1();
  }

  mu_angles_t const& get_angle_matrix() const
  {
    return mu_angles;
  }

  mu_container_t const& get_mus() const
  {
    return mus;
  }
  
  //! Updates the count of a particular mu
  //!
  //! @param[in] delta_count positive count means added reference, while negative count means
  //!            removed reference
  void update_count(size_t iteration_index, int delta_count)
  {
    assert(mapping_indices.count(iteration_index) > 0);
    size_t index_in_list = mapping_indices[iteration_index];
    
    typename mu_container_t::iterator it(mus.begin());
    std::advance(it, index_in_list);
    it->count += delta_count;
  }


  //! Adds a mu to the managed list of mus
  //!
  //! The mu is stored in the first available place
  mu_reference const& add_mu(const data_t& new_mu, size_t iteration_index, size_t nb_references)
  {
    using namespace std; // bringing acos
    
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
      // the position in the matrix stays the same
    }
    else
    {
      // there is nothing to prune here anyway, so we add a new column/row to the 
      // matrix
      
      mus.push_back(mu_reference(new_mu));
      it = --mus.end();
      it->index_matrix = mu_angles.size1();
      it->iteration_index   = iteration_index;
      index_in_list    = mus.size()-1;
      it->count = nb_references;      
      
      // adding a row / column
      std::cout << "before" << std::endl;
      print_matrix(this->get_angle_matrix());
      mu_angles.resize(mu_angles.size1() + 1, true);
      mu_angles(mu_angles.size1()-1, mu_angles.size1()-1) = 0;
      std::cout << "after" << std::endl;
      print_matrix(this->get_angle_matrix());
            
    }
    
    mapping_indices[iteration_index] = index_in_list;
    
    
    for(typename mu_container_t::const_iterator it2(mus.begin()); it != ite; ++it)
    {
      // computing the inner product and storing into the appropriate place of
      // the "distance" matrix
      // we avoid the computation of these angles for the elements that will get pruned
    
      if(it2 == it)
        continue;
        
      if(it2->count == 0)
        continue;
    
      mu_angles(it2->index_matrix, it->index_matrix) = angle(it->mu, it2->mu);
      std::cout << " -- " << mu_angles(it2->index_matrix, it->index_matrix) << std::endl;
    }
    
    return *it;
  }

};




BOOST_AUTO_TEST_CASE(test_add_mu)
{
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

    instance.add_mu(current, i+1, 10);
  }
  
  print_matrix(instance.get_angle_matrix());
  
}


#if 0



  //std::cout << instance.get_angle_matrix() << std::endl;
  BOOST_CHECK_EQUAL(instance.prune(), 1);

  // checking the removal
  mu_angle_helper_t::mu_container_t const& mus = instance.get_mus();

  BOOST_CHECK_EQUAL(mus.size(), nb_elements-1);


  // checking the first element
  size_t value = 0;
  size_t index = 0;
  for(mu_angle_helper_t::mu_container_t::const_iterator it(mus.begin()), ite(mus.end()); 
      it != ite;
      ++it, index++)
  {
    BOOST_CHECK_EQUAL(it->mu(0), value);
    if(index != 3)
    {
      value++;
    }
  }

  mu_angle_helper_t::mu_angles_t const& mu_angles = instance.get_angle_matrix();
  BOOST_CHECK_EQUAL(mu_angles.size1(), nb_elements-1);
  BOOST_CHECK_EQUAL(mu_angles.size2(), nb_elements-1);


  // checking the angles
  value = 0;
  index = 0;
  for(size_t i = 0, to_skip = 0; i < instance.get_nb_mus(); i++)
  {
    for(size_t j = i, to_skip = 0; i < instance.get_nb_mus(); i++)
    {
      //instance.get_angle_matrix()(i, 0) = 
    }
  }

}
#endif