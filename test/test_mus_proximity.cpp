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
    size_t index_list;    //!< the index of the mu in the list
    data_t mu;            //!< the data
    
    mu_reference(const data_t& mu_) : mu(mu_) {}
  };

  typedef typename data_t::value_type scalar_t;

  // we need to keep the mus in order to compute the inner product
  typedef std::list<mu_reference> mu_container_t;
  mu_container_t mus;
  
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
    return acos(out / (sqrt(norm_l) * sqrt(norm_r)));
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
      for(typename mu_container_t::const_iterator it(mus.begin()), ite(mus.end()); it != ite; ++it, index++)
      {
        if(!it->count)
        {
          index_to_remove.insert(index);
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


  void add_mu(const data_t& new_mu, size_t nb_references)
  {
    using namespace std; // bringing acos
    
    mus.push_back(mu_reference(new_mu));
    mu_reference &last = mus.back();
    last.index_matrix = mu_angles.size1();
    last.index_list   = mus.size()-1;
    last.count = nb_references;
    
    
    // here we can test if there is something to prune, if so, replace the 
    // first available column by this vector and do not forget to remove it from
    // the list
    
    
    // adding a row / column
    mu_angles.resize(mu_angles.size1() + 1, true);
    mu_angles(mu_angles.size1()-1, mu_angles.size1()-1) = 0;
    
    typename mu_container_t::iterator it(mus.begin()), ite(mus.end());
    --ite;
    for(; it != ite; ++it)
    {
      // computing the inner product and storing into the appropriate place of
      // the "distance" matrix
      mu_angles(it->index_matrix, last.index_matrix) = angle(last.mu, it->mu);
    }
  }

};


BOOST_AUTO_TEST_CASE(test_prune)
{
  namespace ub = boost::numeric::ublas;

  typedef ub::vector<double> data_t;
  typedef mus_proximity<data_t> mu_angle_helper_t;
  mu_angle_helper_t instance;

  const int nb_elements = 10;
  const int dimensions = 5;

  for(int i = 0; i < nb_elements; i++)
  {
    data_t current(dimensions);
    for(int j = 0; j < dimensions; j++)
    {
      current(j) = i;
    }

    instance.add_mu(current, (i + 3) % nb_elements); // the 3rd should be pruned
  }


  std::cout << instance.get_angle_matrix() << std::endl;
  BOOST_CHECK_EQUAL(instance.prune(), 1);

  // checking the removal
  mu_angle_helper_t::mu_container_t const& mus = instance.get_mus();

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


  // checking the inner product
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