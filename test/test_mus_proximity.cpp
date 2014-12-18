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


/*!brief Handles the changes between two consecutive estimations of mus.
 *
 */
template <class data_t>
struct mus_proximity
{
  struct mu_reference
  {
    size_t count; //!< the number of vectors referencing the mu
    size_t index; //!< the index of the mu in the matrix
    data_t mu;    //!< the data
    
    mu_reference(const data_t& mu_) : mu(mu_) {}
  };

  std::list<mu_reference> mus;
  typedef typename std::list<mu_reference>::iterator it_t;
  
  typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower> mu_angles_t;
  
  mu_angles_t mu_angles;
  
  mus_proximity() : mus(), mu_angles()
  {}
  
  // may be optimized
  double inner_product(const data_t& left, const data_t& right) const
  {
    typedef typename data_t::value_type scalar_t;
    double out(0);
    
    scalar_t const *const left_p = left.data()[0];
    scalar_t const *const right_p= right.data()[0];
    for(size_t i = 0, j = left.size(); i < j; i++)
    {
      out += left_p[i] * right_p[i];
    }
    return out;
  }
  
  //! Removes any unnecessary element from the list of mus
  //! in order to lower the computation of the symmetric matrix to its bare minimum
  void prune()
  {
    std::set<size_t> index_to_remove;

    {
      size_t index(0);
      for(it_t it(mus.begin()), ite(mus.end()); it != ite; ++it, index++)
      {
        if(!it->count)
        {
          index_to_remove.insert(index);
        }
      }
    }
    
    mu_angles_t new_matrix(mu_angles.size1() - index_to_remove.size());
    
    for(size_t index = 0, index2 = 0; index < mu_angles.size1(); index++)
    {
      if(index_to_remove.count(index))
        continue;
      
      
      index2++;
      
    }
    
    
    new_matrix.swap(mu_angles);
  }

  void add_mu(const data_t& new_mu)
  {
    using namespace std; // bringing acos
    
    mus.push_back(mu_reference(new_mu));
    mu_reference &last = mus.back();
    last.index = mu_angles.size1();
    last.count = 0; // todo: count
    
    
    
    for(it_t it(mus.begin()), ite(mus.end()); it != ite; ++it)
    {
      // computing the inner product and storing into the appropriate place of
      // the "distance" matrix
      double inner = inner_product(last.mu, it->mu);
      mu_angles(it->index, last.index) = acos(inner);
    }
  }

};