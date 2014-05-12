// Copyright 2014, Max Planck Institute for Intelligent Systems.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

#ifndef ROBUST_PCA_STABLE_TRIMMING_HPP__
#define ROBUST_PCA_STABLE_TRIMMING_HPP__

/*!@file
 * Robust PCA functions, following the paper of Soren Hauberg.
 *
 * This file contains the implementation of the trimmed version, stable version.
 */

// boost heaps
#include <boost/heap/pairing_heap.hpp>

// for the thread pools
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/signals2.hpp>



// utilities
#include <include/private/utilities.hpp>


namespace robust_pca
{
  
  
  
  /*!@brief Robust PCA subspace algorithm, with stable trimming
   *
   * This class implements the robust PCA using the Grassmanian averaging. This is the naive implementation which is
   * suitable for small datasets.
   *
   * @author Soren Hauberg, Raffi Enficiaud
   */
  template <class data_t, class norm_mu_t = details::norm2>
  struct robust_pca_with_stable_trimming_impl
  {
  private:
    details::random_data_generator<data_t> random_init_op;
    norm_mu_t norm_op;
    double trimming_percentage;
    
    template <class vector_acc_t>
    void apply_quantile_to_vector(const data_t& current_data, vector_acc_t& v_acc) const
    {
      typename vector_acc_t::iterator it_acc_features(v_acc.begin());
      for(typename data_t::const_iterator it(current_data.begin()), ite(current_data.end()); it < ite; ++it, ++it_acc_features)
      {
        (*it_acc_features)(*it);
      }
    }
    
    template <class vector_bounds_t, class vector_number_elements_t>
    void selective_acc_to_vector(
       const vector_bounds_t& lower_bounds, const vector_bounds_t& upper_bounds,
       const data_t &initial_data,
       bool sign,
       data_t &v_selective_accumulator,
       vector_number_elements_t& v_selective_acc_count) const
    {
      for(int i = 0, j = initial_data.size(); i < j; i++)
      {
        typename data_t::value_type const v(initial_data[i]);
        if(v < lower_bounds[i])
          continue;
        if(v > upper_bounds[i])
          continue;
        
        if(sign)
        {
          v_selective_accumulator[i] += v;
        }
        else
        {
          v_selective_accumulator[i] -= v;
        }
        v_selective_acc_count[i]++;
      }
    }
    
    
    // the count does not change 
    template <class vector_quantile_t>
    void selective_update_acc_to_vector(
      const vector_quantile_t& quantiles1, const vector_quantile_t& quantiles2,
      const data_t &initial_data,
      bool sign,
      data_t &v_selective_accumulator) const
    {
      for(int i = 0, j = initial_data.size(); i < j; i++)
      {
        typename data_t::value_type const v(initial_data[i]);
        if(v < quantiles1[i])
          continue;
        if(v > quantiles2[i])
          continue;
        if(sign)
        {
          v_selective_accumulator[i] += 2*v;
        }
        else
        {
          v_selective_accumulator[i] -= 2*v;
        }
      }
    }    
    
    
  public:
    robust_pca_with_stable_trimming_impl(double trimming_percentage_ = 0) :
    random_init_op(details::fVerySmallButStillComputable, details::fVeryBigButStillComputable),
    trimming_percentage(trimming_percentage_)
    {
      assert(trimming_percentage_ >= 0 && trimming_percentage_ <= 1);
    }
    
    
    
    /*! Performs the computation of the current subspace on the elements given by the two iterators.
     *  @tparam it_t an input forward iterator to input vectors points. Each element pointed by the underlying iterator should be iterable and
     *   should provide a vector point.
     *  @tparam it_o_projected_vectors output forward iterator pointing on a container of vector points.
     *  @tparam it_norm_t an output iterator on weights/norms of the vectors. The output elements should be numerical (norm output)
     *
     * @param[in] max_iterations the maximum number of iterations at each dimension.
     * @param[in] max_dimension_to_compute the maximum number of dimensions to compute in the PCA (only the first max_dimension_to_compute will be
     *            computed).
     * @param[in] it input iterator at the beginning of the data
     * @param[in] ite input iterator at the end of the data
     * @param[in, out] it_norm_out input read-write iterator at the beginning of the computed norms. The iterator should be able to address
     *            as many element as there is in between it and ite (ie. @c std::distance(it, ite)
     * @param[in, out] it_projected
     * @param[in] initial_guess if provided, the initial vector will be initialized to this value.
     * @param[out] it_eigenvectors an iterator on the beginning of the area where the detected eigenvectors will be stored. The space should be at least max_dimension_to_compute.
     *
     * @returns true on success, false otherwise
     * @pre
     * - @c !(it >= ite)
     * - all the vectors given by the iterators pair should be of the same size (no check is performed).
     */
    template <class it_t, class it_o_projected_vectors, class it_o_eigenvalues_t>
    bool batch_process(
      const int max_iterations,
      int max_dimension_to_compute,
      it_t const it,
      it_t const ite,
      it_o_projected_vectors const it_projected,
      it_o_eigenvalues_t it_eigenvectors,
      std::vector<data_t> const * initial_guess = 0)
    {
      using namespace boost::accumulators;
      
      
      // add some log information
      if(it >= ite)
      {
        return false;
      }
      
      
      
      // estimation of the quantiles using the P^2 algorithm.
      typedef accumulator_set<double, stats< tag::min, tag::max, tag::extended_p_square_quantile(quadratic)> > accumulator_t;
      
      
      size_t size_data(0);
      
      const int number_of_dimensions = static_cast<int>(it->size());
      max_dimension_to_compute = std::min(max_dimension_to_compute, number_of_dimensions);
      
      
      
      
      // the first element is used for the init guess because for dynamic std::vector like element, the size is needed.
      data_t mu(initial_guess != 0 ? (*initial_guess)[0] : random_init_op(*it));
      assert(mu.size() == number_of_dimensions);
      
      // normalizing
      typename norm_mu_t::result_type norm_mu(norm_op(mu));
      mu /= norm_mu;
      
      
      
      // Accumulator object containing the required quantiles. 
      std::vector<double> probs = { trimming_percentage/2, 1-trimming_percentage/2 };
      const accumulator_t quantil_obj(extended_p_square_probabilities = probs);
      
      // vector of accumulator objects
      std::vector<accumulator_t> v_acc(number_of_dimensions, quantil_obj);
      
      // vector of valid bounds
      std::vector<double> v_min_threshold(number_of_dimensions);
      std::vector<double> v_max_threshold(number_of_dimensions);
      
      
      
      // first pass on the data, we compute the bounds
      v_acc.assign(number_of_dimensions, quantil_obj);
      
      
      // copy of the vectors. This should be avoided if the containers are of the same type and for the first iteration. 
      // The vectors are copied into the temporary container. During the copy, the percentiles are computed.
      it_o_projected_vectors it_tmp_projected(it_projected);
      for(it_t it_copy(it); it_copy != ite; ++it_copy, ++it_tmp_projected, size_data++)
      {
        *it_tmp_projected = *it_copy;
        apply_quantile_to_vector(*it_tmp_projected , v_acc);
      }
      
      // extracting the bounds
      for(int i = 0; i < number_of_dimensions; i++)
      {
        v_min_threshold[i] = quantile(v_acc[i], quantile_probability = trimming_percentage/2);
        v_max_threshold[i] = quantile(v_acc[i], quantile_probability = 1-trimming_percentage/2);
      }
      
      
      
      // initial iterator on the output eigenvectors
      it_o_eigenvalues_t const it_output_eigen_vector_beginning(it_eigenvectors);
      
      // storing the signs for fast updates
      std::vector<bool> v_signs(size_data);
      
      // for each dimension
      for(int current_dimension = 0; current_dimension < max_dimension_to_compute; current_dimension++, ++it_eigenvectors)
      {
        
        details::convergence_check<data_t> convergence_op(mu);
        
        data_t acc = data_t(number_of_dimensions, 0);
        std::vector<size_t> acc_counts(number_of_dimensions, 0);
        
        data_t previous_mu = mu;
        
        // fill the accumulator
        for(size_t s = 0; s < size_data; ++it_tmp_projected, s++)
        {
          typename it_o_projected_vectors::reference current_data = *it_tmp_projected;
          
          bool sign = boost::numeric::ublas::inner_prod(current_data, previous_mu) >= 0;
          v_signs[s] = sign;
          
          // trimming is applied to the accumulator + taking into account the sign.
          selective_acc_to_vector(v_min_threshold, v_max_threshold, current_data, sign, acc, acc_counts);
          
        }
        
        
        // compute the first mu
        for(int i = 0; i < number_of_dimensions; i++)
        {
          assert(acc_counts[i]);
          mu[i] = acc[i] / acc_counts[i];
        }
        mu /= norm_op(mu);
        
        
        // performs the iterations
        for(int iterations = 0; (!convergence_op(mu) && iterations < max_iterations); iterations++)
        {
          it_o_projected_vectors it_tmp_projected(it_projected);
          
          for(size_t s = 0; s < size_data; ++it_tmp_projected, s++)
          {
            typename it_o_projected_vectors::reference current_data = *it_tmp_projected;
            
            bool sign = boost::numeric::ublas::inner_prod(current_data, previous_mu) >= 0;
            
            if(sign != v_signs[s])
            {
              v_signs[s] = sign;
              selective_update_acc_to_vector(v_min_threshold, v_max_threshold, current_data, sign, acc);
              
            }
          }
          
          for(int i = 0; i < number_of_dimensions; i++)
          {
            assert(acc_counts[i]);
            mu[i] = acc[i] / acc_counts[i];
          }
          mu /= norm_op(mu);
          
          
        }
        
        // orthogonalisation 
        for(it_o_eigenvalues_t it_mus(it_output_eigen_vector_beginning); it_mus < it_eigenvectors; ++it_mus)
        {
          mu -= boost::numeric::ublas::inner_prod(mu, *it_mus) * (*it_mus);
        }
        mu /= norm_op(mu);
        
        // mu is the eigenvector of the current dimension, we store it in the output vector
        *it_eigenvectors = mu;
        
        // projection onto the orthogonal subspace
        if(current_dimension < max_dimension_to_compute - 1)
        {
          it_o_projected_vectors it_tmp_projected(it_projected);
          
          // update of vectors in the orthogonal space
          // for that to work, we need to reproject on the orthogonal subspace of all the previous eigenvalues
          for(size_t s(0); s < size_data; ++it_tmp_projected, s++)
          {
            typename it_o_projected_vectors::reference current_vector = *it_tmp_projected;
            current_vector -= boost::numeric::ublas::inner_prod(mu, current_vector) * mu;
          }
          
          mu = initial_guess != 0 ? (*initial_guess)[current_dimension + 1] : random_init_op(*it);
          
          
          // orthogonalisation 
          for(it_o_eigenvalues_t it_mus(it_output_eigen_vector_beginning); it_mus <= it_eigenvectors; ++it_mus)
          {
            mu -= boost::numeric::ublas::inner_prod(mu, *it_mus) * (*it_mus);
          }
          mu /= norm_op(mu);
          
          
          
        }
        
        
        
        
      }
      
      return true;
    }
  };
  
  
  
  


}


#endif /* ROBUST_PCA_STABLE_TRIMMING_HPP__ */
