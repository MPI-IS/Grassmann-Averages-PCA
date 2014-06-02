// Copyright 2014, Max Planck Institute for Intelligent Systems.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

#ifndef ROBUST_PCA_HPP__
#define ROBUST_PCA_HPP__

/*!@file
 * Robust PCA functions, following the paper of Soren Hauberg.
 *
 * @note These implementations assume the existence of boost somewhere. 
 */

#include <vector>


#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector.hpp>


// for the thread pools
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/signals2.hpp>

// utilities
#include <include/private/utilities.hpp>


namespace robust_pca
{

  namespace ub = boost::numeric::ublas;

  /*!@brief Robust PCA with Grassmann Average algorithm
   *
   * This class implements the robust PCA using the Grassmann average. Its purpose is to compute the PCA of a dataset @f$\{X_i\}@f$, where each @f$X_i@f$ is a vector of dimension
   * D. The particularity of the Grassman average scheme is to be more stable than other algorithms.
   * 
   * The algorithm is the following:
   * - pick a random or a given @f$\mu_{k, 0}@f$, where @f$k@f$ is the current eigen-vector being computed and @f$0@f$ is the current iteration number (0). 
   * - until the sequence @f$(\mu_{k, t})_t@f$ converges, do:
   *   - computes the sign @f$s_{j, t}@f$ of the projection of the input vectors @f$X_j@f$ onto @f$\mu_{i, t}@f$. We have @f[s_{j, t} = X_j \cdot \mu_{k, t} \geq 0@f]
   *   - compute the update of @f$\mu{k, .}@f$: @f[\mu_{k, t+1} = \frac{\sum_j s_{j, t} X_j}{\left\|\sum_j s_{j, t} X_j\right\|}@f]
   * - project the @f$X_j@f$'s onto the orthogonal subspace of @f$\mu_{k} = \lim_{t \rightarrow +\infty} \mu_{k, t}@f$: @f[\forall j, X_{j} = X_{j} - X_{j}\cdot\mu_{k} @f]
   *
   * The range taken by @f$k@f$ is a parameter of the algorithm: @c max_dimension_to_compute (see robust_pca_impl::batch_process). 
   * The range taken by @f$t@f$ is also a parameter of the algorithm: @c max_iterations (see robust_pca_impl::batch_process).
   * The test for convergence is delegated to the class details::convergence_check.
   *
   * The computation is distributed among several threads. The multithreading strategy is 
   * - to split the computation of @f$\sum_j s_{j, t} X_j@f$ among several independant chunks. This computation involves the inner product and the sign. Each chunk addresses 
   *   a subset of the data @f$\{X_j\}@f$ without any overlap with other chunks. The maximal size of a chunk can be configured through the function robust_pca_impl::set_max_chunk_size.
   *   By default, the size of the chunk would be the size of the data divided by the number of threads.
   * - to split the computation of the projection onto the orthogonal subspace of @f$\mu_{k}@f$.
   * - to split the computation of the regular PCA algorithm (if any) into several independant chunks.
   *
   * The number of threads can be configured through the function robust_pca_impl::set_nb_processors.
   * 
   * @note
   * The algorithm may also perform a few "regular PCA" steps, which is the computation of the eigen-vector with highest eigen-value. This can be configured through the function
   * robust_pca_impl::set_nb_steps_pca.
   *
   * @tparam data_t type of vectors used for the computation. 
   * @tparam norm_mu_t norm used to normalize the eigen-vector and project them onto the unit circle.
   *
   * @author Soren Hauberg, Raffi Enficiaud
   */
  template <class data_t, class norm_mu_t = details::norm2>
  struct robust_pca_impl
  {
  private:
    //! Random generator for initialising @f$\mu@f$ at each dimension. 
    details::random_data_generator<data_t> random_init_op;

    //! Norm used for normalizing @f$\mu@f$.
    norm_mu_t norm_op;

    //! Number of parallel tasks that will be used for computing.
    size_t nb_processors;

    //! Maximal size of a chunk (infinity by default).
    size_t max_chunk_size;

    //! Number of steps for the initial PCA like algorithm (default to 3).
    size_t nb_steps_pca;

    //!@internal
    //!@brief Contains the logic for processing part of the accumulator
    struct asynchronous_chunks_processor
    {
    private:
      size_t nb_elements;
      data_t accumulator;
      std::vector<bool> v_signs;
      size_t data_dimension;
      
      // this is to send an update of the value of mu to all listeners
      // the connexion should be managed externally
      typedef boost::signals2::signal<void (data_t const&)>
        connector_accumulator_t;
      connector_accumulator_t signal_acc;

      typedef boost::signals2::signal<void ()>
        connector_counter_t;
      connector_counter_t signal_counter;


      
      typedef typename data_t::value_type scalar_t;
      
      //! The matrix containing a copy of the data
      scalar_t *p_c_matrix;
      
      std::vector<double> inner_prod_results;
    public:
      asynchronous_chunks_processor() : nb_elements(0), data_dimension(0), p_c_matrix(0)
      {
      }
      
      ~asynchronous_chunks_processor()
      {
        delete [] p_c_matrix;
      }


      //! Sets the data range
      template <class container_iterator_t>
      bool set_data_range(container_iterator_t const &b, container_iterator_t const& e)
      {
        //begin = b;
        //end = e;
        nb_elements = std::distance(b, e);
        assert(nb_elements > 0);
        v_signs.resize(nb_elements);
        inner_prod_results.resize(nb_elements);
        
        delete [] p_c_matrix;
        p_c_matrix = new scalar_t[nb_elements*data_dimension];
        
        container_iterator_t bb(b);

        for(int column = 0; column < nb_elements; column++, ++bb)
        {
          double* current_line = p_c_matrix + column;
          for(int line = 0; line < data_dimension; line ++, current_line += nb_elements)
          {         
            *current_line = (*bb)(line);
          }
          
        }
        
        return true;
      }

      //! Sets the dimension of the problem
      //! @pre data_dimensions_ strictly positive
      void set_data_dimensions(size_t data_dimensions_)
      {
        data_dimension = data_dimensions_;
        assert(data_dimension > 0);
      }

      //! Returns the connected object that will receive the notification of the updated accumulator.
      connector_accumulator_t& connector_accumulator()
      {
        return signal_acc;
      }

      //! Returns the connected object that will receive the notification of the end of the current process.
      connector_counter_t& connector_counter()
      {
        return signal_counter;
      }
      
      //! Inner product computation
      void compute_inner_products(data_t const &mu)
      {
        double *out = &inner_prod_results[0];
        double mu_element = mu(0);
        double *current_line = p_c_matrix;
        
        for(int column = 0; column < nb_elements; column++)
        {
          out[column] = mu_element * current_line[column];
        }
        current_line += nb_elements;
        
        for(int line = 1; line < data_dimension; line ++, current_line += nb_elements)
        {
          mu_element = mu(line);    
          for(int column = 0; column < nb_elements; column++)
          {
            out[column] += mu_element * current_line[column];
          }               
        }
      }

      //! PCA steps
      void pca_accumulation(data_t const &mu)
      {
        accumulator = data_t(data_dimension, 0);
        compute_inner_products(mu);
               
        double const * const p_inner_product = &inner_prod_results[0];

        typename data_t::iterator it(accumulator.begin());
        for(size_t i = 0; i < data_dimension; i++, ++it)
        {
          double const * const current_line = p_c_matrix + i*nb_elements;
          double acc = 0;
          for(size_t s = 0; s < nb_elements; s++)
          {
            acc += p_inner_product[s] * current_line[s];
          }
          *it = acc;
        }


        // posts the new value to the listeners
        signal_acc(accumulator);
        signal_counter();
      }




      //! Initialises the accumulator and the signs vector from the first mu
      void initial_accumulation(data_t const &mu)
      {
        accumulator = data_t(data_dimension, 0);
        std::vector<bool>::iterator itb(v_signs.begin());

        // first iteration, we store the signs
        compute_inner_products(mu);
        for(size_t s = 0; s < nb_elements; ++itb, s++)
        {
          bool sign = inner_prod_results[s] >= 0;

          *itb = sign;
          double* current_line = p_c_matrix + s;
          if(sign)
          {
            typename data_t::iterator it(accumulator.begin());
            for(int i = 0; i < data_dimension; i++, current_line += nb_elements, ++it)
            {
              *it += *current_line;              
            }
          }
          else 
          {
            typename data_t::iterator it(accumulator.begin());
            for(int i = 0; i < data_dimension; i++, current_line += nb_elements, ++it)
            {
              *it -= *current_line;              
            }
          }
        }


        // posts the new value to the listeners
        signal_acc(accumulator);
        signal_counter();
      }


      //! Update the accumulator and the signs vector from an upated mu
      void update_accumulation(data_t const& mu)
      {
        std::vector<bool>::iterator itb(v_signs.begin());
        
        compute_inner_products(mu);

        for(size_t s = 0; s < nb_elements; ++itb, s++)
        {
          bool sign = inner_prod_results[s] >= 0;
          if(sign != *itb)
          {
            // update the value of the accumulator according to sign change
            *itb = sign;
            double* current_line = p_c_matrix + s;
            if(sign)
            {
              typename data_t::iterator it(accumulator.begin());
              for(int i = 0; i < data_dimension; i++, current_line += nb_elements, ++it)
              {
                *it += 2* (*current_line);              
              }
            }
            else 
            {
              typename data_t::iterator it(accumulator.begin());
              for(int i = 0; i < data_dimension; i++, current_line += nb_elements, ++it)
              {
                *it -= 2* (*current_line);              
              }
            }
          }
        }

        // posts the new value to the listeners
        signal_acc(accumulator);
        signal_counter();
      }

      //! Project the data onto the orthogonal subspace of the provided vector
      void project_onto_orthogonal_subspace(data_t const &mu)
      {
        // update of vectors in the orthogonal space, and update of the norms at the same time. 
        compute_inner_products(mu);
        double *current_line = p_c_matrix;

        for(int line = 0; line < data_dimension; line ++, current_line += nb_elements)
        {
          double mu_element = mu(line);    
          for(int column = 0; column < nb_elements; column++)
          {
            current_line[column] -= mu_element * inner_prod_results[column];
          }               
        }
  
        signal_counter();
      }

    };


    /*!@internal
     * @brief Accumulation gathering the result of all workers.
     *
     * The purpose of this class is to add the computed accumulator of each thread to the final result
     * which contains the sum of all accumulators. 
     *
     */
    struct asynchronous_results_merger : 
      details::threading::asynchronous_results_merger<
        data_t,
        details::threading::merger_addition<data_t>,
        details::threading::initialisation_vector_specific_dimension<data_t>
        >
    {
    private:
      typedef details::threading::initialisation_vector_specific_dimension<data_t> data_init_type;
      typedef details::threading::merger_addition<data_t> merger_type;
      typedef details::threading::asynchronous_results_merger<data_t, merger_type, data_init_type> parent_type;
    
    public:

      /*!Constructor
       *
       * @param dimension_ the number of dimensions of the vector to accumulate
       */
      asynchronous_results_merger(size_t data_dimension_) : parent_type(data_init_type(data_dimension_))
      {}


    };







  public:

    /*!@brief Constructor
     * 
     * By default, the number of available processors is 1 and the maximum size of the chunks is the maximal size
     */
    robust_pca_impl() : 
      random_init_op(details::fVerySmallButStillComputable, details::fVeryBigButStillComputable), 
      nb_processors(1),
      max_chunk_size(std::numeric_limits<size_t>::max()),
      nb_steps_pca(3)
    {}


    //! Sets the number of parallel tasks used for computing.
    bool set_nb_processors(size_t nb_processors_)
    {
      nb_processors = nb_processors_;
      return true;
    }

    //! Sets the maximum chunk size. 
    //!
    //! By default, the chunk size is the size of the data divided by the number of processing threads.
    //! Lowering the chunk size should provid better granularity in the overall processing time at the end 
    //! of the processing.
    bool set_max_chunk_size(size_t chunk_size)
    {
      if(chunk_size == 0)
      {
        return false;
      }
      max_chunk_size = chunk_size;
      return true;
    }

    //! Sets the number of iterations for the "regular PCA" algorithm. 
    bool set_nb_steps_pca(size_t nb_steps)
    {
      nb_steps_pca = nb_steps;
      return true;
    }



    /*!Performs the computation of the eigen-vectors of the provided dataset.
     *
     * @tparam it_t an input random iterator. Each element pointed by the iterator should be convertible to data_t.
     * @tparam it_eigenvectors an output iterator for storing the computed eigenvalues. This iterator should model a forward output iterator.
     *
     * @param[in] max_iterations the maximum number of iterations in order to compute each eigen-vector. 
     * @param[in] max_dimension_to_compute the maximum number of eigen-vectors to compute.
     * @param[in] it an (input) iterator pointing on the beginning of the data
     * @param[in] ite an (input) iterator pointing on the end of the data
     * @param[out] it_eigenvectors an iterator on the beginning of the area where the computed eigen-vectors will be stored. The space should be at least @c max_dimension_to_compute.
     * @param[in] initial_guess if provided, the initial vectors will be initialized to this value. The size of the pointed container should be at least @c max_dimension_to_compute.
     *
     * @returns true on success, false otherwise
     * @pre 
     * - @c !(it >= ite)
     * - all the vectors given by the iterators pair should be of the same size (no check is performed).
     * - @c std::next(it_eigenvectors, i) should yield a valid iterator pointing on a valid storage area, for @c i in [0, max_dimension_to_compute[.
     *
     */
    template <class it_t, class it_o_eigenvalues_t>
    bool batch_process(
      const size_t max_iterations,
      size_t max_dimension_to_compute,
      it_t const it, 
      it_t const ite, 
      it_o_eigenvalues_t it_eigenvectors,
      std::vector<data_t> const * initial_guess = 0)
    {

      // add some log information
      if(it >= ite)
      {
        return false;
      }

      // preparing the thread pool, to avoid individual thread creation/deletion at each step.
      // we perform the init here because it might take some time for the thread to really start.
      boost::asio::io_service ioService;
      boost::thread_group threadpool;


      // in case of non clean exit (or even in case of clean one).
      details::threading::safe_stop worker_lock_guard(ioService, threadpool);

      // this is exactly the number of processors
      boost::asio::io_service::work work(ioService);
      for(int i = 0; i < nb_processors; i++)
      {
        threadpool.create_thread(boost::bind(&boost::asio::io_service::run, &ioService));
      }

      // contains the number of elements. In case the iterator is random access, could be deduced simply 
      // by a call to distance.
      size_t size_data(std::distance(it, ite));

      // size of the chunks.
      const size_t chunks_size = std::min(max_chunk_size, static_cast<size_t>(size_data/nb_processors));
      const size_t nb_chunks = (size_data + chunks_size - 1) / chunks_size;

      // number of dimensions of the data vectors
      const size_t number_of_dimensions = it->size();
      

      // the first element is used for the init guess because for dynamic std::vector like element, the size is needed.
      data_t mu(initial_guess != 0 ? (*initial_guess)[0] : random_init_op(*it));
      assert(mu.size() == number_of_dimensions);
      typename norm_mu_t::result_type norm_mu(norm_op(mu)); // normalizing
      mu /= norm_mu;

      max_dimension_to_compute = std::min(max_dimension_to_compute, number_of_dimensions);
      
      size_t iterations = 0;


      // preparing the ranges on which each processing thread will run.
      // the number of objects can be much more than the current number of processors, in order to
      // avoid waiting too long for a thread (better granularity) but involving a slight overhead in memory and
      // processing at the synchronization point.
      typedef asynchronous_chunks_processor async_processor_t;
      std::vector<async_processor_t> v_individual_accumulators(nb_chunks);

      asynchronous_results_merger async_merger(number_of_dimensions);

      {
        bool b_result;
        it_t it_current_begin(it);
        for(int i = 0; i < nb_chunks; i++)
        {
          // setting the range
          it_t it_current_end;
          if(i == nb_chunks - 1)
          {
            // just in case the division giving the chunk has some rounding (the parenthesis are important
            // otherwise it is a + followed by a -, which can be out of range after the first +)
            it_current_end = it_current_begin + (size_data - chunks_size*(nb_chunks - 1));
          }
          else
          {
            it_current_end = it_current_begin + chunks_size;
          }

          async_processor_t &current_acc_object = v_individual_accumulators[i];

          // updating the dimension of the problem
          current_acc_object.set_data_dimensions(number_of_dimensions);

          b_result = current_acc_object.set_data_range(it_current_begin, it_current_end);
          if(!b_result)
          {
            return b_result;
          }


          // attaching the update object callbacks
          current_acc_object.connector_accumulator().connect(
            boost::bind(&asynchronous_results_merger::update, &async_merger, _1));
          current_acc_object.connector_counter().connect(
            boost::bind(&asynchronous_results_merger::notify, &async_merger));

          // updating the next 
          it_current_begin = it_current_end;
        }
      }



      // for each dimension
      for(size_t current_dimension = 0; current_dimension < max_dimension_to_compute; current_dimension++, ++it_eigenvectors)
      {


        // PCA like initial steps
        if(nb_steps_pca)
        {
          for(size_t pca_it = 0; pca_it < nb_steps_pca; pca_it++)
          {
            // reseting the final accumulator
            async_merger.init();

            // pushing the initialisation of the mu and sign vectors to the pool
            for(int i = 0; i < v_individual_accumulators.size(); i++)
            {
              ioService.post(
                boost::bind(
                  &async_processor_t::pca_accumulation, 
                  boost::ref(v_individual_accumulators[i]), 
                  boost::cref(mu)));
            }

            // waiting for completion (barrier)
            async_merger.wait_notifications(v_individual_accumulators.size());

            // gathering the first mu
            mu = async_merger.get_merged_result();
            mu /= norm_op(mu);
          }
        }



        data_t previous_mu(mu);
        details::convergence_check<data_t> convergence_op(mu);

        // reseting the final accumulator
        async_merger.init();

        // pushing the initialisation of the mu and sign vectors to the pool
        for(int i = 0; i < v_individual_accumulators.size(); i++)
        {
          ioService.post(
            boost::bind(
              &async_processor_t::initial_accumulation, 
              boost::ref(v_individual_accumulators[i]), 
              boost::cref(previous_mu)));
        }

        // waiting for completion (barrier)
        async_merger.wait_notifications(v_individual_accumulators.size());

        // gathering the first mu
        mu = async_merger.get_merged_result();
        mu /= norm_op(mu);


        // other iterations as usual
        for(iterations = 1; !convergence_op(mu) && iterations < max_iterations; iterations++)
        {
          previous_mu = mu;

          // reseting the final accumulator
          async_merger.init();

          // pushing the update of the mu (and signs)
          for(int i = 0; i < v_individual_accumulators.size(); i++)
          {
            ioService.post(boost::bind(&async_processor_t::update_accumulation, boost::ref(v_individual_accumulators[i]), boost::cref(previous_mu)));
          }

          // waiting for completion (barrier)
          async_merger.wait_notifications(v_individual_accumulators.size());

          // gathering the mus
          mu = async_merger.get_merged_result();
          mu /= norm_op(mu);
        }
        

        // mu is the eigenvector of the current dimension, we store it in the output vector
        *it_eigenvectors = mu;

        // projection onto the orthogonal subspace
        if(current_dimension < max_dimension_to_compute - 1)
        {

          async_merger.init_notifications();

          // pushing the update of the mu (and signs)
          for(int i = 0; i < v_individual_accumulators.size(); i++)
          {
            ioService.post(
              boost::bind(
                &async_processor_t::project_onto_orthogonal_subspace, 
                boost::ref(v_individual_accumulators[i]), 
                boost::cref(*it_eigenvectors))); // this is not mu, since we are changing it before the process ends here
          }

          mu = initial_guess != 0 ? (*initial_guess)[current_dimension+1] : random_init_op(*it);

          async_merger.wait_notifications(v_individual_accumulators.size());

        }
        
      }


      // stopping the pool is done in the destruction of worker_lock_guard



      return true;
    }
  };

}

#endif /* ROBUST_PCA_HPP__ */
