// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

#ifndef SIMPLE_AVERAGES_PCA_HPP__
#define SIMPLE_AVERAGES_PCA_HPP__

/*!@file
 * Simple PCA functions, for comparison with the Grassmann averages.
 *
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


namespace grassmann_averages_pca
{

  namespace ub = boost::numeric::ublas;
   
  /*!@brief PCA algorithm.
   *
   * This class is intended to be used in comparison to the Grassmann algorithms.
   * Its purpose is to compute the PCA of a dataset @f$\{X_i\}@f$, where each @f$X_i@f$ is a vector of dimension
   * D. 
   * 
   * The algorithm is the following:
   * - pick a random or a given @f$\mu_{k, 0}@f$, where @f$k@f$ is the current eigen-vector being computed and @f$0@f$ is the current iteration number (0). 
   * - until the sequence @f$(\mu_{k, t})_t@f$ converges, do:
   *   - computes the projection of the input vectors @f$X_j@f$ onto @f$\mu_{i, t}@f$
   *   - accumulate these projections and update @f$\mu_{k, t}@f$: @f[\mu_{k, t+1} = \frac{\sum_j \left(\mu_{i, t}^T \cdot X_j \right) X_j}{\left\|\sum_j \left(\mu_{i, t}^T \cdot X_j \right) X_j\right\|}@f]
   * - project the @f$X_j@f$'s onto the orthogonal subspace of @f$\mu_{k} = \lim_{t \rightarrow +\infty} \mu_{k, t}@f$: @f[\forall j, X_{j} = X_{j} - X_{j}\cdot\mu_{k} @f]
   *
   * The range taken by @f$k@f$ is a parameter of the algorithm: @c max_dimension_to_compute (see em_pca::batch_process). 
   * The range taken by @f$t@f$ is also a parameter of the algorithm: @c max_iterations (see em_pca::batch_process).
   * The test for convergence is delegated to the class details::convergence_check.
   *
   * The computation is distributed among several threads. The multithreading strategy is 
   * - to split the computation of @f$\left( \mu_{i, t}^T \cdot X_j \right) X_j@f$ among several independant chunks. This computation involves the inner product and the addition. Each chunk addresses 
   *   a subset of the data @f$\{X_j\}@f$ without any overlap with other chunks. The maximal size of a chunk can be configured through the function em_pca::set_max_chunk_size.
   *   By default, the size of the chunk would be the size of the data divided by the number of threads.
   * - to split the computation of the projection onto the orthogonal subspace of @f$\mu_{k}@f$.
   *
   * The number of threads can be configured through the function em_pca::set_nb_processors.
   * 
   * @note
   * The API is chosen to be the same as for the other two Grassmann implementations, in order to be able to easily switch from one implementation to the other.
   *
   * @tparam data_t type of vectors used for the computation. 
   * @tparam observer_t an observer type following the signature of the class grassmann_trivial_callback.
   * @tparam norm_mu_t norm used to normalize the eigen-vector and project them onto the unit circle.
   *
   * @author Soren Hauberg, Raffi Enficiaud
   */
  template <class data_t, 
            class observer_t = grassmann_trivial_callback<data_t>,
            class norm_mu_t = details::norm2>
  struct em_pca
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

    //! Indicates that the incoming data is not centered and a centering should be performed prior
    //! to the computation of the PCA.
    bool need_centering;

    //! An instance observing the steps of the algorithm
    observer_t *observer;    

    //!@internal
    //!@brief Contains the logic for processing part of the accumulator
    struct asynchronous_chunks_processor
    {
    private:
      //! Type of the elements contained in the vectors
      typedef typename data_t::value_type scalar_t;

      //! Number of vectors contained in this chunk.
      size_t nb_elements;

      //! Dimension of the vectors.
      size_t data_dimension;

      //! Internal accumulator.
      //! Should live beyond the scope of update and init, as required by the merger.
      data_t accumulator;

      // this is to send an update of the value of mu to one listener
      // the connexion should be managed externally
      typedef boost::function<void (data_t const*)> connector_accumulator_t;
      connector_accumulator_t signal_acc;

      typedef boost::function<void ()> connector_counter_t;
      connector_counter_t signal_counter;

      //! The matrix containing a copy of the data.
      //! The vectors are stored per row in this matrix.
      scalar_t *p_c_matrix;
      
      //! Padding for one line of the matrix
      size_t data_padding;

      //! "Optimized" inner product.
      //! This one has the particularity to be more cache/memory bandwidth friendly. More efficient
      //! implementations may be used, but it turned out that the memory bandwidth is saturated when
      //! many threads are processing different data.
      scalar_t inner_product(scalar_t const* p_mu, scalar_t const* current_line) const
      {
        scalar_t const * const current_line_end = current_line + data_dimension;

        const int _64_elements = static_cast<int>(data_dimension >> 6);
        scalar_t acc(0);

        for(int j = 0; j < _64_elements; j++, current_line += 64, p_mu += 64)
        {
          for(int i = 0; i < 64; i++)
          {
            acc += current_line[i] * p_mu[i];
          }
        }
        for(; current_line < current_line_end; current_line++, p_mu++)
        {
          acc += (*current_line) * (*p_mu);
        }
        return acc;
      }

      //! "Optimized" inner product
      scalar_t inner_product(scalar_t const* p_mu, size_t element_index) const
      {
        return inner_product(p_mu, p_c_matrix + element_index * data_padding);
      }



    public:
      asynchronous_chunks_processor() : nb_elements(0), data_dimension(0), p_c_matrix(0), data_padding(0)
      {
      }
      
      ~asynchronous_chunks_processor()
      {
        delete [] p_c_matrix;
      }


      //! Sets the data range
      template <class container_iterator_t>
      void set_data_range(container_iterator_t const &b, container_iterator_t const& e)
      {
        nb_elements = std::distance(b, e);
        assert(nb_elements > 0);

        // aligning on 32 bytes = 1 << 5
        data_padding = (data_dimension*sizeof(scalar_t) + (1<<5) - 1) & (~((1<<5)-1));
        data_padding /= sizeof(scalar_t);
        
        delete [] p_c_matrix;
        p_c_matrix = new scalar_t[data_padding*nb_elements];
        
        container_iterator_t bb(b);

        scalar_t *current_line = p_c_matrix;
        for(int line = 0; line < nb_elements; line ++, current_line += data_padding, ++bb)
        {         
          for(int column = 0; column < data_dimension; column++)
          {
            current_line[column] = (*bb)(column);
          }
        }
        
        signal_counter();
      }

      //! Sets the dimension of each vectors
      //! @pre data_dimensions_ strictly positive
      void set_data_dimensions(size_t data_dimensions_)
      {
        data_dimension = data_dimensions_;
        assert(data_dimension > 0);
      }

      //! Returns the callback object that will be called to signal an updated accumulator.
      connector_accumulator_t& connector_accumulator()
      {
        return signal_acc;
      }

      //! Returns the callback object that will be called to signal the end of the current computation.
      connector_counter_t& connector_counter()
      {
        return signal_counter;
      }


      //! Centering the data in case it was not possible to do it beforehand
      void data_centering_first_phase(size_t full_dataset_size)
      {
        scalar_t const * current_line = p_c_matrix;
        accumulator = data_t(data_dimension, 0);
        scalar_t * const p_acc_begin = &accumulator.data()[0];
        scalar_t const * const p_acc_end = p_acc_begin + data_dimension;
        
        
        for(size_t current_element = 0; 
            current_element < nb_elements; 
            current_element++, current_line+= data_padding - data_dimension)
        {
          scalar_t * p_acc = p_acc_begin;
          for(; p_acc < p_acc_end; p_acc++, current_line++)
          {
            *p_acc += *current_line;
          }
        }


        for(scalar_t * p_acc = p_acc_begin; p_acc < p_acc_end; p_acc++)
        {
          *p_acc /= full_dataset_size;
        }


        // posts the new value to the listeners for the current dimension
        signal_acc(&accumulator);
        signal_counter();

      }

      //! Project the data onto the orthogonal subspace of the provided vector
      void data_centering_second_phase(data_t const &mean_value)
      {
        scalar_t * current_line = p_c_matrix;
        scalar_t const * const p_mean_begin = &mean_value.data()[0];
        scalar_t const * const p_mean_end = p_mean_begin + data_dimension;
        
        
        for(size_t current_element = 0; 
            current_element < nb_elements; 
            current_element++, current_line+= data_padding - data_dimension)
        {
          const scalar_t * p_mean = p_mean_begin;
          for(; p_mean < p_mean_end; p_mean++, current_line++)
          {
            *current_line -= *p_mean;
          }
        }

        // posts the new value to the listeners for the current dimension
        signal_counter();

      }


      //! PCA steps
      void pca_accumulation(data_t const &mu)
      {
        accumulator = data_t(data_dimension, 0);
        scalar_t const * const p_mu = &mu.data()[0];
        scalar_t * const p_acc = &accumulator.data()[0];
                
        scalar_t const * current_line = p_c_matrix;

        for(size_t s = 0; s < nb_elements; s++, current_line += data_padding)
        {
          const scalar_t inner_prod = inner_product(p_mu, current_line);
          for(size_t d = 0; d < data_dimension; d++)
          {
            p_acc[d] += inner_prod * current_line[d];
          }
        }

        // posts the new value to the listeners
        signal_acc(&accumulator);
        signal_counter();
      }

      //! Project the data onto the orthogonal subspace of the provided vector
      void project_onto_orthogonal_subspace(data_t const &mu)
      {
        // update of vectors in the orthogonal space, and update of the norms at the same time. 
        
        scalar_t const * const p_mu = &mu.data()[0];
        scalar_t * current_line = p_c_matrix;
        
        for(size_t line = 0; line < nb_elements; line ++, current_line += data_padding)
        {
          scalar_t const inner_prod = inner_product(p_mu, current_line);
          for(int column = 0; column < data_dimension; column++)
          {
            current_line[column] -= inner_prod * p_mu[column];
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
     * @note By default the number of processors used for computation is set to 1.
     * The maximum size of the chunks is "infinite": each chunk will receive in that case the size of the data
     * divided by the number of running threads.
     */
    em_pca() : 
      random_init_op(details::fVerySmallButStillComputable, details::fVeryBigButStillComputable), 
      nb_processors(1),
      max_chunk_size(std::numeric_limits<size_t>::max()),
      need_centering(false),
      observer(0)
    {}

    //! Sets the observer of the algorithm. 
    //!
    //! The lifetime of the observer is not managed by this class. Set to 0 to disable
    //! observation.
    bool set_observer(observer_t* observer_)
    {
      observer = observer_;
      return true;
    }

    //! Sets the number of parallel tasks used for computing.
    bool set_nb_processors(size_t nb_processors_)
    {
      nb_processors = nb_processors_;
      return true;
    }

    /*!@brief Sets the maximum chunk size. 
     *
     * By default, the chunk size is the size of the data divided by the number of processing threads.
     * Lowering the chunk size should provid better granularity in the overall processing time at the end 
     * of the processing.
     */
    bool set_max_chunk_size(size_t chunk_size)
    {
      if(chunk_size == 0)
      {
        return false;
      }
      max_chunk_size = chunk_size;
      return true;
    }

    //! Sets the centering flags.
    //!
    //! If set to true, a centering will be performed before applying any computation. 
    bool set_centering(bool need_centering_)
    {
      need_centering = need_centering_;
      return true;
    }



    /*!@brief Performs the computation of the eigen-vectors of the provided dataset.
     *
     * @tparam it_t an input random iterator. Each element pointed by the iterator should be convertible to data_t.
     * @tparam it_o_basisvectors_t an output iterator for storing the computed eigenvalues. This iterator should model a forward output iterator.
     *
     * @param[in] max_iterations the maximum number of iterations in order to compute each eigen-vector. 
     * @param[in] max_dimension_to_compute the maximum number of eigen-vectors to compute.
     * @param[in] it an (input) iterator pointing on the beginning of the data
     * @param[in] ite an (input) iterator pointing on the end of the data
     * @param[out] it_basisvectors an iterator on the beginning of the area where the computed eigen-vectors will be stored. The space should be at least @c max_dimension_to_compute.
     * @param[in] initial_guess if provided, the initial vectors will be initialized to this value. The size of the pointed container should be at least @c max_dimension_to_compute.
     *
     * @returns true on success, false otherwise
     * @pre 
     * - @c !(it >= ite)
     * - all the vectors given by the iterators pair should be of the same size (no check is performed).
     * - @c std::next(it_eigenvectors, i) should yield a valid iterator pointing on a valid storage area, for @c i in [0, max_dimension_to_compute[.
     *
     */
    template <class it_t, class it_o_basisvectors_t>
    bool batch_process(
      const size_t max_iterations,
      size_t max_dimension_to_compute,
      it_t const it, 
      it_t const ite, 
      it_o_basisvectors_t it_basisvectors,
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
      const size_t chunks_size = std::min(max_chunk_size, static_cast<size_t>(ceil(double(size_data)/nb_processors)));
      const size_t nb_chunks = (size_data + chunks_size - 1) / chunks_size;

      // number of dimensions of the data vectors
      const size_t number_of_dimensions = it->size();
      
      // pointer to the begining of the output vectors, for orthonormalisation
      it_o_basisvectors_t const it_basisvectors_begin(it_basisvectors);

      // the first element is used for the init guess because for dynamic std::vector like element, the size is needed.
      data_t mu(initial_guess != 0 ? (*initial_guess)[0] : random_init_op(*it));
      mu *= typename data_t::value_type(1./norm_op(mu)); // normalizing
      assert(mu.size() == number_of_dimensions);

      max_dimension_to_compute = std::min(max_dimension_to_compute, number_of_dimensions);
      
      size_t iterations = 0;


      // preparing the ranges on which each processing thread will run.
      // the number of objects can be much more than the current number of processors, in order to
      // avoid waiting too long for a thread (better granularity) but involving a slight overhead in memory and
      // processing at the synchronization point.
      typedef asynchronous_chunks_processor async_processor_t;
      std::vector<async_processor_t> v_individual_accumulators(nb_chunks);

      asynchronous_results_merger async_merger(number_of_dimensions);
      async_merger.init_notifications();

      {
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

          // attaching the update object callbacks
          current_acc_object.connector_accumulator() = boost::bind(&asynchronous_results_merger::update, &async_merger, _1);
          current_acc_object.connector_counter() = boost::bind(&asynchronous_results_merger::notify, &async_merger);

          // updating the dimension of the problem
          current_acc_object.set_data_dimensions(number_of_dimensions);

          // pushing the asynchronous copy
          ioService.post(
            boost::bind(
              &async_processor_t::template set_data_range<it_t>, 
              boost::ref(v_individual_accumulators[i]), 
              it_current_begin, it_current_end));

          //bool b_result = current_acc_object.set_data_range(it_current_begin, it_current_end);
          //if(!b_result)
          //{
          //  return b_result;
          //}


          // updating the next 
          it_current_begin = it_current_end;
        }
        
        // waiting for completion (barrier)
        async_merger.wait_notifications(v_individual_accumulators.size());
        
      }


      // Centering the data if needed: 
      // - first run the accumulation and gather all results in a multithreaded manner
      // - second center the data with the collected mean
      if(need_centering)
      {
        // Computing the accumulation
        async_merger.init();

        for(int i = 0; i < v_individual_accumulators.size(); i++)
        {
          ioService.post(
            boost::bind(
              &async_processor_t::data_centering_first_phase, 
              boost::ref(v_individual_accumulators[i]),
              size_data)); // size of the dataset to perform division and avoid doing accumulation over big numerical values
        }

        // waiting for completion (barrier)
        async_merger.wait_notifications(v_individual_accumulators.size());

        // gathering the accumulated, already divided by the size 
        data_t mean_vector = async_merger.get_merged_result();

        // sending result to observer
        if(observer)
        {
          observer->signal_mean(mean_vector);
        }


        // centering the data
        async_merger.init();

        for(int i = 0; i < v_individual_accumulators.size(); i++)
        {
          ioService.post(
            boost::bind(
              &async_processor_t::data_centering_second_phase, 
              boost::ref(v_individual_accumulators[i]),
              boost::cref(mean_vector)
              ));
        }

        // waiting for completion (barrier)
        async_merger.wait_notifications(v_individual_accumulators.size());


      }


      // for each dimension
      for(size_t current_subspace_index = 0; 
          current_subspace_index < max_dimension_to_compute; 
          current_subspace_index++, ++it_basisvectors)
      {

        details::convergence_check<data_t> convergence_op(mu);

        // other iterations as usual
        for(iterations = 0; (!convergence_op(mu) && (iterations < max_iterations)) || (iterations == 0); iterations++)
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
            
          double norm_mu = norm_op(mu);
          if(norm_mu < 1E-12)
          {
            if(observer)
            {
              std::ostringstream o;
              o << "The result of the PCA is null for subspace " 
                << current_subspace_index
                << " at iteration @ "
                << iterations;
              observer->log_error_message(o.str().c_str());
            }
            return false;
          }            
            
          mu *= typename data_t::value_type(1./norm_mu);


          // sending result to observer
          if(observer)
          {
            observer->signal_intermediate_result(mu, current_subspace_index, iterations);
          }
        }
          
        // mu is the eigenvector of the current dimension, we store it in the output vector
        *it_basisvectors = mu;

        // sending result to observer
        if(observer)
        {
          observer->signal_eigenvector(*it_basisvectors, current_subspace_index);
        }   

        // projection onto the orthogonal subspace
        if(current_subspace_index < max_dimension_to_compute - 1)
        {

          async_merger.init_notifications();

          // pushing the update of the mu (and signs)
          for(int i = 0; i < v_individual_accumulators.size(); i++)
          {
            ioService.post(
              boost::bind(
                &async_processor_t::project_onto_orthogonal_subspace, 
                boost::ref(v_individual_accumulators[i]), 
                boost::cref(*it_basisvectors))); // this is not mu, since we are changing it before the process ends here
          }

          mu = initial_guess != 0 ? (*initial_guess)[current_subspace_index+1] : random_init_op(*it);

          // project onto the orthogonal subspace
          for(it_o_basisvectors_t it_orthonormalised_element(it_basisvectors_begin); 
              it_orthonormalised_element <= it_basisvectors; 
              ++it_orthonormalised_element)
          {
            mu -= boost::numeric::ublas::inner_prod(mu, *it_orthonormalised_element) * (*it_orthonormalised_element);
          }
          mu *= typename it_t::value_type::value_type(1./norm_op(mu));

          async_merger.wait_notifications(v_individual_accumulators.size());

        }
        
      }


      // stopping the pool is done in the destruction of worker_lock_guard



      return true;
    }
  };

}

#endif /* SIMPLE_AVERAGES_PCA_HPP__ */
