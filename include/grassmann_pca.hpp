// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

#ifndef GRASSMANN_AVERAGES_PCA_HPP__
#define GRASSMANN_AVERAGES_PCA_HPP__

/*!@file
 * Grassmann PCA functions, following the paper of Soren Hauberg.
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


// utilities
#include <include/private/utilities.hpp>
#include <include/private/chunks_processor.hpp>


namespace grassmann_averages_pca
{

  namespace ub = boost::numeric::ublas;

  /*!@brief Grassmann Averages for scalable PCA algorithm
   *
   * This class implements the Grassmann average for computing the PCA in a robust manner. 
   * Its purpose is to compute the PCA of a dataset @f$\{X_i\}@f$, where each @f$X_i@f$ is a vector of dimension
   * D. 
   * 
   * The algorithm is the following:
   * - pick a random or a given @f$\mu_{k, 0}@f$, where @f$k@f$ is the current eigen-vector being computed and @f$0@f$ is the current iteration number (0). 
   * - until the sequence @f$(\mu_{k, t})_t@f$ converges, do:
   *   - computes the sign @f$s_{j, t}@f$ of the projection of the input vectors @f$X_j@f$ onto @f$\mu_{i, t}@f$. We have @f[s_{j, t} = X_j \cdot \mu_{k, t} \geq 0@f]
   *   - compute the update of @f$\mu_{k, .}@f$: @f[\mu_{k, t+1} = \frac{\sum_j s_{j, t} X_j}{\left\|\sum_j s_{j, t} X_j\right\|}@f]
   * - project the @f$X_j@f$'s onto the orthogonal subspace of @f$\mu_{k} = \lim_{t \rightarrow +\infty} \mu_{k, t}@f$: @f[\forall j, X_{j} = X_{j} - X_{j}\cdot\mu_{k} @f]
   *
   * The range taken by @f$k@f$ is a parameter of the algorithm: @c max_dimension_to_compute (see grassmann_pca::batch_process). 
   * The range taken by @f$t@f$ is also a parameter of the algorithm: @c max_iterations (see grassmann_pca::batch_process).
   * The test for convergence is delegated to the class details::convergence_check.
   *
   * The computation is distributed among several threads. The multithreading strategy is 
   * - to split the computation of @f$\sum_j s_{j, t} X_j@f$ among several independant chunks. This computation involves the inner product and the sign. Each chunk addresses 
   *   a subset of the data @f$\{X_j\}@f$ without any overlap with other chunks. The maximal size of a chunk can be configured through the function grassmann_pca::set_max_chunk_size.
   *   By default, the size of the chunk would be the size of the data divided by the number of threads.
   * - to split the computation of the projection onto the orthogonal subspace of @f$\mu_{k}@f$.
   * - to split the computation of the regular PCA algorithm (if any) into several independant chunks.
   *
   * The number of threads can be configured through the function grassmann_pca::set_nb_processors.
   * 
   * @note
   * The algorithm may also perform a few "regular PCA" steps, which is the computation of the eigen-vector with highest eigen-value. This can be configured through the function
   * grassmann_pca::set_nb_steps_pca.
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
  struct grassmann_pca
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

    //! Indicates that the incoming data is not centered and a centering should be performed prior
    //! to the computation of the PCA or the trimmed grassmann average.
    bool need_centering;

    //! An instance observing the steps of the algorithm
    observer_t *observer;    

    //! This is an object that is shared among the three entities of the algorithm:
    //! - the chunk processor
    //! - the main algorithm
    typedef details::mus_distance_optimisation_helper<data_t> inner_product_optimisation_t;


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


      //! Type of the associative map that maps the index of the mu in the algorithm to the delta of vectors
      //! pointing to it. 
      typedef std::map<size_t, int> map_reference_delta_t;
      map_reference_delta_t map_reference_delta;

      /*!Constructor
       *
       * @param dimension_ the number of dimensions of the vector to accumulate
       */
      asynchronous_results_merger(size_t data_dimension_) : parent_type(data_init_type(data_dimension_))
      {}

      //! Override init 
      virtual void init()
      {
        parent_type::init();
        map_reference_delta.clear();
      }

      /*! Receives the deltas of reference count
        * 
        * This is called with the delta of reference count for each of the mus. 
        * @param[in] update_ref_count a map that maps the index of mus with the delta of vectors referencing it
        * @note The call is thread safe.
        */
      void update_mu_reference_counts(std::map<size_t, int> const* update_ref_count)
      {
        boost::unique_lock<typename parent_type::mutex_t> lock(this->internal_mutex);
        for(std::map<size_t, int>::const_iterator it(update_ref_count->begin()), ite(update_ref_count->end());
            it != ite;
            ++it)
        {
          map_reference_delta[it->first] += it->second;
        }
      }


      //! Returns the current merged results.
      //! @warning the call is not thread safe (intended to be called once the wait_notifications returned and no
      //! other thread is working). 
      map_reference_delta_t const& get_merged_reference_counts() const
      {
        return map_reference_delta;
      }

    };



    //! Transfers the updates of the mu counts to the dedicated class instance. 
    void update_mus_reference(asynchronous_results_merger const& async_merger, inner_product_optimisation_t& inner_product_optimisation)
    {
      typename asynchronous_results_merger::map_reference_delta_t const& map_ref_updates = async_merger.get_merged_reference_counts();
      for(typename asynchronous_results_merger::map_reference_delta_t::const_iterator it_reference_mu(map_ref_updates.begin()), ite_reference_mu(map_ref_updates.end());
          it_reference_mu != ite_reference_mu;
          ++it_reference_mu)
      {
        inner_product_optimisation.update_count(it_reference_mu->first, it_reference_mu->second);
      }

    }





  public:

    /*!@brief Constructor
     * 
     * @note By default the number of processors used for computation is set to 1.
     * The maximum size of the chunks is "infinite": each chunk will receive in that case the size of the data
     * divided by the number of running threads.
     */
    grassmann_pca() : 
      random_init_op(details::fVerySmallButStillComputable, details::fVeryBigButStillComputable), 
      nb_processors(1),
      max_chunk_size(std::numeric_limits<size_t>::max()),
      nb_steps_pca(3),
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

    //! Sets the number of iterations for the "regular PCA" algorithm. 
    bool set_nb_steps_pca(size_t nb_steps)
    {
      nb_steps_pca = nb_steps;
      return true;
    }

    //! Sets the centering flags.
    //!
    //! If set to true, a centering will be performed before applying any computation (PCA and Grassmann averages). 
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
      

      // the first element is used for the init guess because for dynamic std::vector like element, the size is needed.
      data_t mu(initial_guess != 0 ? (*initial_guess)[0] : random_init_op(*it));
      mu *= typename data_t::value_type(1./norm_op(mu)); // normalizing
      assert(mu.size() == number_of_dimensions);

      max_dimension_to_compute = std::min(max_dimension_to_compute, number_of_dimensions);
      

      // this object would allow to minimise the number of inner products computed inside the chunk processors
      // by keeping track of the angles. The chunk processors need in turn keep track of some indices and angles
      // the refer to the mus 
      inner_product_optimisation_t inner_product_optimisation;


      size_t iterations = 0;
      
      // preparing the ranges on which each processing thread will run.
      // the number of objects can be much more than the current number of processors, in order to
      // avoid waiting too long for a thread (better granularity) but involving a slight overhead in memory and
      // processing at the synchronization point.
      typedef asynchronous_chunks_processor<data_t> async_processor_t;
      std::vector<async_processor_t> v_individual_accumulators(nb_chunks, async_processor_t(inner_product_optimisation));

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
          current_acc_object.connector_reference_count() = boost::bind(&asynchronous_results_merger::update_mu_reference_counts, &async_merger, _1);

          // updating the dimension of the problem
          current_acc_object.set_data_dimensions(number_of_dimensions);

          // pushing the asynchronous copy
          ioService.post(
            boost::bind(
              &async_processor_t::template set_data_range<it_t>, 
              boost::ref(v_individual_accumulators[i]), 
              it_current_begin, it_current_end));


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
            
            double norm_mu = norm_op(mu);
            if(norm_mu < 1E-12)
            {
              if(observer)
              {
                std::ostringstream o;
                o << "The result of the PCA is null for subspace " << current_subspace_index;
                observer->log_error_message(o.str().c_str());
              }
              return false;
            }            
            
            mu *= typename data_t::value_type(1./norm_mu);
          }
          
          // sending result to observer
          if(observer)
          {
            observer->signal_pca(mu, current_subspace_index);
          }          
        }


        // init the convergence object
        details::convergence_check<data_t> convergence_op(mu);

        // register the initial mu
        inner_product_optimisation.add_mu(mu, 0, 0);


        // reseting the accumulator and the notifications
        async_merger.init();

        // pushing the initialisation of the mu and sign vectors to the pool
        for(int i = 0; i < v_individual_accumulators.size(); i++)
        {
          ioService.post(
            boost::bind(
              &async_processor_t::initial_accumulation, 
              boost::ref(v_individual_accumulators[i]), 
              boost::cref(mu)));
        }

        
        // waiting for completion (barrier)
        async_merger.wait_notifications(v_individual_accumulators.size());

        // gathering the first mu
        mu = async_merger.get_merged_result();
        mu *= typename data_t::value_type(1./norm_op(mu));


        // updating the references, for this case, the initial mu should be referenced by all
        update_mus_reference(async_merger, inner_product_optimisation);


        // other iterations as usual
        for(iterations = 1; !convergence_op(mu) && iterations < max_iterations; iterations++)
        {
          // new mu is added to the optimised, and anything not needed is removed
          inner_product_optimisation.add_mu(mu, iterations, 0);



          // reseting the notifications (just the notifications, no need for init)
          async_merger.init_notifications();

          // pushing the update of the mu (and signs)
          for(int i = 0; i < v_individual_accumulators.size(); i++)
          {
            ioService.post(
              boost::bind(
                &async_processor_t::update_accumulation, 
                boost::ref(v_individual_accumulators[i]), // instance
                boost::cref(mu), iterations)); // params
          }

          // waiting for completion (barrier)
          async_merger.wait_notifications(v_individual_accumulators.size());

          // gathering the mus
          mu = async_merger.get_merged_result();
          mu *= typename data_t::value_type(1./norm_op(mu));

          // updating the references
          update_mus_reference(async_merger, inner_product_optimisation);
          inner_product_optimisation.prune();


          // sending result to observer
          if(observer)
          {
            observer->signal_intermediate_result(mu, current_subspace_index, iterations);
          }
        }

        // clear the optimisation part, not needed anymore
        inner_product_optimisation.clear();


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
                &async_processor_t::template project_onto_orthogonal_subspace<typename it_o_basisvectors_t::value_type>,
                boost::ref(v_individual_accumulators[i]), 
                *it_basisvectors)); // this is not mu, since we are changing it before the process ends here
          }

          mu = initial_guess != 0 ? (*initial_guess)[current_subspace_index+1] : random_init_op(*it);

          async_merger.wait_notifications(v_individual_accumulators.size());

        }
        
      }


      // stopping the pool is done in the destruction of worker_lock_guard



      return true;
    }
  };

}

#endif /* GRASSMANN_AVERAGES_PCA_HPP__ */
