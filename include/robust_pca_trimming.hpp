// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

#ifndef GRASSMANN_AVERAGES_PCA_TRIMMED_HPP__
#define GRASSMANN_AVERAGES_PCA_TRIMMED_HPP__



/*!@file
 * Robust PCA functions, following the paper of Soren Hauberg.
 *
 * This file contains the implementation of the trimmed version. 
 */


// for the thread pools
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <numeric>

#include <boost/scoped_array.hpp>
#include <boost/function.hpp>

// utilities
#include <include/private/utilities.hpp>


namespace grassmann_averages_pca
{


  namespace details
  {

    //!@internal
    //!Helper object for the updates
    template <class scalar_t>
    struct s_dimension_update
    {
      size_t dimension;
      scalar_t value;
    };

    /*!@internal
     * @brief Adaptation of merger_addition concept for accumulator and count at the same time.
     * 
     * The trimmed version of the robust pca algorithm may strip some element along each dimension. 
     * In order to compute the @f$\mu@f$ properly, the count should also be transfered.
     */
    template <class data_t>
    struct merger_update_specific_dimension
    {
      typedef data_t input_t;
      bool operator()(input_t &current_state, s_dimension_update<typename data_t::value_type> const& update_value) const
      {
        current_state(update_value.dimension) += update_value.value;
        return true;
      }
    };
    
  }


  /*!@brief Grassmann Average algorithm for robust PCA computation, with trimming of outliers.
   *
   * This class implements the robust PCA using the Grassmann average with trimming of the outliers. 
   * Its purpose is to compute the PCA of a dataset @f$\mathbf{X} = \{X_i\}@f$, where each @f$X_i@f$ is a vector of dimension
   * D, and also by being "more" robust to outliers that might occur in the original data.
   * The particularity of the Grassman average scheme is to be more stable than other algorithms against the dimension of the data.
   * 
   * The algorithm is the following:
   * - pick a random or a given @f$\mu_{k, 0}@f$, where @f$k@f$ is the current basis vector being computed and @f$0@f$ is the current iteration number (0). 
   * - ensure this @f$\mu_{k, 0}@f$ is orthogonal to the previous detected @f$\mu_{k', 0},\, \forall k' \in [0, k)@f$
   * - until the sequence @f$(\mu_{k, t})_t@f$ converges, do:
   *   - computes the sign @f$s_{j, t}@f$ of the projection of the input vectors @f$X_j@f$ onto @f$\mu_{i, t}@f$. We have @f[s_{j, t} = X_j \cdot \mu_{k, t} \geq 0@f]
   *   - for each dimension @f$0 \lt d \leq D@f$, do:
   *     - For all @f$j@f$, consider the multiplied data set projected onto dimension @f$d@f$: 
   *       @f[\mathbf{X_\mu^{(d)}} = \left\{proj_d \left(s_{j, t} \cdot X_j\right) \right\} = \left\{s_{j, t} \cdot X_j^{(d)}\right\}@f]
   *       which is a 1-D sequence
   *     - compute the indexes @f$J_d@f$ of the @f$\frac{K}{2}@f$ lowest and biggest points of this 1-D sequence @f$\mathbf{X_\mu^{(d)}}@f$
   *     - compute the update of @f$\mu_{k, .}@f$ for dimension @f$d@f$ : 
   *       @f[\mu_{k, t+1}^{(d)} = \frac{\sum_{j \notin J_d} proj_d \left(s_{j, t} \cdot X_j \right)}{\# J - \# J_d} @f]
   *   - normalize @f[\mu_{k, t+1} = \frac{\left(\mu_{k, t+1}^{(1)}, \ldots \mu_{k, t+1}^{(D)}\right)^t}{\left\|\left(\mu_{k, t+1}^{(1)}, \ldots \mu_{k, t+1}^{(D)}\right)^t\right\|}@f]
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
   * - to compute the @f$\frac{K}{2}@f$ extremal points along each dimension in parallel
   * - to split the computation of the projection onto the orthogonal subspace of @f$\mu_{k}@f$.
   * - to split the computation of the regular PCA algorithm (if any) into several independant chunks.
   *
   * The number of threads can be configured through the function grassmann_pca::set_nb_processors.
   * 
   * @note
   * The algorithm may also perform a few "regular PCA" steps, which is the computation of the basis vector with highest "eigen-value". This can be configured through the function
   * grassmann_pca::set_nb_steps_pca.
   *
   * @tparam data_t type of vectors used for the computation. 
   * @tparam norm_mu_t norm used to normalize the basis vectors and project them onto the unit circle.
   *
   * @author Soren Hauberg, Raffi Enficiaud
   */
  template <class data_t, class norm_mu_t = details::norm2>
  struct grassmann_pca_with_trimming
  {
  private:
    //! Random generator for initialising @f$\mu@f$ at each dimension. 
    details::random_data_generator<data_t> random_init_op;

    //! Norm used for normalizing @f$\mu@f$.
    norm_mu_t norm_op;
    

    //! The percentage of the data that should be trimmed.
    //! The trimming is performed symmetrically in the upper and lower distribution of the data, hence
    //! each side is trimmed by trimming_percentage/2.
    double trimming_percentage;

    //! Number of parallel tasks that will be used for computing.
    size_t nb_processors;

    //! Maximal size of a chunk (infinity by default).
    size_t max_chunk_size;

    //! Number of steps for the initial PCA like algorithm (defaults to 3).
    size_t nb_steps_pca;

    //! Type of the element of data_t. 
    typedef typename data_t::value_type scalar_t;

    //! Type of the vector used for counting the element falling into the non-trimmed range.
    typedef boost::numeric::ublas::vector<size_t> count_vector_t;



    //!@internal
    //!@brief Contains the logic for processing part of the accumulator
    struct s_robust_pca_trimmed_processor_inner_products
    {
    private:
      //! Scalar type of the data
      typedef typename data_t::value_type scalar_t;
      
      typedef details::s_dimension_update<scalar_t> accumulator_element_t;

      
      size_t nb_elements;                 //!< The size of the current dataset
      size_t data_dimension;              //!< The dimension of the data
      size_t nb_elements_to_keep;         //!< The number of elements to keep.
      
      typedef boost::function<void ()> connector_counter_t;
      connector_counter_t signal_counter;

      typedef boost::function<void (accumulator_element_t const*)> connector_accumulator_dimension_t;
      connector_accumulator_dimension_t signal_acc_dimension;

      std::vector<accumulator_element_t> v_accumulated_per_dimension;

 
      //! The matrix containing a copy of the data
      scalar_t *p_c_matrix;
      
      //! The result of the inner products
      std::vector<scalar_t> inner_prod_results;
      std::vector<scalar_t> current_dimension_sort;
  
      void compute_inner_products(data_t const &mu)
      {
        scalar_t *out = &inner_prod_results[0];
        scalar_t mu_element = mu(0);
        scalar_t *current_line = p_c_matrix;
        
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


    public:
      s_robust_pca_trimmed_processor_inner_products() : 
        nb_elements(0), 
        data_dimension(0), 
        nb_elements_to_keep(0),
        p_c_matrix(0)
      {
      }

      ~s_robust_pca_trimmed_processor_inner_products()
      {
        delete [] p_c_matrix;
      }
      
      
      //! Returns the connected object that will receive the notification of the end of the current process.
      connector_counter_t& connector_counter()
      {
        return signal_counter;
      }     
      
      connector_accumulator_dimension_t& connector_accumulator()
      {
        return signal_acc_dimension;
      }
      
      //! Sets the data range
      template <class container_iterator_t>
      bool set_data_range(container_iterator_t const &b, container_iterator_t const& e)
      {
        if(data_dimension <= 0)
        {
          return false;
        }
        
        nb_elements = std::distance(b, e);

        assert(data_dimension > 0);
        assert(nb_elements > 0);
        
        delete [] p_c_matrix;
        p_c_matrix = new scalar_t[nb_elements*data_dimension];
        
        container_iterator_t bb(b);
        
        for(int column = 0; column < nb_elements; column++, ++bb)
        {
          scalar_t* current_line = p_c_matrix + column;
          for(int line = 0; line < data_dimension; line ++, current_line += nb_elements)
          {         
            *current_line = (*bb)(line);
          }
          
        }
        
        inner_prod_results.resize(nb_elements);
        current_dimension_sort.reserve(nb_elements);
        return true;
      }

      //! Sets the dimension of the data vectors
      //! @pre data_dimensions_ is strictly positive
      void set_data_dimensions(size_t data_dimensions_)
      {
        assert(data_dimensions_ > 0);
        data_dimension = data_dimensions_;
        v_accumulated_per_dimension.resize(data_dimension);
      }

      //! Sets the number of element to keep in the upper and lower distributions.
      void set_nb_elements_to_keep(int nb_elements_to_keep_)
      {
        assert(nb_elements_to_keep_ >= 0);
        nb_elements_to_keep = nb_elements_to_keep_;
      }



      //! PCA steps
      void pca_accumulation(data_t const &mu)
      {
        compute_inner_products(mu);
               
        scalar_t const * const p_inner_product = &inner_prod_results[0];

        for(size_t dimension = 0; dimension < data_dimension; dimension++)
        {
          scalar_t const * const current_line = p_c_matrix + dimension*nb_elements;
          scalar_t acc = 0;
          for(size_t s = 0; s < nb_elements; s++)
          {
            acc += p_inner_product[s] * current_line[s];
          }

          // posts the new value to the listeners for the current dimension
          accumulator_element_t &result = v_accumulated_per_dimension[dimension];
          result.dimension = dimension;
          result.value = acc;
          signal_acc_dimension(&result);
        }

        signal_counter();
      }

      
      void compute_data_matrix(data_t const &mu, scalar_t* p_out, size_t padding)
      {
        // updates the internal inner products
        compute_inner_products(mu);
        
        // the current line is spans a particular dimension
        scalar_t *current_line = p_c_matrix;

        // this spans the inner product results for all dimensions

        std::vector<int> v_mult(nb_elements);
        for(size_t element(0); element < nb_elements; element++)
        {
          v_mult[element] = inner_prod_results[element] >= 0 ? 1 : -1;
        }
        int const * const out = &v_mult[0];

        for(size_t current_dimension = 0; 
            current_dimension < data_dimension; 
            current_dimension++, current_line += nb_elements, p_out+= padding)
        {
          for(size_t element(0); element < nb_elements; element++)
          {
            p_out[element] = out[element] * current_line[element];
          }
        }
        
        // signals the main merger
        signal_counter();
      }
      
      
      void compute_bounded_accumulation(size_t dimension, size_t nb_total_elements, scalar_t* p_data)
      {
        std::nth_element(p_data, p_data + nb_elements_to_keep, p_data + nb_total_elements);
        std::nth_element(p_data + nb_elements_to_keep+1, p_data + nb_total_elements - nb_elements_to_keep-1, p_data + nb_total_elements);
        scalar_t acc = std::accumulate(p_data + nb_elements_to_keep, p_data + nb_total_elements - nb_elements_to_keep, scalar_t(0.));
        
        accumulator_element_t &result = v_accumulated_per_dimension[dimension];
        result.dimension = dimension;

        result.value = acc / (nb_total_elements - 2*nb_elements_to_keep);
        // signals the update
        signal_acc_dimension(&result);
        // signals the main merger
        signal_counter();
      }
      
      //! Project the data onto the orthogonal subspace of the provided vector
      void project_onto_orthogonal_subspace(data_t const &mu)
      {
        compute_inner_products(mu);
        scalar_t *current_line = p_c_matrix;
        
        for(int line = 0; line < data_dimension; line ++, current_line += nb_elements)
        {
          scalar_t mu_element = mu(line);    
          for(int column = 0; column < nb_elements; column++)
          {
            current_line[column] -= mu_element * inner_prod_results[column];
          }               
        }

        signal_counter();
      }

    };


    /*!@internal
     * @brief Merges the result of all workers and signals the results to the main thread.
     *
     * The purpose of this class is to add the computed accumulator of each thread to the final result
     * which contains the sum of all accumulators. 
     *
     */
    struct asynchronous_results_merger : 
      details::threading::asynchronous_results_merger<
        data_t,
        details::merger_update_specific_dimension<data_t>,
        details::threading::initialisation_vector_specific_dimension<data_t>,
        details::s_dimension_update<typename data_t::value_type>
      >
    {
    public:
      typedef data_t result_t;

    private:
      typedef details::threading::initialisation_vector_specific_dimension<data_t> data_init_type;
      typedef details::merger_update_specific_dimension<data_t> merger_type;


      const size_t data_dimension;
      
    public:
      typedef details::threading::asynchronous_results_merger<
        result_t,
        merger_type, 
        data_init_type,
        details::s_dimension_update<typename data_t::value_type>
      > parent_type;
      typedef typename parent_type::lock_t lock_t;

      /*!Constructor
       *
       * @param dimension_ the number of dimensions of the vector to accumulate
       */
      asynchronous_results_merger(size_t data_dimension_) : 
        parent_type(data_init_type(data_dimension_)),
        data_dimension(data_dimension_)
      {}

    };





  public:
    /*!@brief Constructor
     *
     * Constructs an instance of the RobustPCA with trimming with the provided percentage of trimming.
     *
     * @param[in] trimming_percentage_ the percentage of data that should be trimmed from the lower and upper distributions.
     * @note By default the number of processors used for computation is set to 1.
     * The maximum size of the chunks is "infinite": each chunk will receive in that case the size of the data
     * divided by the number of running threads.
     */
    grassmann_pca_with_trimming(double trimming_percentage_ = 0) :
      random_init_op(details::fVerySmallButStillComputable, details::fVeryBigButStillComputable),
      trimming_percentage(trimming_percentage_),
      nb_processors(1),
      max_chunk_size(std::numeric_limits<size_t>::max()),
      nb_steps_pca(3)
    {
      assert(trimming_percentage_ >= 0 && trimming_percentage_ <= 1);
    }

    //! Sets the number of parallel tasks used for computing.
    bool set_nb_processors(size_t nb_processors_)
    {
      assert(nb_processors_ >= 1);
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

    //! Sets the number of iterations for the initial PCA like algorithm. 
    bool set_nb_steps_pca(size_t nb_steps)
    {
      nb_steps_pca = nb_steps;
      return true;
    }



    
    


    /*!@brief Performs the computation of the current subspace on the elements given by the two iterators.
     *
     * @tparam it_t an input forward iterator to input vectors points. Each element pointed by the underlying iterator should be iterable and
     *   should provide a vector point.
     * @tparam it_o_basisvectors_t an output iterator for storing the computed basis vectors. This iterator should model a forward output iterator.
     *
     * @param[in] max_iterations the maximum number of iterations at each dimension.
     * @param[in] max_dimension_to_compute the maximum number of data_dimension to compute in the PCA (only the first max_dimension_to_compute will be
     *            computed).
     * @param[in] it input iterator at the beginning of the data
     * @param[in] ite input iterator at the end of the data
     * @param[in] initial_guess if provided, the initial vectors will be initialized to this value.
     * @param[out] it_basisvectors an iterator on the beginning of the area where the detected basis vectors will be stored. The space should be at least max_dimension_to_compute.
     *
     * @returns true on success, false otherwise
     * @pre
     * - @c !(it >= ite)
     * - all the vectors given by the iterators pair should be of the same size (no check is performed).
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
      const size_t size_data(std::distance(it, ite));

      // size of the chunks.
      const size_t chunks_size = std::min(max_chunk_size, static_cast<size_t>(ceil(double(size_data)/nb_processors)));
      const size_t nb_chunks = (size_data + chunks_size - 1) / chunks_size;


      // number of dimensions of the data vectors
      const size_t number_of_dimensions = it->size();
      max_dimension_to_compute = std::min(max_dimension_to_compute, number_of_dimensions);


      // initial iterator on the output basis vectors
      it_o_basisvectors_t const it_output_basis_vector_beginning(it_basisvectors);
      it_o_basisvectors_t it_output_basis_vector_end(it_output_basis_vector_beginning);
      std::advance(it_output_basis_vector_end, max_dimension_to_compute);

      // the initialisation of mus
      {
        it_o_basisvectors_t it_basis(it_output_basis_vector_beginning);
        for(int i = 0; it_basis != it_output_basis_vector_end; ++it_basis, ++i)
        {
          *it_basis = initial_guess != 0 ? (*initial_guess)[i] : random_init_op(*it);
        }
      }
      if(!details::gram_schmidt_orthonormalisation(it_output_basis_vector_beginning, it_output_basis_vector_end, it_output_basis_vector_beginning, norm_op))
      {
        return false;
      }


      // preparing mu
      data_t mu(*it_basisvectors);
      assert(mu.size() == number_of_dimensions);
      


      // number of elements to keep
      const int K_elements = static_cast<int>(std::ceil(trimming_percentage*size_data/2));


      // preparing the ranges on which each processing thread will run.
      // the number of objects can be much more than the current number of processors, in order to
      // avoid waiting too long for a thread (better granularity) but involving a slight overhead in memory and
      // processing at the synchronization point.
      typedef s_robust_pca_trimmed_processor_inner_products async_processor_t;
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

          // the processor object for this new range
          async_processor_t &current_acc_object = v_individual_accumulators[i];

          // updating the dimension of the problem
          current_acc_object.set_data_dimensions(number_of_dimensions);

          // setting the number of elements to keep
          current_acc_object.set_nb_elements_to_keep(K_elements);

          b_result = current_acc_object.set_data_range(it_current_begin, it_current_end);
          if(!b_result)
          {
            return b_result;
          }
          
          // attaching the update object callbacks
          current_acc_object.connector_accumulator() = boost::bind(
              &asynchronous_results_merger::update, 
              &async_merger, 
              _1);

          current_acc_object.connector_counter() = boost::bind(&asynchronous_results_merger::notify, &async_merger);

          // updating the next 
          it_current_begin = it_current_end;
        }
      }




      // this matrix is a copy of the data. It is a convenient structure for storing the 
      // flipped vectors that will then be trimmed and accumulated. But this is a copy of the initial data:
      // if the memory pressure is too big, then something else should be found (like a vector of vector, being
      // potentially scattered in memory).
      boost::scoped_array<scalar_t> matrix_temp(new scalar_t[number_of_dimensions*size_data]);






      // for each dimension
      for(int current_dimension = 0; current_dimension < max_dimension_to_compute; current_dimension++, ++it_basisvectors)
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




        details::convergence_check<data_t> convergence_op(mu);

        int iterations = 0;
        for(; (!convergence_op(mu) && iterations < max_iterations) || iterations == 0; iterations++)
        {

          // reseting the merger object
          async_merger.init();

          // pushing the computation of the bounds
          for(int i = 0; i < v_individual_accumulators.size(); i++)
          {
            ioService.post(
              boost::bind(
                &async_processor_t::compute_data_matrix, 
                boost::ref(v_individual_accumulators[i]), 
                boost::cref(mu),
                matrix_temp.get() + i*chunks_size,
                size_data));
          }


          // waiting for completion (barrier)
          async_merger.wait_notifications(v_individual_accumulators.size());

          // clearing the notifications
          async_merger.init_notifications();

          // pushing the computation of the trimmed accumulation on each dimension
          for(int dim_to_compute = 0; dim_to_compute < number_of_dimensions; dim_to_compute++)
          {
            ioService.post(
              boost::bind(
                &async_processor_t::compute_bounded_accumulation, 
                boost::ref(v_individual_accumulators[0]), 
                dim_to_compute,
                size_data,
                matrix_temp.get() + dim_to_compute*size_data));
          }


          // waiting for completion (barrier)
          async_merger.wait_notifications(number_of_dimensions);
          

          // gathering the mus
          mu = async_merger.get_merged_result();
          
          // normalize mu on the sphere
          mu /= norm_op(mu);

        }



        // orthogonalisation against previous basis vectors
        bool renormalise(false);
        for(it_o_basisvectors_t it_mus(it_output_basis_vector_beginning); it_mus < it_basisvectors; ++it_mus)
        {
          mu -= boost::numeric::ublas::inner_prod(mu, *it_mus) * (*it_mus);
          renormalise = true;
        }
        if(renormalise)
        {
          mu /= norm_op(mu);
        }
        

        // mu is the basis vector of the current dimension, we store it in the output vector
        *it_basisvectors = mu;

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
                boost::cref(*it_basisvectors))); // this is not mu, since we are changing it before the process ends here
          }

          // this is to follow the matlab implementation, but the idea is the following:
          // each time we pick a new candidate vector, we project it to the orthogonal subspace of the previously computed 
          // basis vectors. This can be done in two ways:
          // 1. compute the projection on the orthogonal subspace of the current (or next) candidate
          // 2. compute the projection on the orthogonal subspace of the remainder elements
          //
          // in order to be consistent with the matlab implementation, the second choice is implemented here
          if(current_dimension+1 < max_dimension_to_compute)
          {
            it_o_basisvectors_t remainder(it_basisvectors);
            ++remainder;

            if(!details::gram_schmidt_orthonormalisation(it_output_basis_vector_beginning, it_output_basis_vector_end, remainder, norm_op))
            {
              return false;
            }

            mu = *remainder;
          }


          // wait for the workers
          async_merger.wait_notifications(v_individual_accumulators.size());

          
        }

      }

      return true;
    }

  };

}


#endif /* GRASSMANN_AVERAGES_PCA_TRIMMED_HPP__ */
