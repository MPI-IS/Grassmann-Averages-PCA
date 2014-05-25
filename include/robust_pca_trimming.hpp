// Copyright 2014, Max Planck Institute for Intelligent Systems.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

#ifndef ROBUST_PCA_TRIMMED_HPP__
#define ROBUST_PCA_TRIMMED_HPP__



/*!@file
 * Robust PCA functions, following the paper of Soren Hauberg.
 *
 * This file contains the implementation of the trimmed version. 
 */

// boost heaps
#include <boost/heap/pairing_heap.hpp>

// for the thread pools
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/signals2.hpp>
#include <numeric>

#include <boost/scoped_array.hpp>

// utilities
#include <include/private/utilities.hpp>


namespace robust_pca
{


  namespace details
  {
    
    //! Adaptation of merger_addition concept for accumulator and count at the same time.
    //! 
    //! The trimmed version of the robust pca algorithm may strip some element along each dimension. 
    //! In order to compute the @f$\mu@f$ properly, the count should also be transfered.
    template <class data_t, class count_t>
    struct merger_addition_with_count
    {
      typedef std::pair<data_t, count_t> input_t;
      bool operator()(input_t &current_state, input_t const& update_value) const
      {
        current_state.first += update_value.first;
        current_state.second += update_value.second;
        return true;
      }
    };

    struct s_dimension_update
    {
      size_t dimension;
      double value;
      size_t nb_elements;
    };

    //! Adaptation of merger_addition concept for accumulator and count at the same time.
    //! 
    //! The trimmed version of the robust pca algorithm may strip some element along each dimension. 
    //! In order to compute the @f$\mu@f$ properly, the count should also be transfered.
    template <class data_t, class count_t>
    struct merger_addition_specific_dimension
    {
      typedef std::pair<data_t, count_t> input_t;
      bool operator()(input_t &current_state, s_dimension_update const& update_value) const
      {
        current_state.first(update_value.dimension) += update_value.value;
        current_state.second(update_value.dimension) += update_value.nb_elements;
        return true;
      }
    };


    //! Adaptation of initialisation_vector_specific_dimension concept for accumulation and count.
    //! 
    //! See merger_addition_with_count.
    template <class data_t, class count_t>
    struct initialisation_vector_specific_dimension_with_count
    {
      size_t data_dimension;
      typedef typename data_t::value_type scalar_t;
      typedef typename count_t::value_type count_scalar_t;
      typedef std::pair<data_t, count_t> input_t;

      initialisation_vector_specific_dimension_with_count(size_t dimension) : data_dimension(dimension)
      {}

      bool operator()(input_t & current_state) const
      {
        current_state.first = boost::numeric::ublas::scalar_vector<scalar_t>(data_dimension, 0);
        current_state.second = boost::numeric::ublas::scalar_vector<count_scalar_t>(data_dimension, 0);
        return true;
      }
    };

    
  }


  /*!@brief Robust PCA subspace algorithm, with trimming
   *
   * This class implements the robust PCA using the Grassmanian averaging with trimming.
   *
   * @author Soren Hauberg, Raffi Enficiaud
   */
  template <class data_t, class norm_mu_t = details::norm2>
  struct robust_pca_with_trimming_impl
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

    //! Type of the element of data_t. 
    typedef typename data_t::value_type scalar_t;

    //! Type of the vector used for counting the element falling into the non-trimmed range.
    typedef boost::numeric::ublas::vector<size_t> count_vector_t;









    //! Performs a selective update of the accumulator, given the bounds computed previously
    template <class vector_bounds_t, class vector_number_elements_t>
    static void selective_acc_to_vector(
      const vector_bounds_t& lower_bounds, const vector_bounds_t& upper_bounds,
      const data_t &initial_data,
      bool sign,
      data_t &v_selective_accumulator,
      vector_number_elements_t& v_selective_acc_count)
    {
      for(size_t i = 0, j = initial_data.size(); i < j; i++)
      {
        typename data_t::value_type const v(sign ? initial_data[i] : -initial_data[i]);
        if(v < lower_bounds[i])
          continue;
        if(v > upper_bounds[i])
          continue;
        

        v_selective_accumulator[i] += v;
        v_selective_acc_count[i]++;
      }
    }






    //!@internal
    //!@brief Contains the logic for processing part of the accumulator
    struct s_robust_pca_trimmed_processor_inner_products
    {
    private:
      //! Scalar type of the data
      typedef typename data_t::value_type scalar_t;
      
      typedef details::s_dimension_update accumulator_element_t;

      
      size_t nb_elements;                 //!< The size of the current dataset
      size_t data_dimension;              //!< The dimension of the data
      size_t nb_elements_to_keep;         //!< The number of elements to keep.
      
      typedef boost::signals2::signal<void ()>
        connector_counter_t;
      connector_counter_t signal_counter;

      typedef boost::signals2::signal<void (accumulator_element_t const&)> 
        connector_accumulator_dimension_t;
      connector_accumulator_dimension_t signal_acc_dimension;

 
      //! The matrix containing a copy of the data
      scalar_t *p_c_matrix;
      
      //! The result of the inner products
      std::vector<double> inner_prod_results;
      std::vector<scalar_t> current_dimension_sort;
  
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
          double* current_line = p_c_matrix + column;
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
      }

      //! Sets the number of element to keep in the upper and lower distributions.
      void set_nb_elements_to_keep(int nb_elements_to_keep_)
      {
        assert(nb_elements_to_keep_ >= 0);
        nb_elements_to_keep = nb_elements_to_keep_;
      }


      
      
      void compute_data_matrix(data_t const &mu, typename data_t::value_type* p_out, size_t padding)
      {
        // updates the internal inner products
        compute_inner_products(mu);
        
        // the current line is spans a particular dimension
        double *current_line = p_c_matrix;

        // this spans the inner product results for all dimensions

        std::vector<double> v_mult(nb_elements);
        for(size_t element(0); element < nb_elements; element++)
        {
          v_mult[element] = inner_prod_results[element] >= 0 ? 1 : -1;
        }
        double const * const out = &v_mult[0];

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
      
      
      void compute_bounded_accumulation(size_t dimension, size_t nb_total_elements, typename data_t::value_type* p_data)
      {

#if 0
        std::nth_element(p_data, p_data + nb_elements_to_keep, p_data + nb_total_elements);
        const typename data_t::value_type min_value = p_data[nb_elements_to_keep];
        std::nth_element (p_data, p_data + nb_elements_to_keep, p_data + nb_total_elements, std::greater<typename data_t::value_type>());
        const typename data_t::value_type max_value = p_data[nb_elements_to_keep];


        size_t nb_elements_acc = nb_total_elements - 2*nb_elements_to_keep;

        double acc = 0;
        for(int i = 0; i < nb_total_elements; i++)
        {
          if(p_data[i] >= min_value && p_data[i] <= max_value)
            acc += p_data[i];
        }
#endif

        std::nth_element(p_data, p_data + nb_elements_to_keep, p_data + nb_total_elements);
        std::nth_element(p_data + nb_elements_to_keep+1, p_data + nb_total_elements - nb_elements_to_keep-1, p_data + nb_total_elements);
        double acc = std::accumulate(p_data + nb_elements_to_keep, p_data + nb_total_elements - nb_elements_to_keep, 0.);
        
        accumulator_element_t result;
        result.dimension = dimension;
        result.value = acc;
        result.nb_elements = nb_total_elements - 2*nb_elements_to_keep;
        // signals the update
        signal_acc_dimension(result);
        // signals the main merger
        signal_counter();
      }
      
      //! Project the data onto the orthogonal subspace of the provided vector
      void project_onto_orthogonal_subspace(data_t const &mu)
      {
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
     * @brief Merges the result of all workers and signals the results to the main thread.
     *
     * The purpose of this class is to add the computed accumulator of each thread to the final result
     * which contains the sum of all accumulators. 
     *
     */
    struct asynchronous_results_merger : 
      details::threading::asynchronous_results_merger<
        std::pair<data_t, count_vector_t>,
        details::merger_addition_specific_dimension<data_t, count_vector_t>,
        details::initialisation_vector_specific_dimension_with_count<data_t, count_vector_t>,
        details::s_dimension_update
      >
    {
    public:
      typedef std::pair<data_t, count_vector_t> result_t;

    private:
      typedef details::initialisation_vector_specific_dimension_with_count<data_t, count_vector_t> data_init_type;
      typedef details::merger_addition_specific_dimension<data_t, count_vector_t> merger_type;


      const size_t data_dimension;
      
    public:
      typedef details::threading::asynchronous_results_merger<
        result_t,
        merger_type, 
        data_init_type,
        details::s_dimension_update> parent_type;
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
    robust_pca_with_trimming_impl(double trimming_percentage_ = 0) :
      random_init_op(details::fVerySmallButStillComputable, details::fVeryBigButStillComputable),
      trimming_percentage(trimming_percentage_),
      nb_processors(1),
      max_chunk_size(std::numeric_limits<size_t>::max())
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
    
    


    /*! Performs the computation of the current subspace on the elements given by the two iterators.
     *  @tparam it_t an input forward iterator to input vectors points. Each element pointed by the underlying iterator should be iterable and
     *   should provide a vector point.
     *  @tparam it_o_projected_vectors output forward iterator pointing on a container of vector points.
     *  @tparam it_norm_t an output iterator on weights/norms of the vectors. The output elements should be numerical (norm output)
     *
     * @param[in] max_iterations the maximum number of iterations at each dimension.
     * @param[in] max_dimension_to_compute the maximum number of data_dimension to compute in the PCA (only the first max_dimension_to_compute will be
     *            computed).
     * @param[in] it input iterator at the beginning of the data
     * @param[in] ite input iterator at the end of the data
     * @param[in, out] it_projected
     * @param[in] initial_guess if provided, the initial vector will be initialized to this value.
     * @param[out] it_eigenvectors an iterator on the beginning of the area where the detected eigenvectors will be stored. The space should be at least max_dimension_to_compute.
     *
     * @returns true on success, false otherwise
     * @pre
     * - @c !(it >= ite)
     * - all the vectors given by the iterators pair should be of the same size (no check is performed).
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
      const size_t size_data(std::distance(it, ite));

      // size of the chunks.
      const size_t chunks_size = std::min(max_chunk_size, static_cast<size_t>(size_data/nb_processors));
      const size_t nb_chunks = (size_data + chunks_size - 1) / chunks_size;


      // number of dimensions of the data vectors
      const size_t number_of_dimensions = it->size();
      max_dimension_to_compute = std::min(max_dimension_to_compute, number_of_dimensions);


      // initial iterator on the output eigenvectors
      it_o_eigenvalues_t const it_output_eigen_vector_beginning(it_eigenvectors);
      it_o_eigenvalues_t it_output_eigen_vector_end(it_output_eigen_vector_beginning);
      std::advance(it_output_eigen_vector_end, max_dimension_to_compute);

      // the initialisation of mus
      {
        it_o_eigenvalues_t it_eigen(it_output_eigen_vector_beginning);
        for(int i = 0; it_eigen != it_output_eigen_vector_end; ++it_eigen, ++i)
        {
          *it_eigen = initial_guess != 0 ? (*initial_guess)[i] : random_init_op(*it);
        }
      }
      if(!details::gram_schmidt_orthonormalisation(it_output_eigen_vector_beginning, it_output_eigen_vector_end, it_output_eigen_vector_beginning, norm_op))
      {
        return false;
      }


      // preparing mu
      data_t mu(*it_eigenvectors);
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
          current_acc_object.connector_accumulator().connect(
            boost::bind(
              &asynchronous_results_merger::update, 
              &async_merger, 
              _1));

          current_acc_object.connector_counter().connect(
            boost::bind(&asynchronous_results_merger::notify, &async_merger));

          // updating the next 
          it_current_begin = it_current_end;
        }
      }




      // this matrix is a copy of the data. It is a convenient structure for storing the 
      // flipped vectors that will then be trimmed and accumulated. But this is a copy of the initial data:
      // if the memory pressure is too big, then something else should be found (like a vector of vector, being
      // potentially scattered in memory).
      boost::scoped_array<double> matrix_temp(new double[number_of_dimensions*size_data]);






      // for each dimension
      for(int current_dimension = 0; current_dimension < max_dimension_to_compute; current_dimension++, ++it_eigenvectors)
      {

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
          mu = async_merger.get_merged_result().first;
          count_vector_t const& count_vector = async_merger.get_merged_result().second;
          // divide each of the acc by acc_counts, and then take the norm
          for(int i = 0; i < number_of_dimensions; i++)
          {
            assert(count_vector[i]);
            mu[i] /= count_vector[i];
          }

          // normalize mu on the sphere
          mu /= norm_op(mu);

        }



        // orthogonalisation against previous eigenvectors
        bool renormalise(false);
        for(it_o_eigenvalues_t it_mus(it_output_eigen_vector_beginning); it_mus < it_eigenvectors; ++it_mus)
        {
          mu -= boost::numeric::ublas::inner_prod(mu, *it_mus) * (*it_mus);
          renormalise = true;
        }
        if(renormalise)
        {
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

          // this is to follow the matlab implementation, but the idea is the following:
          // each time we pick a new candidate vector, we project it to the orthogonal subspace of the previously computed 
          // eigenvectors. This can be done in two ways:
          // 1. compute the projection on the orthogonal subspace of the current (or next) candidate
          // 2. compute the projection on the orthogonal subspace of the remainder elements
          //
          // in order to be consistent with the matlab implementation, the second choice is implemented here
          if(current_dimension+1 < max_dimension_to_compute)
          {
            it_o_eigenvalues_t remainder(it_eigenvectors);
            ++remainder;

            if(!details::gram_schmidt_orthonormalisation(it_output_eigen_vector_beginning, it_output_eigen_vector_end, remainder, norm_op))
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


#endif /* ROBUST_PCA_TRIMMED_HPP__ */