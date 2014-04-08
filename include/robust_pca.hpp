// Copyright 2014, Max Planck Institute for Intelligent Systems.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

#ifndef ROBUST_PCA_HPP__
#define ROBUST_PCA_HPP__

/*!@file
 * Robust PCA functions, following the paper of Soren Hauberg.
 *
 * @note These implementations assume the existence of boost somewhere. Currently, it is just used for 
 * computing the norms of differences and to generate some random data (which is now TR1 and C++0X).
 */


#include <boost/numeric/conversion/bounds.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <vector>


// for the trimmed version
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/extended_p_square_quantile.hpp>


// for the thread pools
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/signals2.hpp>

namespace robust_pca
{
 
  //! Wrapper object for infinity/max @f$\ell_\infty@f$ norm.
  struct norm_infinity
  {
    template <class vector_t>
    double operator()(vector_t const& v) const
    {
      return boost::numeric::ublas::norm_inf(v);
    }
  };


  //! Returns the square of the @f$\ell_2@f$ norm.
  struct norm_ell2_square
  {
    template <class vector_t>
    double operator()(vector_t const& v) const
    {
      double acc(0);
      for(typename vector_t::const_iterator it(v.begin()), ite(v.end());
          it < ite;
          ++it)
      {
        typename vector_t::const_reference v(*it);
        acc += v * v;
      }
      return acc;
    }
  };

  //! Returns the @f$\ell_2@f$ norm of a vector.
  struct norm2
  {
    typedef double result_type;
    norm_ell2_square op;
    template <class vector_t>
    result_type operator()(vector_t const& v) const
    {
      return std::sqrt(op(v));
    }
  };


  /* @brief Checks the convergence of a numerical scheme.
   *
   * @tparam data_t should be default constructible and constructible with one parameter.
   *
   * The convergence is assumed when the norm between two subsequent states is less than a certain @f$epsilon@f$.
   */
  template <class data_t, class norm_t = norm_infinity>
  struct convergence_check
  {
    const double epsilon;
    norm_t norm_comparaison;
    data_t previous_state;
    convergence_check(data_t const& prev, double epsilon_ = 1E-5) : 
      epsilon(epsilon_), 
      previous_state(prev)
    {}

    //! Returns true on convergence
    bool operator()(data_t const& current_state)
    {
      bool ret = norm_comparaison(current_state - previous_state) < epsilon;
      previous_state = current_state;
      return ret;
    }

  };


  // some issues with the random number generator
  const double fVeryBigButStillComputable = 1E10;
  const double fVerySmallButStillComputable = -1E10;


  /*@brief Used to initialize the initial guess to some random data.
   *
   *@tparam data_t the type of the data returned. It should model a vector:
   * - default, copy constructible and constructible with a size
   * - has a member @c size returning a value of type @c data_t::Size_type
   * - has a field @c data_t::Value_type
   */
  template <
    class data_t, 
    class random_distribution_t = 
      typename boost::mpl::if_<
        boost::is_floating_point<typename data_t::value_type>,
        boost::random::uniform_real_distribution<typename data_t::value_type>,
        boost::random::uniform_int_distribution<typename data_t::value_type>
      >::type,
    class random_number_generator_t = boost::random::mt19937
  >
  struct random_data_generator
  {
    typedef typename data_t::value_type value_type;
    mutable random_number_generator_t rng;
    value_type min_bound;
    value_type max_bound;
    
    //! Default construction
    random_data_generator(
      value_type min_value_ = boost::numeric::bounds<value_type>::lowest(),
      value_type max_value_ = boost::numeric::bounds<value_type>::highest()) : 
      rng(), 
      min_bound(min_value_), 
      max_bound(max_value_)
    {}

    //! Constructs from the specified seeded random number generator.
    random_data_generator(
      random_number_generator_t const &r,
      value_type min_value_ = boost::numeric::bounds<value_type>::lowest(),
      value_type max_value_ = boost::numeric::bounds<value_type>::highest()) : 
      rng(r),
      min_bound(min_value_), 
      max_bound(max_value_)
    {}


    data_t operator()(const data_t& v) const
    {
      data_t out(v.size());
      random_distribution_t dist(min_bound, max_bound);
      for(typename data_t::size_type i(0), j(v.size()); i < j; i++)
      {
        out[i] = dist(rng);
      }
      return out;
    }
  };








  /*!@brief Robust PCA subspace algorithm
   *
   * This class implements the robust PCA using the Grassmanian averaging. This is the naive implementation which is
   * suitable for small datasets. 
   * @todo add reference
   * @author Soren Hauberg, Raffi Enficiaud
   */
  template <class data_t, class norm_2_t = norm2>
  struct robust_pca_impl
  {
  private:
    random_data_generator<data_t> random_init_op;
    norm_2_t norm_op;

    //! Number of parallel tasks that will be used for computing.
    int nb_processors;

    //! Maximal size of a chunk (infinity by default).
    size_t max_chunk_size;

    //!@internal
    //!@brief Contains the logic for processing part of the accumulator
    template <class container_iterator_t>
    struct s_accumulator_processor
    {
    private:
      container_iterator_t begin, end;
      size_t nb_elements;
      data_t accumulator;
      std::vector<bool> v_signs;
      int data_dimension;

      // this is to send an update of the value of mu to all listeners
      // the connexion should be managed externally
      typedef boost::signals2::signal<void (data_t const&)> connector_accumulator_t;
      typedef boost::signals2::signal<void ()>              connector_counter_t;
      connector_accumulator_t signal_acc;
      connector_counter_t signal_counter;


    public:
      s_accumulator_processor() : nb_elements(0), data_dimension(0)
      {
      }


      //! Sets the data range
      bool set_data_range(container_iterator_t const &b, container_iterator_t const& e)
      {
        begin = b;
        end = e;
        nb_elements = std::distance(b, e);
        assert(nb_elements > 0);
        v_signs.resize(nb_elements);
        return true;
      }

      //! Sets the dimension of the problem
      //! @pre data_dimensions_ strictly positive
      void set_data_dimensions(int data_dimensions_)
      {
        data_dimension = data_dimensions_;
        assert(data_dimension > 0);
      }

      //! Returns the connected object that will receive the notification of the updated accumulator.
      connector_accumulator_t& accumulator_connector()
      {
        return signal_acc;
      }

      //! Returns the connected object that will receive the notification of the end of the current process.
      connector_counter_t& counter_connector()
      {
        return signal_counter;
      }
      

      //! Initialises the accumulator and the signs vector from the first mu
      void initial_accumulation(data_t const &mu)
      {
        accumulator = data_t(data_dimension, 0);
        std::vector<bool>::iterator itb(v_signs.begin());
        
        container_iterator_t it_data(begin);

        // first iteration, we store the signs
        for(size_t s = 0; s < nb_elements; ++it_data, ++itb, s++)
        {
          typename container_iterator_t::reference current_data = *it_data;
          bool sign = boost::numeric::ublas::inner_prod(current_data, mu) >= 0;
          *itb = sign;

          if(sign)
          {
            accumulator += current_data;
          }
          else 
          {
            accumulator -= current_data;
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
        
        container_iterator_t it_data(begin);

        // first iteration, we store the signs
        for(size_t s = 0; s < nb_elements; ++it_data, ++itb, s++)
        {
          typename container_iterator_t::reference current_data = *it_data;

          bool sign = boost::numeric::ublas::inner_prod(current_data, mu) >= 0;
          if(sign != *itb)
          {
            // update the value of the accumulator according to sign change
            *itb = sign;

            if(sign)
            {
              accumulator += 2 * current_data;
            }
            else
            {
              accumulator -= 2 * current_data;
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
        container_iterator_t it_data(begin);

        // update of vectors in the orthogonal space, and update of the norms at the same time. 
        for(size_t s = 0; s < nb_elements; ++it_data, s++)
        {
          typename container_iterator_t::reference current_vector = *it_data;
          current_vector -= boost::numeric::ublas::inner_prod(mu, current_vector) * mu;
        }

        signal_counter();
      }

    };


    /*! Accumulation gathering the result of all workers.
     *
     * When the vector are of high dimension, it should be interesting to also have
     * the update of the final accumulator within each threads. 
     */
    struct asynchronous_addition : boost::noncopyable
    {
    private:
      mutable boost::mutex internal_mutex;
      data_t current_value;
      volatile int nb_updates;
      const int data_dimension;

    public:

      /*!Constructor
       *
       * @param dimension_ the number of dimensions of the vector to accumulate
       */
      asynchronous_addition(int data_dimension_) : data_dimension(data_dimension_)
      {}

      void init()
      {
        clear_accumulation();
        clear_counter();
      }

      //! Initialises the internal state of the accumulator
      void clear_accumulation()
      {
        current_value = data_t(data_dimension, 0);
      }

      //! Initialises the internal state of the counter
      void clear_counter()
      {
        nb_updates = 0;
      }

      /*! Function receiving the updated value of the vectors to accumulate
       *  from different threads.
       * 
       *  @note The call is thread safe.
       */
      void update(data_t const& updated_value)
      {
        boost::lock_guard<boost::mutex> guard(internal_mutex);
        current_value += updated_value;
      }

      void update()
      {
        boost::lock_guard<boost::mutex> guard(internal_mutex);
        nb_updates ++;
      }

      //! Returns the number of updates received so far.
      int get_nb_updates() const
      {
        // to avoid possibly corrupted data, but I think for int this is unnecessary
        boost::lock_guard<boost::mutex> guard(internal_mutex); 
        return nb_updates;
      }

      //! Returns the current value 
      data_t const& get_accumulated_data() const
      {
        return current_value;
      }
    };


    //! Ensures the proper stop of the processing pool
    struct safe_stop
    {
      boost::asio::io_service& io_service;
      boost::thread_group& thread_group;

      safe_stop(boost::asio::io_service& ios, boost::thread_group& tg) : io_service(ios), thread_group(tg)
      {
      }

      ~safe_stop()
      {
        io_service.stop();
        thread_group.join_all();
      }
    };




  public:
    robust_pca_impl() : 
      random_init_op(fVerySmallButStillComputable, fVeryBigButStillComputable), 
      nb_processors(1),
      max_chunk_size(std::numeric_limits<size_t>::max())
    {}


    //! Sets the number of parallel tasks used for computing.
    bool set_nb_processors(int nb_processors_)
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
      max_chunk_size = chunk_size;
      return true;
    }



    /*! Performs the computation of the current subspace on the elements given by the two iterators.
     *
     * @tparam it_t an input forward iterator to input vectors points. Each element pointed by the underlying iterator should be iterable and
     *  should provide a vector point.
     * @tparam it_o_projected_vectors output random access iterator pointing on a container of vector points. 
     * @tparam it_norm_t an output iterator on weights/norms of the vectors. The output elements should be numerical (norm output)
     *
     * @param[in] max_iterations the maximum number of iterations at each dimension. 
     * @param[in] max_dimension_to_compute the maximum number of data_dimension to compute in the PCA (only the first max_dimension_to_compute will be 
     *            computed).
     * @param[in] it input iterator at the beginning of the data
     * @param[in] ite input iterator at the end of the data
     * @param[in, out] it_norm_out input read-write iterator at the beginning of the computed norms. The iterator should be able to address
     *            as many element as there is in between it and ite (ie. @c std::distance(it, ite)
     * @param[in, out] it_projected
     * @param[in] initial_guess if provided, the initial vectors will be initialized to this value. The size of the pointed container should be at least max_dimension_to_compute.
     * @param[out] it_eigenvectors an iterator on the beginning of the area where the detected eigenvectors will be stored. The space should be at least max_dimension_to_compute.
     *
     * @returns true on success, false otherwise
     * @pre 
     * - @c !(it >= ite)
     * - all the vectors given by the iterators pair should be of the same size (no check is performed).
     *
     * @note the iterator it_o_projected_vectors should implement random access to the elements.
     */
    template <class it_t, class it_o_projected_vectors, class it_o_eigenvalues_t/*, class it_norm_t*/>
    bool batch_process(
      const int max_iterations,
      int max_dimension_to_compute,
      it_t const it, 
      it_t const ite, 
      it_o_projected_vectors const it_projected,
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
      safe_stop worker_lock_guard(ioService, threadpool);

      // this is exactly the number of processors
      boost::asio::io_service::work work(ioService);
      for(int i = 0; i < nb_processors; i++)
      {
        threadpool.create_thread(boost::bind(&boost::asio::io_service::run, &ioService));
      }

      // contains the number of elements. In case the iterator is random access, could be deduced simply 
      // by a call to distance.
      size_t size_data(0);


      // The vectors are copied into the temporary container. 
      // TODO this is not necessary in case if the it_o_projected_vectors and it_t return both the same type.
      it_o_projected_vectors it_tmp_projected(it_projected);
      for(it_t it_copy(it); it_copy != ite; ++it_copy, ++it_tmp_projected, size_data++)
      {
        *it_tmp_projected = *it_copy;
      }

      assert(size_data == std::distance(it, ite));

      // size of the chunks.
      const size_t chunks_size = std::min(max_chunk_size, static_cast<size_t>(size_data/nb_processors));
      const size_t nb_chunks = (size_data + chunks_size - 1) / chunks_size;

      

      // the first element is used for the init guess because for dynamic std::vector like element, the size is needed.
      data_t mu(initial_guess != 0 ? (*initial_guess)[0] : random_init_op(*it));
      typename norm_2_t::result_type norm_mu(norm_op(mu)); // normalizing
      mu /= norm_mu;

      const int number_of_dimensions = static_cast<int>(mu.size());
      max_dimension_to_compute = std::min(max_dimension_to_compute, number_of_dimensions);
      
      int iterations = 0;


      // preparing the ranges on which each processing thread will run.
      // the number of objects can be much more than the current number of processors, in order to
      // avoid waiting too long for a thread (better granularity) but involving a slight overhead in memory and
      // processing at the synchronization point.
      typedef s_accumulator_processor<it_o_projected_vectors> individual_accumulators_t;
      std::vector<individual_accumulators_t> v_individual_accumulators(nb_chunks);

      asynchronous_addition async_add_object(number_of_dimensions);

      {
        bool b_result;
        it_o_projected_vectors it_current_begin(it_projected);
        for(int i = 0; i < nb_chunks; i++)
        {
          // setting the range
          it_o_projected_vectors it_current_end;
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

          individual_accumulators_t &current_acc_object = v_individual_accumulators[i];

          b_result = current_acc_object.set_data_range(it_current_begin, it_current_end);
          if(!b_result)
          {
            return b_result;
          }

          // updating the dimension of the problem
          current_acc_object.set_data_dimensions(number_of_dimensions);

          // attaching the update object callbacks
          current_acc_object.accumulator_connector().connect(boost::bind(&asynchronous_addition::update, &async_add_object, _1));
          current_acc_object.counter_connector().connect(boost::bind(&asynchronous_addition::update, &async_add_object));

          // updating the next 
          it_current_begin = it_current_end;
        }
      }



      // for each dimension
      for(int current_dimension = 0; current_dimension < max_dimension_to_compute; current_dimension++, ++it_eigenvectors)
      {

        convergence_check<data_t> convergence_op(mu);

        data_t previous_mu(mu);

        // reseting the final accumulator
        async_add_object.init();

        // pushing the initialisation of the mu and sign vectors to the pool
        for(int i = 0; i < v_individual_accumulators.size(); i++)
        {
          ioService.post(boost::bind(&individual_accumulators_t::initial_accumulation, boost::ref(v_individual_accumulators[i]), boost::cref(previous_mu)));
        }

        // waiting for completion (barrier)
        while(async_add_object.get_nb_updates() < v_individual_accumulators.size())
        {
          boost::this_thread::yield();
        }

        // gathering the first mu
        mu = async_add_object.get_accumulated_data();
        mu /= norm_op(mu);


        // other iterations as usual
        for(iterations = 1; !convergence_op(mu) && iterations < max_iterations; iterations++)
        {
          previous_mu = mu;
          it_o_projected_vectors it_tmp_projected(it_projected);

          // reseting the final accumulator
          async_add_object.init();

          // pushing the update of the mu (and signs)
          for(int i = 0; i < v_individual_accumulators.size(); i++)
          {
            ioService.post(boost::bind(&individual_accumulators_t::update_accumulation, boost::ref(v_individual_accumulators[i]), boost::cref(previous_mu)));
          }

          // waiting for completion (barrier)
          while(async_add_object.get_nb_updates() < v_individual_accumulators.size())
          {
            boost::this_thread::yield();
          }

          // gathering the mus
          mu = async_add_object.get_accumulated_data();
          mu /= norm_op(mu);
        }


        // mu is the eigenvector of the current dimension, we store it in the output vector
        *it_eigenvectors = mu;

        // projection onto the orthogonal subspace
        if(current_dimension < max_dimension_to_compute - 1)
        {

          async_add_object.clear_counter();

          // pushing the update of the mu (and signs)
          for(int i = 0; i < v_individual_accumulators.size(); i++)
          {
            ioService.post(boost::bind(&individual_accumulators_t::project_onto_orthogonal_subspace, boost::ref(v_individual_accumulators[i]), boost::cref(*it_eigenvectors)));
          }

          mu = initial_guess != 0 ? (*initial_guess)[current_dimension+1] : random_init_op(*it);

          while(async_add_object.get_nb_updates() < v_individual_accumulators.size())
          {
            boost::this_thread::yield();
          }

        }
        
      }


      // stopping the pool is done in the destruction of worker_lock_guard



      return true;
    }
  };












  /*!@brief Robust PCA subspace algorithm, with trimming
  *
  * This class implements the robust PCA using the Grassmanian averaging. This is the naive implementation which is
  * suitable for small datasets.
  *
  * @author Soren Hauberg, Raffi Enficiaud
  */
  template <class data_t, class norm_2_t = norm2>
  struct robust_pca_with_trimming_impl
  {
  private:
    random_data_generator<data_t> random_init_op;
    norm_2_t norm_op;
    double min_trim_percentage;
    double max_trim_percentage;

    template <class vector_acc_t>
    void apply_quantile_to_vector(const data_t& current_data, bool sign, vector_acc_t& v_acc) const
    {
      typename vector_acc_t::iterator it_acc_features(v_acc.begin());
      for(typename data_t::const_iterator it(current_data.begin()), ite(current_data.end()); it < ite; ++it, ++it_acc_features)
      {
        typename data_t::value_type v(sign ? *it : -(*it));
        (*it_acc_features)(v);
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
        typename data_t::value_type const v(sign ? initial_data[i] : -initial_data[i]);
        if(v < lower_bounds[i])
          continue;
        if(v > upper_bounds[i])
          continue;
        

        v_selective_accumulator[i] += v;
        v_selective_acc_count[i]++;
      }
    }

  public:
    robust_pca_with_trimming_impl(double min_trim_percentage_ = 0, double max_trim_percentage_ = 1.) :
      random_init_op(fVerySmallButStillComputable, fVeryBigButStillComputable),
      min_trim_percentage(min_trim_percentage_),
      max_trim_percentage(max_trim_percentage_)
    {
      assert(min_trim_percentage_ < max_trim_percentage_);
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
      typename norm_2_t::result_type norm_mu(norm_op(mu));
      mu /= norm_mu;


      
      // Accumulator object containing the required quantiles. 
      std::vector<double> probs(2);
      probs[0] = min_trim_percentage;
      probs[1] = max_trim_percentage;
      const accumulator_t quantil_obj(extended_p_square_probabilities = probs);
      
      // vector of accumulator objects
      std::vector<accumulator_t> v_acc(number_of_dimensions, quantil_obj);

      // vector of valid bounds
      std::vector<double> v_min_threshold(number_of_dimensions);
      std::vector<double> v_max_threshold(number_of_dimensions);


      // copy of the vectors. This should be avoided if the containers are of the same type and for the first iteration. 
      // The vectors are copied into the temporary container. During the copy, the percentiles are computed.
      it_o_projected_vectors it_tmp_projected(it_projected);
      for(it_t it_copy(it); it_copy != ite; ++it_copy , ++it_tmp_projected, size_data++)
      {
        *it_tmp_projected = *it_copy;
      }



      // initial iterator on the output eigenvectors
      it_o_eigenvalues_t const it_output_eigen_vector_beginning(it_eigenvectors);



      // for each dimension
      for(int current_dimension = 0; current_dimension < max_dimension_to_compute; current_dimension++, ++it_eigenvectors)
      {

        convergence_check<data_t> convergence_op(mu);



        for(int iterations = 0; (!convergence_op(mu) && iterations < max_iterations) || iterations == 0; iterations++)
        {

          // first pass on the data, we compute the bounds
          v_acc.assign(number_of_dimensions, quantil_obj);

          {
            it_o_projected_vectors it_tmp_projected(it_projected);

            for(size_t s = 0; s < size_data; ++it_tmp_projected, s++)
            {
              typename it_o_projected_vectors::reference current_data = *it_tmp_projected;

              bool sign = boost::numeric::ublas::inner_prod(current_data, mu) >= 0;

              // updating the bounds
              apply_quantile_to_vector(current_data, sign, v_acc);
            }
          }


          // extracting the bounds
          for(int i = 0; i < number_of_dimensions; i++)
          {
            v_min_threshold[i] = quantile(v_acc[i], quantile_probability = min_trim_percentage); //extract_result</*tag::p_square_quantile*/tag::min>(v_acc[i]);
            v_max_threshold[i] = quantile(v_acc[i], quantile_probability = max_trim_percentage); //extract_result</*tag::p_square_quantile*/tag::max>(v_acc[i]);// quantile(v_acc[i], quantile_probability = max_trim_percentage);
          }

          //v_acc.clear();
          


          data_t previous_mu = mu;
          data_t acc = data_t(number_of_dimensions, 0);
          std::vector<size_t> acc_counts(number_of_dimensions, 0);

          it_o_projected_vectors it_tmp_projected(it_projected);

          for(size_t s = 0; s < size_data; ++it_tmp_projected, s++)
          {
            typename it_o_projected_vectors::reference current_data = *it_tmp_projected;

            bool sign = boost::numeric::ublas::inner_prod(current_data, previous_mu) >= 0;

            // computing the bounds for the next round
            //apply_quantile_to_vector(current_data, sign, v_acc);

            // trimming is applied to the accumulator + taking into account the sign.
            selective_acc_to_vector(v_min_threshold, v_max_threshold, current_data, sign, acc, acc_counts);

          }


          // divide each of the acc by acc_counts, and then take the norm
          for(int i = 0; i < number_of_dimensions; i++)
          {
            assert(acc_counts[i]);
            mu[i] = acc[i] / acc_counts[i];
          }

          // normalize mu on the sphere
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

          mu = initial_guess != 0 ? (*initial_guess)[current_dimension+1] : random_init_op(*it);


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






  /*!@brief Robust PCA subspace algorithm, with trimming
  *
  * This class implements the robust PCA using the Grassmanian averaging. This is the naive implementation which is
  * suitable for small datasets.
  *
  * @author Soren Hauberg, Raffi Enficiaud
  */
  template <class data_t, class norm_2_t = norm2>
  struct robust_pca_with_stable_trimming_impl
  {
  private:
    random_data_generator<data_t> random_init_op;
    norm_2_t norm_op;
    double min_trim_percentage;
    double max_trim_percentage;

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
    robust_pca_with_stable_trimming_impl(double min_trim_percentage_ = 0, double max_trim_percentage_ = 1.) :
      random_init_op(fVerySmallButStillComputable, fVeryBigButStillComputable),
      min_trim_percentage(min_trim_percentage_),
      max_trim_percentage(max_trim_percentage_)
    {
      assert(min_trim_percentage_ < max_trim_percentage_);
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
      typename norm_2_t::result_type norm_mu(norm_op(mu));
      mu /= norm_mu;



      // Accumulator object containing the required quantiles. 
      std::vector<double> probs = { min_trim_percentage, max_trim_percentage };
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
        v_min_threshold[i] = quantile(v_acc[i], quantile_probability = min_trim_percentage); //extract_result</*tag::p_square_quantile*/tag::min>(v_acc[i]);
        v_max_threshold[i] = quantile(v_acc[i], quantile_probability = max_trim_percentage); //extract_result</*tag::p_square_quantile*/tag::max>(v_acc[i]);// quantile(v_acc[i], quantile_probability = max_trim_percentage);
      }
      
      

      // initial iterator on the output eigenvectors
      it_o_eigenvalues_t const it_output_eigen_vector_beginning(it_eigenvectors);

      // storing the signs for fast updates
      std::vector<bool> v_signs(size_data);

      // for each dimension
      for(int current_dimension = 0; current_dimension < max_dimension_to_compute; current_dimension++, ++it_eigenvectors)
      {

        convergence_check<data_t> convergence_op(mu);

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

#endif /* ROBUST_PCA_HPP__ */
