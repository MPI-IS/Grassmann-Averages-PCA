// Copyright 2014, Max Planck Institute for Intelligent Systems.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

#ifndef ROBUST_PCA_UTILITIES_HPP__
#define ROBUST_PCA_UTILITIES_HPP__

/*!@file
 * Robust PCA companion functions.
 *
 * This file contains some utility function for multithreading, norm computation, convergence check
 * 
 */


#include <boost/numeric/conversion/bounds.hpp>


#include <boost/asio/io_service.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/recursive_mutex.hpp>


#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector.hpp>


namespace robust_pca
{


  namespace details
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
    
    
    template <class T>
    double inner_prod(const T& v1, const T& v2)
    {
      double prod(0);
      assert(v1.size() == v2.size());
      double const *p1 = &v1(0);
      double const *p2 = &v1(0) + v1.size();
      double const *p3 = &v2(0);
      for(; p1 < p2; ++p1, ++p3)
      {
        prod += (*p1) * (*p3);
      }
      return prod;
    }



    /*!@brief Gram Schmidt orthonormalisation of a collection of vectors.
     * @tparam it_t iterator on the collection of vectors. Should model a forward input iterator.
     * @tparam norm_t type of the norm operator.
     * 
     * @param it beginning of the collection of vectors
     * @param ite end of the collection of vectors
     * @param start first element of the collection to be orthonormalized. start should be inside the range given by it and ite. 
     */
    template <class it_t, class norm_t>
    bool gram_schmidt_orthonormalisation(it_t it, it_t ite, it_t start, norm_t const &norm_op)
    {
      
      if(start == it)
      {
        *start /= norm_op(*start);
        ++start;
      }

      it_t previous(start);
              
      for(; start != ite; ++previous, ++start)
      {
        typename it_t::reference current = *start;
        for(it_t it_orthonormalised_element(it); it_orthonormalised_element < previous; ++it_orthonormalised_element)
        {
          current -= boost::numeric::ublas::inner_prod(current, *it_orthonormalised_element) * (*it_orthonormalised_element);
        }
        current /= norm_op(current);
              
      }
      return true;
    }



    /*!@brief Checks the convergence of a numerical scheme.
     *
     * @tparam data_t: type of the data.
     * @tparam norm_t: the norm used in order to compare the closeness of two successive results.
     *
     * The type of the data should meet the following requirements:
     * - data_t should be copy constructible and assignable.
     * - operator- is defined between two instances of data_t and return a type compatible with the input of the norm operator (usually a data_t).
     *
     * The convergence is assumed as soon as the norm between two subsequent states is less than a certain @f$\epsilon@f$, that is
     * the functor returns true if:
     * @f[\left\|v_t - v_{t-1}\right\| < \epsilon@f]
     *
     * @note When the convergence is reached, the previous states are not updated anymore (the calling algorithm is supposed to stop).
     */
    template <class data_t, class norm_t = norm_infinity>
    struct convergence_check
    {
      const double epsilon;

      //! Holds an instance of the norm used for checking the convergence.
      norm_t norm_comparison;

      data_t previous_state;

      //! Initialise the instance with the initial state of the vector.
      convergence_check(data_t const& current_state, double epsilon_ = 1E-5) : 
        epsilon(epsilon_), 
        previous_state(current_state)
      {}

      //! Returns true on convergence.
      bool operator()(data_t const& current_state)
      {
        bool ret = norm_comparison(current_state - previous_state) < epsilon;
        if(!ret)
        {
          previous_state = current_state;
        }
        return ret;
      }

    };


    // some issues with the random number generator
    const double fVeryBigButStillComputable = 1E10;
    const double fVerySmallButStillComputable = -1E10;


    /*!@brief Used to initialize the initial guess to some random data.
     *
     * @tparam data_t the type of the data returned. It should model a vector:
     *  - default, copy constructible and constructible with a size
     *  - has a member @c size returning a value of type @c data_t::Size_type
     *  - has a field @c data_t::Value_type
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



    namespace threading
    {

      //! Ensures the proper stop of the processing pool and the finalisation of all threads.
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


      //! Helper structure for managing additions of standard uBlas vectors.
      //! 
      //! This class is intended to be used with asynchronous_results_merger.
      template <class data_t>
      struct merger_addition
      {
        bool operator()(data_t &current_state, data_t const& update_value) const
        {
          current_state += update_value;
          return true;
        }
      };


      //! Helper structure for managing initialisations of standard uBlas vectors.
      //! 
      //! This class is intended to be used with asynchronous_results_merger.
      template <class data_t>
      struct initialisation_vector_specific_dimension
      {
        const size_t data_dimension;
        typedef typename data_t::value_type scalar_t;

        initialisation_vector_specific_dimension(size_t dimension) : data_dimension(dimension)
        {}

        bool operator()(data_t & current_state) const
        {
          current_state = boost::numeric::ublas::scalar_vector<scalar_t>(data_dimension, 0);
          return true;
        }
      };




      /*!@brief Merges the result of all workers and signals the results to the main thread.
       *
       * The purpose of this class is to add the computed accumulator of each thread to the final result
       * which contains the sum of all accumulators. 
       *
       */
      template <class result_type_, class merger_type, class init_result_type>
      struct asynchronous_results_merger : boost::noncopyable
      {
      public:
        typedef result_type_ result_type;

      protected:
        typedef boost::recursive_mutex mutex_t;
        typedef boost::lock_guard<mutex_t> lock_t;
        mutable mutex_t internal_mutex;        //!< mutex for thread safety
        //! Holds the current value of the merge.
        //! This variable is constantly updated as chunk processed finish. 
        result_type current_value;                  

        //! Holds the instance of the class responsible for merging new values (updates) to the
        //! current instance (current_value).
        merger_type merger_instance;

        //! Holds the instance of the class responsible for initialising the current value to
        //! an initial state (before any merge arrives).
        init_result_type initialisation_instance;

        //! Number of updates after the initialisation
        volatile int nb_updates;

        //! Thread synchronisation (event sent after an update, for counting).
        boost::condition_variable_any condition_;

      public:

        /*!Constructor
         *
         * @param initialisation_instance_ an instance of the class initialising the current state.
         */
        asynchronous_results_merger(init_result_type const &initialisation_instance_) : 
          initialisation_instance(initialisation_instance_)
        {}

        //! Initializes the internal states
        void init()
        {
          init_results();
          init_notifications();
        }

        //! Initialises the internal state of the accumulator
        void init_results()
        {
          initialisation_instance(current_value);
        }


        //! Initialises the number of notifications.
        //! Also called by init.
        void init_notifications()
        {
          nb_updates = 0;
        }

        /*! Receives the updated value of the vectors to accumulate from each worker.
         * 
         *  @note The call is thread safe.
         */
        void update(result_type const& updated_value)
        {
          lock_t guard(internal_mutex);
          merger_instance(current_value, updated_value);
        }



        /*! Function receiving the update notification.
         * 
         *  @note The call is thread safe.
         */
        void notify()
        {
          lock_t guard(internal_mutex);
          nb_updates ++;
          condition_.notify_one();
        }
     
        //! Returns once the number of updates reaches the number in argument.
        //!
        //!@warning if an inappropriate number is given, the method might never return.
        bool wait_notifications(size_t nb_notifications)
        {
          boost::unique_lock<mutex_t> lock(internal_mutex);
          while (nb_updates < nb_notifications)
          {
            // when entering wait, the lock is unlocked and made available to other threads.
            // when awakened, the lock is locked before wait returns. 
            condition_.wait(lock);
          }
        
          return true;
        
        }


        //! Returns the current merged results.
        //! @warning the call is not thread safe.
        result_type const& get_merged_result() const
        {
          return current_value;
        }

      };



    } // namespace threading
  } // namespace details
} // namespace robust_pca


#endif /* ROBUST_PCA_UTILITIES_HPP__*/ 
