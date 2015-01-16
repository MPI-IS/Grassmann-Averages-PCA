// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

#ifndef GRASSMANN_AVERAGES_PCA_UTILITIES_HPP__
#define GRASSMANN_AVERAGES_PCA_UTILITIES_HPP__

/*!@file
 * Grassmann averages for robust PCA, companion functions.
 *
 * This file contains some utility function for multithreading, norm computation, convergence check
 * 
 */

#include <cmath>
#include <boost/numeric/conversion/bounds.hpp>


#include <boost/asio/io_service.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/recursive_mutex.hpp>


#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#include <numeric>
#include <set>

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif


// lock free queue, several producers, one consumer
//#include <boost/lockfree/queue.hpp>

namespace grassmann_averages_pca
{

  //! A callback class for monitoring the advance of the algorithm
  //!
  //! All calls are made in the main thread: there is no thread-safety issue.
  template <class data_t>
  struct grassmann_trivial_callback
  {

    //! Called to provide important messages/logs
    void log_error_message(const char* message) const
    {
      std::cout << message << std::endl;
    }

    //! This is called after centering the data in order to keep track 
    //! of the mean of the dataset
    void signal_mean(const data_t& mean) const
    {}

    //! Called after the computation of the PCA
    void signal_pca(const data_t& mean,
                    size_t current_eigenvector_dimension) const
    {}

    //! Called each time a new eigenvector is computed
    void signal_eigenvector(const data_t& current_eigenvector, 
                            size_t current_eigenvector_dimension) const
    {}

    //! Called at every step of the algorithm, at the end of the step
    void signal_intermediate_result(
      const data_t& current_eigenvector_state, 
      size_t current_eigenvector_dimension,
      size_t current_iteration_step) const
    {}

  };


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
    

    //! A trivial version of acos that avoids NaN
    template <class precision = double>
    struct safe_acos
    {
      typedef precision result_type;

      template <class T>
      result_type operator()(T const& v) const
      {
        using namespace std; // bringing acos

        if(v > 1) 
        {
          return 0;
        }
        else if(v < -1)
        {
          return precision(M_PI);
        }
        return acos(v);
      }
    };

    //! An faster approximate acos function with linear interpolation
    //!
    //! This implementation just interpolate the function acos on pivot points. The number of 
    //! pivot points is given at construction time
    //!
    //! @tparam precision precision of the internal representation of the values. 
    template <class precision = double>
    struct safe_acos_linear_interpolation
    {
      typedef precision result_type;

    private:
      std::vector<precision> v_table;
      precision step;
      precision step_inv;

    public:

      //!@brief Construct the internal interpolation table.
      //! 
      //! @param[in] nb_coefficients number of coefficients stored (which is also the size of the internal table minus 1).
      safe_acos_linear_interpolation(size_t nb_coefficients = 100) : 
         v_table(nb_coefficients + 1)
      {
        using namespace std; // bringing acos

        step = static_cast<precision>(2./nb_coefficients);
        step_inv = static_cast<precision>(nb_coefficients / 2.);

        for(size_t s = 0; s < nb_coefficients; s++)
        {
          v_table[s] = acos(-1 + s * step);
        }
        v_table[nb_coefficients] = 0;
      }

      //! Evaluation of the acos function
      template <class T>
      result_type operator()(T v) const
      {
        // this branch is killing everything
        if(v >= 1) 
        {
          return 0;
        }
        else if(v <= -1)
        {
          return precision(M_PI);
        }

        v = (v + 1) * step_inv;
        size_t index = static_cast<size_t>(v);
        v -= index;

        return (1-v) * v_table[index] + v * v_table[index + 1];

      }
    };


    /*!brief Handles the changes between two consecutive estimations of mus in order to avoid
     *       the computation of the inner product if possible.
     *
     */
    template <class data_t>
    struct mus_distance_optimisation_helper
    {

    private:
      //! Structure intended to store the mus and the number 
      //! of vectors referencing them. 
      struct mu_reference_count_data
      {
        size_t count;               //!< the number of vectors referencing the mu
        size_t iteration_index;     //!< the index of the iteration on which this mu appeared
        data_t mu;                  //!< the data
    
        mu_reference_count_data(const data_t& mu_) : mu(mu_) {}
      };

      typedef typename data_t::value_type scalar_t;

      // we need to keep the mus for computing the inner product
      typedef std::list<mu_reference_count_data> mu_container_t;
      mu_container_t mus;
  
      // This will map the index of the mu from the algorithm to the index of the list @c mus
      typedef std::map<size_t, size_t> map_index_to_list_index_t;
      map_index_to_list_index_t mapping_indices;
  
      // The matrix containing the angles between the vectors stored in the mu container. 
      // This is a symetric matrix and the indices follow the one of the list @c mus
      typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower> mu_angles_t;
      mu_angles_t mu_angles;
 

    public:

      mus_distance_optimisation_helper() : mus(), mapping_indices(), mu_angles()
      {}


      //! Clears the content 
      void clear()
      {
        mus.clear();
        mapping_indices.clear();
        mu_angles.clear();
      }
  
      //! Computes the angles between two vectors. 
      //!
      //! This function is meant to be used by general vectors (non normalized). If the vectors are 
      //! normalised already, consider using angle_normalized_vectors. 
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

        out /= sqrt(norm_l) * sqrt(norm_r);
        if(out > 1) 
        {
          out = 1;
        }
        else if(out < -1)
        {
          out = -1;
        }
        return acos(out);
      }
  
      //! Computes the angles between two vectors. 
      //!
      //! This version supposes the vectors are already normalized
      //! and the corresponding computations are simplified compared to @ref angle
      double angle_normalized_vectors(const data_t& left, const data_t& right) const
      {
        using namespace std; // bringing acos

        typedef typename data_t::value_type scalar_t;
        double out(0);
    
        scalar_t const *const left_p = &left.data()[0];
        scalar_t const *const right_p= &right.data()[0];
        for(size_t i = 0, j = left.size(); i < j; i++)
        {
          out += left_p[i] * right_p[i];
        }

        if(out > 1) 
        {
          out = 1;
        }
        else if(out < -1)
        {
          out = -1;
        }

        return acos(out);  
      }
  
      //! Removes any unnecessary element from the list of mus
      //! in order to lower the computation of the symmetric matrix to its bare minimum
      //! @return the number of pruned elements
      size_t prune()
      {
        std::set<size_t> index_to_remove;

        // identifies the rows/columns to be pruned
        {
          size_t index(0);
          for(typename mu_container_t::iterator it(mus.begin()), ite(mus.end()); it != ite; ++it, index++)
          {
            if(!it->count)
            {
              index_to_remove.insert(index);
              it = mus.erase(it);
              if(it == mus.end())
                break;
            }
          }
        }
    

        // reconstructs the matrix
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


        // reconstruct the mapping index
        mapping_indices.clear();
        {
          size_t index(0);
          for(typename mu_container_t::iterator it(mus.begin()), ite(mus.end()); it != ite; ++it, index++)
          {
            mapping_indices[it->iteration_index] = index;
          }
        }

        return index_to_remove.size();
      }


      //! Returns the number of vectors that are currently held by the instance
      size_t get_nb_mus() const
      {
        return mu_angles.size1();
      }

      //! Returns the angles matrix. 
      //!
      //! This should not be used for accessing the angles directly as the indices of the elements
      //! do not follow the indices of the algorithm. The mapping provided by mu_reference_count_data::iteration_index
      //! should be used jointly to the access to this matrix.
      mu_angles_t const& get_angle_matrix() const
      {
        return mu_angles;
      }

      //! Returns the list of elements that are stored into this container. 
      mu_container_t const& get_mus() const
      {
        return mus;
      }


      //! Gets the value of the angle between two iteration index vectors
      double get_angle_from_indices(size_t index1, size_t index2) const
      {
        assert(mapping_indices.count(index1) > 0);
        assert(mapping_indices.count(index2) > 0);
        return mu_angles(mapping_indices.at(index1), mapping_indices.at(index2));
      }
  
      //! Updates the count of a particular mu
      //!
      //! @param[in] iteration_index index manipulated by the algorithm, indicating the time of arrival of the mu that
      //!            is targetted by the update.
      //! @param[in] delta_count positive count means added reference, while negative count means
      //!            removed reference. 
      //! 
      //! @pre 
      //! This value of delta_count + current_count should be positive.
      void update_count(size_t iteration_index, int delta_count)
      {
        assert(mapping_indices.count(iteration_index) > 0);
        size_t index_in_list = mapping_indices[iteration_index];
    
        typename mu_container_t::iterator it(mus.begin());
        std::advance(it, index_in_list);
        assert(it->count + delta_count >= 0);
        it->count += delta_count;
      }


      //! Adds a mu to the managed list of mus
      //!
      //! The mu is stored in the first available place. 
      //!
      //! @param[in] new_mu new value of mu to be stored
      //! @param[in] iteration_index iteration on the algorithm side at which new_mu arrived. This is the value that should
      //!            be consecutevely used for accessing the angles. 
      //! @param[in] nb_references number of references to this mu
      //! 
      //! @note 
      //! It is possible to have a nb_references equal to 0 since this is what happens for a newly computed mu. 
      mu_reference_count_data const& add_mu(const data_t& new_mu, size_t iteration_index, size_t nb_references)
      {
    
        // finding the first available place to put this new mu:
        size_t index_in_list(0);
        typename mu_container_t::iterator it(mus.begin());
        typename mu_container_t::iterator ite(mus.end());
    
        for(; it != ite; ++it, index_in_list++)
        {
          if(!it->count)
          {
            break;
          }
        }
    
        if(it != ite)
        {
          it->mu = new_mu;
          it->count= nb_references;
          it->iteration_index = iteration_index;
        }
        else
        {
          // there is nothing to prune here anyway, so we add a new column/row to the 
          // matrix
      
          mus.push_back(mu_reference_count_data(new_mu));
          it = --mus.end();
          it->iteration_index   = iteration_index;
          index_in_list    = mus.size()-1;
          it->count = nb_references;      
      
          // adding a row / column with preservation of the content
          mu_angles.resize(mu_angles.size1() + 1, true);
          mu_angles(mu_angles.size1()-1, mu_angles.size1()-1) = 0;
            
        }
    
        mapping_indices[iteration_index] = index_in_list;
    
    
        size_t current_index(0);
        for(typename mu_container_t::const_iterator it2(mus.begin()); it2 != ite; ++it2, current_index++)
        {
          // computing the inner product and storing into the appropriate place of
          // the "distance" matrix
    
          if(it2 == it)
            continue;
       
          // not updating against the one that will be pruned
          if(it2->count == 0)
            continue;
    
          mu_angles(current_index, index_in_list) = angle(it->mu, it2->mu);
        }
    
        return *it;
      }

    };






    /*!@brief Gram Schmidt orthonormalisation of a collection of vectors.
     * @tparam it_t iterator on the collection of vectors. Should model a forward input iterator.
     * @tparam norm_t type of the norm operator.
     * 
     * @param it beginning of the collection of vectors
     * @param ite end of the collection of vectors
     * @param start first element of the collection to be orthonormalized. start should be inside the range given by it and ite. 
     * @param norm_op the norm used to normalise the vectors
     */
    template <class it_t, class norm_t>
    bool gram_schmidt_orthonormalisation(it_t it, it_t ite, it_t start, norm_t const &norm_op)
    {
      
      if(start == it)
      {
        *start *= typename it_t::value_type::value_type(1./norm_op(*start));
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
        current *= typename it_t::value_type::value_type(1./norm_op(current));
              
      }
      return true;
    }


    /*!@brief Computes the mean of a data set after having removed the lower and upper k first elements.
     *
     * This function computes @f[\sum_{k \leq i < N-k} p_{o(i)}@f] where 
     * - @f$p@f$ is the data set of size @f$N@f$, and @f$p_i@f$ is its ith element
     * - @f$o@f$ is a function ordering the data set: @f$\forall i, p_{o(i)} \leq p_{o(i+1)}, 0 \leq i < N @f$

     * @tparam T type of the data set. All internal accumulations will be performed with this type. T should not be const as
     *         the data set will be modified in place.
     *
     * @param p_data the data set composed of nb_total_elements of type T. 
     * @param nb_total_elements number of elements of the data set
     * @param k_first_last number of elements to remove from the lower and upper distributions.
     *
     * @pre @f$ \text{k_first_last} \leq \frac{\text{nb_total_elements}}{2}@f$
     */
    template <class T>
    T compute_mean_within_bounds(T *p_data, size_t nb_total_elements, size_t k_first_last)
    {
      if(k_first_last < nb_total_elements / 2)
      {
        std::nth_element(p_data, p_data + k_first_last, p_data + nb_total_elements);
        std::nth_element(p_data + k_first_last+1, p_data + nb_total_elements - k_first_last-1, p_data + nb_total_elements);
        T acc = std::accumulate(p_data + k_first_last, p_data + nb_total_elements - k_first_last, T(0));
          
        return acc / (nb_total_elements - 2*k_first_last);
      }
      else
      {
        assert(k_first_last == nb_total_elements / 2);
        std::nth_element(p_data, p_data + k_first_last, p_data + nb_total_elements);
          
        if(nb_total_elements & 1)
        {
          return *(p_data + k_first_last);
        }
        else
        {
          return (*(p_data + k_first_last) + *std::max_element(p_data, p_data + k_first_last)) / 2;
        }
      }
    }




    /*!@brief Checks the convergence of a sequence.
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
     * @note Once the convergence is reached, the internal states are not updated anymore (the calling algorithm is supposed to stop).
     */
    template <class data_t, class norm_t = norm_infinity>
    struct convergence_check
    {
      //! The amount of change below which the sequence is considered as having reached a steady point.
      const double epsilon;

      //! Holds an instance of the norm used for checking the convergence.
      norm_t norm_comparison;

      //! The previous value
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


    //! @namespace
    namespace threading
    {


      //! Ensures the proper stop of the processing pool and the finalisation of all threads.
      struct safe_stop
      {
      private:
        boost::asio::io_service& io_service;
        boost::thread_group& thread_group;

      public:
        safe_stop(boost::asio::io_service& ios, boost::thread_group& tg) : io_service(ios), thread_group(tg)
        {}

        ~safe_stop()
        {
          io_service.stop();
          thread_group.join_all();
        }
      };


      //! @brief Helper structure for managing additions on standard uBlas vectors.
      //! 
      //! This class is intended to be used with asynchronous_results_merger. It just adds an update to the current state.
      //! @tparam data_t type of the vectors. It is supposed that data_t implements in-place addition (@c data_t::operator+=).
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
      //! @tparam data_t type of the vectors
      //! @note This implementation supposes that the type is compatible with boost::numeric::ublas::vector
      template <class data_t>
      struct initialisation_vector_specific_dimension
      {
      private:
        const size_t data_dimension;                    //!< Dimension of the vectors
        typedef typename data_t::value_type scalar_t;   //<! Scalar type

      public:
        //! Initialise the instance with the dimension of the data. 
        //! The dimension is fixed. 
        initialisation_vector_specific_dimension(size_t dimension) : data_dimension(dimension)
        {}

        //! Initialise the current state a null (0) vector of the dimension guiven at construction.
        bool operator()(data_t & current_state) const
        {
          current_state = boost::numeric::ublas::scalar_vector<scalar_t>(data_dimension, 0);
          return true;
        }
      };




      /*!@brief Merges the result of all workers and signals the results to the main thread.
       *
       * The purpose of this class is to gather the computation results coming from several threads into one unique result seen by the main calling thread. 
       * Each thread computes a partial update of the final result. These partial update are signalled to this instance via @c asynchronous_results_merger::update (thread safe). 
       * These updates are gathered/merged to the final result through the "merger" instance (of type @c merger_type) in a thread safe manner.
       * The number of updates is also signalled to the main thread via a call to @c asynchronous_results_merger::notify. The main thread supposes the computation over/in sync if it received
       * an amount of notification through the @c asynchronous_results_merger::wait function.
       *
       * @tparam result_type_ the type of the final result.
       * @tparam merger_type the type of the merger. The merger should be a callable with two arguments: result_type_ and update_element_
       * @tparam init_result_type the type of the initialiser. The initialiser should be a callable with one argument of type result_type_.
       * @tparam update_element_ the type of the update. These updates are provided by the several workers to this merger. 
       *
       * @note This implementation supposes that the pointers to the update elements remain after the call to @c asynchronous_results_merger::notify. This is because
       * the implementation tries to avoid any "long" or time consuming lock. If the merge cannot be performed in the asynchronous_results_merger::update call itself,
       * then the update element is queued and the merge is performed in the main calling thread (the wait function). 
       */
      template <class result_type_, class merger_type, class init_result_type, class update_element_ = result_type_>
      struct asynchronous_results_merger : boost::noncopyable
      {
      public:

        typedef result_type_ result_type;         //!< The type returned by asynchronous_results_merger::get_merged_result
        typedef update_element_ update_element;   //!< The type used for the updates.

      protected:
        typedef boost::recursive_mutex mutex_t;     //!< Type of the mutex. This one is re-entrant/recursive in order to allow the same thread locking it several times.
        typedef boost::lock_guard<mutex_t> lock_t;  //!< Exclusive lock

        //! Mutex for critical sections. 
        //!@note This mutex is re-entrant.
        mutable mutex_t internal_mutex;

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
          initialisation_instance(initialisation_instance_),
          nb_updates(0)
        {}

        //! Initializes the internal states
        virtual void init()
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

        /*! Receives the update element from each worker.
         * 
         * The update element is passed to the merger in order to create an updated value of the internal result.
         * @note The call is thread safe.
         */
        void update(update_element const* updated_value)
        {
          boost::unique_lock<mutex_t> lock(internal_mutex);
          merger_instance(current_value, *updated_value);
        }



        /*! Function receiving the update notification.
         * 
         *  @note The call is thread safe.
         */
        void notify()
        {
          lock_t guard(internal_mutex);
          ++nb_updates;

          //
          // The condition functions are not async-signal safe, and should not be called from a signal handler. 
          // In particular, calling pthread_cond_signal or pthread_cond_broadcast from a signal handler may
          // deadlock the calling thread.
          // 
          // Raffi: the notification outside the lock above causes a deadlock, apparently it should be protected
          // from simultaneous access. 
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

            //assert(nb_updates);           // cannot be awakened if there is no update
            assert(lock.owns_lock());
          }

          assert(lock.owns_lock());
          return true;
        }


        //! Returns the current merged results.
        //! @warning the call is not thread safe (intended to be called once the wait_notifications returned and no
        //! other thread is working). 
        result_type const& get_merged_result() const
        {
          return current_value;
        }

        //! Returns the current merged results.
        //! @warning the call is not thread safe (intended to be called once the wait_notifications returned and no
        //! other thread is working). 
        result_type & get_merged_result()
        {
          return current_value;
        }
      };



    } // namespace threading
  } // namespace details
} // namespace grassmann_averages_pca


#endif /* GRASSMANN_AVERAGES_PCA_UTILITIES_HPP__*/ 
