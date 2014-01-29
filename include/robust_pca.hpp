#ifndef ROBUST_PCA_HPP__
#define ROBUST_PCA_HPP__

//Released under the boost licence v1.0

/*!@file
 * Robust PCA functions, following the paper of Soren Hausberg
 * @todo add proper licence

 *
 * @note These implementations assume the existence of boost somewhere. Currently, it is just used for 
 * computing the norms of differences and to generate some random data (which is now TR1 and C++0X).
 */


#include <boost/numeric/conversion/bounds.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>



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
    convergence_check(double epsilon_) : 
      epsilon(epsilon_), 
      previous_state(std::numeric_limits<double>::max())
    {}

    bool operator()(data_t const& current_state)
    {
      return norm_comparaison(current_state - previous_state) < epsilon;
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
   * @author Soren Hausberg, Raffi Enficiaud
   */
  template <class data_t, class norm_2_t = norm2>
  struct robust_pca_impl
  {
  private:
    random_data_generator<data_t> random_init_op;
    norm_2_t norm_op;

  public:
    robust_pca_impl() : random_init_op(fVerySmallButStillComputable, fVeryBigButStillComputable)
    {}



    /*! Performs the computation of the current subspace on the elements given by the two iterators.
     *  @tparam it_t an input iterator to vectors. Each element pointed by the underlying iterator should be iterable and
     *   should provide a vector.
     *  @tparam it_norm_t an output iterator on weights/norms of the vectors. The output elements should be numerical (norm output)
     *
     * @param[in] initial_guess if provided, the initial vector will be initialized to this value. 
     *
     * @pre 
     * - @c !(it >= ite)
     * - all the vectors given by the iterators pair should be of the same size (no check is performed).
     */
    template <class it_t, class it_norm_t>
    bool batch_process(it_t it, it_t const ite, it_norm_t it_norm_out, data_t const * initial_guess = 0)
    {
      // add some log information
      if(it >= ite)
      {
        return false;
      }

      it_norm_t it_norm_out_copy(it_norm_out);
      for(it_t it_copy(it); it_copy != ite; ++it_copy, ++it_norm_out_copy)
      {
        *it_norm_out_copy = norm_op(*it_copy);
      }


      // the first element is used for the init guess because for dynamic std::vector like element, the size is needed.
      data_t mu(initial_guess != 0 ? *initial_guess : random_init_op(*it));

      // normalizing
      typename norm_2_t::result_type norm_mu(norm_op(mu));
      mu /= norm_mu;
      

      return true;
    }

    /*!@brief Performs the projection of the data on the space orthogonal to the current subspace
     *
     * After having found the current subspace, all the data should be projected onto the orthogonal subspace
     * of the current state.
     */
    template <class it_t>
    bool projection_onto_subspace(it_t it, it_t const ite) const
    {
      return true;
    }
  };

}

#endif /* ROBUST_PCA_HPP__ */
