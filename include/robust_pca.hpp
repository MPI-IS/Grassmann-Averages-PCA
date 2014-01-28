#ifndef ROBUST_PCA_HPP__
#define ROBUST_PCA_HPP__

//Released under the boost licence v1.0

/*!@file
 * Robust PCA functions, following the paper of Soren Hausberg
 * @todo add proper licence

 *
 * @note These implementations assume the existence of boost somewhere. Currently, it is just used for 
 * computing the norms of differences.
 */



#include <boost/numeric/ublas/vector_expression.hpp>

namespace robust_pca
{
 
  //! Wrapper object for infinity/max norm
  struct norm_inf
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
        acc += (*it) * (*it);
      }
      return acc;
    }
  };


  /* @brief Checks the convergence of a numerical scheme.
   *
   * @tparam data_t should be default constructible and constructible with one parameter.
   *
   * The convergence is assumed when the norm between two subsequent states is less than a certain @f$epsilon@f$.
   */
  template <class data_t, class norm_t = norm_inf>
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


  /*!@brief Robust PCA subspace algorithm
   *
   * This class implements the robust PCA using the Grassmanian averaging. This is the naive implementation which is
   * suitable for small datasets. 
   * @todo add reference
   * @author Soren Hausberg, Raffi Enficiaud
   */
  template <class data_t>
  struct robust_pca_impl
  {
    robust_pca_impl(){}



    /*! Performs the computation of the current subspace on the elements given by the two iterators.
     *  @tparam it_t an input iterator to vectors. Each element pointed by the underlying iterator should be iterable and
     *   should provide a vector.
     *  @tparam it_norm_t an output iterator on weights/norms of the vectors. The output elements should be numerical (norm output)
     */
    template <class it_t, class it_norm_t>
    bool batch_process(it_t it, it_t const ite, it_norm_t it_norm_out)
    {
      it_norm_t it_norm_out_copy(it_norm_out);
      for(it_t it_copy(it); it_copy != ite; ++it_copy, ++it_norm_out_copy)
      {
        *it_norm_out_copy = boost::numeric::ublas::norm_2(*it_copy);
      }
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
