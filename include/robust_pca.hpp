#ifndef ROBUST_PCA_HPP__
#define ROBUST_PCA_HPP__

//Released under the boost licence v1.0

/*!@file
 * Robust PCA functions, following the paper of Soren Hausberg
 * @todo add proper licence
 */

namespace robust_pca
{


  /*!@brief Robust PCA subspace algorithm
   *
   * This class implements the robust PCA using the Grassmanian averaging.
   * @todo add reference
   * @author Soren Hausberg, Raffi Enficiaud
   */
  template <class data_t>
  struct robust_pca_impl
  {
    robust_pca_impl(){}



    /*! Performs the computation of the current subspace on the elements given by the two iterators
     */
    template <class it_t>
    bool batch_process(it_t it, it_t const ite)
    {
      return true;
    }

    /*!@brief Performs the projection of the data on the space orthogonal to the current subspace
     *
     * After having found the current subspace, all the data should be projected onto the orthogonal subspace
     * of the current state.
     */
    template <class it_t>
    bool project_subspace(it_t it, it_t const ite) const
    {
      return true;
    }
  };

}

#endif /* ROBUST_PCA_HPP__ */