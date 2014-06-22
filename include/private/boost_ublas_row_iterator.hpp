// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

#ifndef BOOST_UBLAS_MATRIX_HELPER_GRASSMANN_AVERAGES_PCA_HPP__
#define BOOST_UBLAS_MATRIX_HELPER_GRASSMANN_AVERAGES_PCA_HPP__

/*!@file
 * Contains utility classes around uBlas matrices(mainly an iterator on matrix rows).
 */

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/iterator/iterator_adaptor.hpp>



namespace grassmann_averages_pca
{
  namespace details
  {
    namespace ublas_helpers
    {


      /*!@brief Iterator on rows of a matrix.
       *
       * This iterator is an adaptor that iterates over the rows of an ublas matrix. The returned element
       * is a matrix proxy that provides an ublas vector semantic. 
       * @author Raffi Enficiaud
       */
      template <class matrix_t>
      class row_iter : 
        public boost::iterator_facade<
          row_iter<matrix_t>                            // Derived
        , boost::numeric::ublas::matrix_row<matrix_t>   // Value
        , std::random_access_iterator_tag               // CategoryOrTraversal
        , boost::numeric::ublas::matrix_row<matrix_t>   // reference
        >
      {
      private:
        typedef row_iter<matrix_t> this_type;

        // this is a nice technic for SFNIAE, taken from the examples of boost.iterator.
        struct enabler {};

        size_t index;
        matrix_t *matrix;

        typedef boost::numeric::ublas::matrix_row<matrix_t> return_t;

      public:

        //! Default constructor
        row_iter() : index(std::numeric_limits<size_t>::max()), matrix(0)
        {}

        /*! Constructs an iterator on the nth rows of the given matrix
         *
         * @param[in] matrix_ the matrix from which the rows will be extracted
         * @param[in] index_ the index of the row of the matrix on which the iterator starts.
         *
         * @pre the provided index is lower than the total number of rows of the matrix.
         */
        row_iter(matrix_t &matrix_, size_t index_) : index(index_), matrix(&matrix_)
        {
          // below the <= is correct since it should be possible to give an index one passed the end
          assert(index_ <= matrix->size1());
        }


        /*! Constructor from an iterator 
         *
         * This constructor is selected when the row_iter is instanciated from a row_iter on another type, and
         * type are convertible from the other matrix type to the current type (example non-const to const). Otherwise
         * the SFNIAE does reveal this constructor in the set of available constructor for this class.
         */
        template <class other_matrix_t>
        row_iter(
          row_iter<other_matrix_t> const& other, 
          typename boost::enable_if<
            boost::is_convertible<typename other_matrix_t::iterator1, typename matrix_t::iterator1>, 
            enabler>::type = enabler()) 
          : 
          index(other.index), matrix(other.matrix)
        {}

      private:
        friend class boost::iterator_core_access;

        //!@name boost::iterator_facade interface
        //!@{

        void increment()
        {
          assert(matrix);
          assert(index < matrix->size1());
          index++;
        }

        bool equal(this_type const& other) const
        {
          assert(matrix == other.matrix);
          return this->index == other.index;
        }

        return_t dereference() const
        {
          assert(matrix);
          return return_t(*matrix, index);
        }

        typename this_type::difference_type distance_to(this_type const& r) const
        {
          assert((matrix != 0) && (r.matrix == matrix));
          return typename this_type::difference_type(r.index) - typename this_type::difference_type(index); // sign promotion
        }

        void advance(typename this_type::difference_type n)
        {
          if(n < 0)
          {
            assert((-n) <= static_cast<typename this_type::difference_type>(index));
            index -= n;
          }
          else
          {
            assert(n + index <= matrix->size1());
            index += n;
          }
          
        }

        //!@}

      };

    }
  }
}

#endif
