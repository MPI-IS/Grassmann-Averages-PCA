// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

#ifndef BOOST_UBLAS_MATLAB_HELPER_GRASSMANN_AVERAGES_PCA_HPP__
#define BOOST_UBLAS_MATLAB_HELPER_GRASSMANN_AVERAGES_PCA_HPP__


#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>


/*!@file
 * This file contains some helper function for wrapping matlab arrays into C++. One of the achievements
 * is the implementation of a specific storage for uBlas wrapping existing array, which avoids the need
 * to copy data, and behaves well with uBlas. 
 */


namespace grassmann_averages_pca
{
  namespace details
  {
    namespace ublas_helpers
    {


      namespace ub = boost::numeric::ublas;


      /*!@brief uBlas allocator for external storage.
       *
       * Adapted from Gunter Winkler tricks http://www.guwi17.de/ublas/index.html. This class should follow the storage concept
       * of uBlas. However, the size of the container cannot change (same memory area), some functions are omitted and should generate a compile-time error.
       * 
       * When used from Matlab, the matrices are column_major. The storage is as-is there, so the matrix should be declared column major.
       * @author Raffi Enficiaud
       *
       * @todo may be replaced by array_adaptor implementation, which is not properly documented in the boost.uBlas.
       */
      template <class T>
      class external_storage_adaptor : public ub::storage_array< external_storage_adaptor<T> >
      {

        typedef external_storage_adaptor<T> this_type;

      public:
        typedef size_t size_type;
        typedef ptrdiff_t difference_type;
        typedef T value_type;

        typedef typename boost::add_reference<T>::type reference;
        typedef typename boost::add_pointer<T>::type pointer;
        typedef typename boost::add_pointer<typename boost::add_const<T>::type>::type const_pointer;


      public:

        //! Default constructor
        external_storage_adaptor() :
          size_(0), 
          data_(0)
        {}


        external_storage_adaptor(size_type size, pointer data) :
          size_(size), 
          data_(data)
        {}


        /*! Constructor from another instance
         *
         * The two instances will share the data
         */
        external_storage_adaptor(const this_type& rhs) : 
          size_(rhs.size_), data_(rhs.data_)
        {}

        ~external_storage_adaptor() {}


        // Resizing
        void resize(size_type size)
        {
          // the size cannot change, otherwise there is a problem with the external storage.
          // it can only be set to 0
          assert(size == 0 || size_ == 0 || size == size_);
          size_ = size;
        }

        void resize(size_type size, pointer data)
        {
          assert(size == 0 || size_ == 0 || size == size_);
          size_ = size;

          assert(data == 0 || data_ == 0 || data == data_);
          data_ = data;
        }

        // Random Access Container
        size_type max_size() const
        {
          return std::numeric_limits<size_type>::max();
        }
        
        //! Returns true if the storage is empty, false otherwise.
        bool empty() const
        {
          return size_ == 0;
        }

        //! Returns the current size contained.
        size_type size() const
        {
          return size_;
        }

        //!@name Random access container concept
        //!@{

        //! Implements the random access container.
        reference operator [] (size_type i) const
        {
          assert(i < size_);
          return data_[i];
        }

        // Iterators simply are pointers. There is no need to
        // provide an iterator for the const version, as the normal pointer is convertible 
        // to the const pointer. There is no ambiguity with the call neither, since the 
        // call should not change the storage itself.
        typedef pointer iterator;
        typedef const_pointer const_iterator;

        iterator begin() const
        {
          return data_;
        }
        iterator end() const
        {
          return data_ + size_;
        }

        // Reverse iterators
        typedef std::reverse_iterator<iterator> reverse_iterator;
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

        reverse_iterator rbegin() const
        {
          return reverse_iterator(end());
        }

        reverse_iterator rend() const
        {
          return reverse_iterator(begin());
        }

        //! @} 
        // Random access container concept.

      private:
        size_type size_;
        pointer data_;
      };




    }
  }
}


namespace boost { namespace numeric { namespace ublas {
  /*! Specializing the vector temporary type policy of uBlas.
   *
   * Sticking to the standard definition of this class would create temporary
   * vectors with the same allocation object as the matrix to which it is applied. 
   * However, external_storage_adaptor does not allow construction with size only. A 
   * way to compile 
   * @code
   * matrix_row< matrix<double, row_major, external_storage_adaptor> >(some_matrix, some_row) = zero_vector<double>(nb_columns)
   * @endcode
   * is precisely to avoid the creation of temporaries with external_storage_adaptor allocator.
   */
  template <class T, class L>
  struct vector_temporary_traits<matrix<T, L, grassmann_averages_pca::details::ublas_helpers::external_storage_adaptor<T> > >
  {
    typedef vector<T> type;
  };
}}}


#endif /* BOOST_UBLAS_MATLAB_HELPER_GRASSMANN_AVERAGES_PCA_HPP__ */
