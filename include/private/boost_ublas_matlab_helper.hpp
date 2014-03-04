#ifndef BOOST_UBLAS_MATLAB_HELPER_ROBUST_PCA_HPP__
#define BOOST_UBLAS_MATLAB_HELPER_ROBUST_PCA_HPP__

#include "mex.h"

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>



namespace robust_pca
{
  namespace ublas_matlab_helper
  {


    namespace ub = boost::numeric::ublas;


    /*!@brief uBlas allocator for external read-only storage
     *
     * Adapted from Gunter Winkler tricks http://www.guwi17.de/ublas/index.html. This class should follow the storage concept
     * of uBlas. However, since the data is read-only, some functions are omitted and should generate a compile-time error.
     * 
     * When used from Matlab, the matrices are column_major. The storage is as-is there, so the matrix should be declared column major.
     * @author Raffi Enficiaud
     */
    template<class T>
    class readonly_array_adaptor : public ub::storage_array< readonly_array_adaptor<T> >
    {

      typedef readonly_array_adaptor<T> this_type;

    public:
      typedef size_t size_type;
      typedef ptrdiff_t difference_type;
      typedef T value_type;

      typedef const T &const_reference;
      typedef const T *const_pointer;


    public:

      //! Default constructor
      readonly_array_adaptor() :
        size_(0), 
        data_(0)
      {}


      readonly_array_adaptor(size_type size, const_pointer data) :
        size_(size), 
        data_(data) {
      }


      /*! Constructor from another instance
       *
       * The two instances will share the data
       */
      readonly_array_adaptor(const this_type& rhs) : 
        size_(rhs.size_), data_(rhs.data_)
      {}

      ~readonly_array_adaptor() {}


      // Resizing
      void resize(size_type size)
      {
        size_ = size;
      }

      void resize(size_type size, const_pointer data)
      {
        size_ = size;
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
      const_reference operator [] (size_type i) const
      {
        assert(i < size_);
        return data_[i];
      }

      // Iterators simply are pointers.
      typedef const_pointer const_iterator;

      const_iterator begin() const
      {
        return data_;
      }
      const_iterator end() const
      {
        return data_ + size_;
      }

      // this typedef is used by vector and matrix classes
      typedef const_pointer iterator;

      // Reverse iterators
      typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
      typedef std::reverse_iterator<iterator> reverse_iterator;

      const_reverse_iterator rbegin() const
      {
        return const_reverse_iterator(end());
      }

      const_reverse_iterator rend() const
      {
        return const_reverse_iterator(begin());
      }

      //! @} 
      // Random access container concept.

    private:
      size_type size_;
      const_pointer data_;
    };




  }



}


#endif /* BOOST_UBLAS_MATLAB_HELPER_ROBUST_PCA_HPP__ */
