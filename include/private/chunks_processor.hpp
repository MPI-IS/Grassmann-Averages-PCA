// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

#ifndef GRASSMANN_AVERAGES_PCA_CHUNK_PROCESSORS_HPP__
#define GRASSMANN_AVERAGES_PCA_CHUNK_PROCESSORS_HPP__

/*!@file
 * Chunks of data multithreaded processing.
 *
 */


#include <cmath>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <include/private/utilities.hpp>


namespace grassmann_averages_pca
{

    //!@internal
    //!@brief Contains the logic for processing part of the accumulator
    template <class data_t>
    struct asynchronous_chunks_processor
    {
    private:
      //! Type of the elements contained in the vectors
      typedef typename data_t::value_type scalar_t;

      //! Number of vectors contained in this chunk.
      size_t nb_elements;

      //! Dimension of the vectors.
      size_t data_dimension;

      //! Internal accumulator.
      //! Should live beyond the scope of update and init, as required by the merger.
      data_t accumulator;

      //! Signs stored for decreasing the number of updates
      std::vector<bool> v_signs;

      //! Index to which a vector is related regarding the computation of the inner product
      std::vector<size_t> v_mu_reference_indices;

      //! Vector of inner products
      //! This vector stores at index k the inner product of the at k and the corresponding 
      //! reference mu (whose index in the algorithm iteration is given by @c v_mu_reference_indices[k]).
      std::vector<scalar_t> v_angles_value;
      
      // this is to send an update of the value of mu to one listener
      // the connexion should be managed externally
      typedef boost::function<void (data_t const*)> connector_accumulator_t;
      connector_accumulator_t signal_acc;

      typedef boost::function<void ()> connector_counter_t;
      connector_counter_t signal_counter;

      // this is to send the delta of vectors referencing a specific index of mus
      typedef std::map<size_t, int> map_reference_delta_t;
      typedef boost::function<void (map_reference_delta_t const*)> connector_reference_delta_update_t;
      connector_reference_delta_update_t signal_reference_count;

      //! The matrix containing a copy of the data.
      //! The vectors are stored per row in this matrix.
      scalar_t *p_c_matrix;
      
      //! Padding for one line of the matrix
      size_t data_padding;

      //! This is the class that will check if the computation of the inner product is needed
      typedef details::mus_distance_optimisation_helper<data_t> inner_product_optimisation_t;
      inner_product_optimisation_t &inner_product_optimisation;

      //! "Optimized" inner product.
      //! This one has the particularity to be more cache/memory bandwidth friendly. More efficient
      //! implementations may be used, but it turned out that the memory bandwidth is saturated when
      //! many threads are processing different data.
      scalar_t inner_product(scalar_t const* p_mu, scalar_t const* current_line) const
      {
        scalar_t const * const current_line_end = current_line + data_dimension;

        const int _64_elements = static_cast<int>(data_dimension >> 6);
        scalar_t acc(0);

        for(int j = 0; j < _64_elements; j++, current_line += 64, p_mu += 64)
        {
          for(int i = 0; i < 64; i++)
          {
            acc += current_line[i] * p_mu[i];
          }
        }
        for(; current_line < current_line_end; current_line++, p_mu++)
        {
          acc += (*current_line) * (*p_mu);
        }
        return acc;
      }

      scalar_t normed_inner_product(scalar_t const* p_mu, scalar_t const* current_line) const
      {
        scalar_t const * const current_line_end = current_line + data_dimension;

        const int _64_elements = static_cast<int>(data_dimension >> 6);
        scalar_t acc(0);
        scalar_t acc_vector(0);

        for(int j = 0; j < _64_elements; j++, current_line += 64, p_mu += 64)
        {
          for(int i = 0; i < 64; i++)
          {
            acc += current_line[i] * p_mu[i];
            acc_vector += current_line[i] * current_line[i];
          }
        }
        for(; current_line < current_line_end; current_line++, p_mu++)
        {
          acc += (*current_line) * (*p_mu);
          acc_vector += (*current_line) * (*current_line);
        }
        return acc / sqrt(acc_vector);
      }


      //! "Optimized" inner product
      scalar_t inner_product(scalar_t const* p_mu, size_t element_index) const
      {
        return inner_product(p_mu, p_c_matrix + element_index * data_padding);
      }

      //! Returns the inner product divided by the l2 norm of the referenced vector
      scalar_t normed_inner_product(scalar_t const* p_mu, size_t element_index) const
      {
        return normed_inner_product(p_mu, p_c_matrix + element_index * data_padding);
      }



    public:
      asynchronous_chunks_processor(inner_product_optimisation_t &inner_product_optimisation_) : 
        nb_elements(0), 
        data_dimension(0), 
        p_c_matrix(0), 
        data_padding(0),
        inner_product_optimisation(inner_product_optimisation_)
      {
      }
      
      ~asynchronous_chunks_processor()
      {
        delete [] p_c_matrix;
      }


      //! Sets the data range
      template <class container_iterator_t>
      void set_data_range(container_iterator_t const &b, container_iterator_t const& e)
      {
        nb_elements = std::distance(b, e);
        assert(nb_elements > 0);
        v_signs.resize(nb_elements);
        v_mu_reference_indices.resize(nb_elements);
        v_angles_value.resize(nb_elements);

        // aligning on 32 bytes = 1 << 5
        data_padding = (data_dimension*sizeof(scalar_t) + (1<<5) - 1) & (~((1<<5)-1));
        data_padding /= sizeof(scalar_t);
        
        delete [] p_c_matrix;
        p_c_matrix = new scalar_t[data_padding*nb_elements];
        
        container_iterator_t bb(b);

        scalar_t *current_line = p_c_matrix;
        for(int line = 0; line < nb_elements; line ++, current_line += data_padding, ++bb)
        {
          for(int column = 0; column < data_dimension; column++)
          {
            current_line[column] = (*bb)(column);
          }
        }
        
        signal_counter();
      }

      //! Sets the dimension of each vectors
      //! @pre data_dimensions_ strictly positive
      void set_data_dimensions(size_t data_dimensions_)
      {
        data_dimension = data_dimensions_;
        assert(data_dimension > 0);
      }

      //! Returns the callback object that will be called to signal an updated accumulator.
      connector_accumulator_t& connector_accumulator()
      {
        return signal_acc;
      }

      //! Returns the callback object that will be called to signal the end of the current computation.
      connector_counter_t& connector_counter()
      {
        return signal_counter;
      }

      //! Returns the callback object that will be called to signal an updated reference count map (deltas).
      connector_reference_delta_update_t& connector_reference_count()
      {
        return signal_reference_count;
      }

      //! Centering the data in case it was not possible to do it beforehand
      void data_centering_first_phase(size_t full_dataset_size)
      {
        scalar_t const * current_line = p_c_matrix;
        accumulator = data_t(data_dimension, 0);
        scalar_t * const p_acc_begin = &accumulator(0);
        scalar_t const * const p_acc_end = p_acc_begin + data_dimension;
        
        
        for(size_t current_element = 0; 
            current_element < nb_elements; 
            current_element++, current_line+= data_padding - data_dimension)
        {
          scalar_t * p_acc = p_acc_begin;
          for(; p_acc < p_acc_end; p_acc++, current_line++)
          {
            *p_acc += *current_line;
          }
        }


        for(scalar_t * p_acc = p_acc_begin; p_acc < p_acc_end; p_acc++)
        {
          *p_acc /= full_dataset_size;
        }


        // posts the new value to the listeners for the current dimension
        signal_acc(&accumulator);
        signal_counter();

      }

      //! Project the data onto the orthogonal subspace of the provided vector
      void data_centering_second_phase(data_t const &mean_value)
      {
        scalar_t * current_line = p_c_matrix;
        scalar_t const * const p_mean_begin = &mean_value(0);
        scalar_t const * const p_mean_end = p_mean_begin + data_dimension;
        
        
        for(size_t current_element = 0; 
            current_element < nb_elements; 
            current_element++, current_line+= data_padding - data_dimension)
        {
          const scalar_t * p_mean = p_mean_begin;
          for(; p_mean < p_mean_end; p_mean++, current_line++)
          {
            *current_line -= *p_mean;
          }
        }

        // posts the new value to the listeners for the current dimension
        signal_counter();

      }


      //! PCA steps
      void pca_accumulation(data_t const &mu)
      {
        accumulator = data_t(data_dimension, 0);
        scalar_t const * const p_mu = &mu.data()[0];
        scalar_t * const p_acc = &accumulator.data()[0];
                
        scalar_t const * current_line = p_c_matrix;

        for(size_t s = 0; s < nb_elements; s++, current_line += data_padding)
        {
          const scalar_t inner_prod = inner_product(p_mu, current_line);
          for(size_t d = 0; d < data_dimension; d++)
          {
            p_acc[d] += inner_prod * current_line[d];
          }
        }

        // posts the new value to the listeners
        signal_acc(&accumulator);
        signal_counter();
      }



      //! Initialises the accumulator and the signs vector from the first mu
      void initial_accumulation(data_t const &mu)
      {
        details::safe_acos acos_function_object;

        accumulator = data_t(data_dimension, 0);
        v_mu_reference_indices.assign(v_mu_reference_indices.size(), 0); // all vectors will bind to 0, which is the current index

        std::vector<bool>::iterator itb(v_signs.begin());
        std::vector<scalar_t>::iterator it_angles(v_angles_value.begin());

        scalar_t const * const p_mu = &mu.data()[0];
        scalar_t * const p_acc = &accumulator.data()[0];

        // first iteration, we store the signs
        for(size_t s = 0; s < nb_elements; ++itb, ++it_angles, s++)
        {
          scalar_t const normed_inner = normed_inner_product(p_mu, s);
          bool const sign = normed_inner >= 0;

          *itb = sign;
          *it_angles = abs(acos_function_object(normed_inner)); // vectors are normalized, so the angle is the acos of the inner product

          scalar_t *p(p_acc);
          scalar_t const * current_line = p_c_matrix + s * data_padding;
          scalar_t const * const current_line_end = current_line + data_dimension;
          if(sign)
          {
            for(; current_line < current_line_end; ++current_line, ++p)
            {
              *p += *current_line;              
            }
          }
          else 
          {
            for(; current_line < current_line_end; ++current_line, ++p)
            {
              *p -= *current_line;              
            }
          }
        }


        // posts the new value to the listeners
        signal_acc(&accumulator);

        // every vector refers to the same mu
        map_reference_delta_t map_update;
        map_update[0] = static_cast<int>(nb_elements);
        signal_reference_count(&map_update);

        signal_counter();
      }


      //! Update the accumulator and the signs vector from an upated mu
      void update_accumulation(data_t const& mu, size_t current_index)
      {
        details::safe_acos acos_function_object;


        map_reference_delta_t map_update;
        accumulator = data_t(data_dimension, 0);
        size_t count_current_index(0); // should be positive only

        bool update = false;

        scalar_t const * const p_mu = &mu.data()[0];
        scalar_t * const p_acc = &accumulator.data()[0];
        scalar_t const * current_line = p_c_matrix;


        std::vector<bool>::iterator itb(v_signs.begin());
        std::vector<size_t>::iterator it_indices(v_mu_reference_indices.begin());
        std::vector<scalar_t>::iterator it_angles(v_angles_value.begin());

        for(size_t s = 0; 
            s < nb_elements; 
            ++itb, ++it_indices, ++it_angles, s++, current_line += data_padding)
        {
          const size_t previous_index = *it_indices;

          double previous_angle = inner_product_optimisation.get_angle_from_indices(previous_index, current_index);
          if(*it_angles + abs(previous_angle) < M_PI / 2)
          {
            // in this case, the scalar product does not need to be recomputed, and the sign does not change
            // against the previous one
            assert(*itb == (inner_product(p_mu, current_line) >= 0));
            continue;
          }

          scalar_t const normed_inner = normed_inner_product(p_mu, current_line);

          map_update[previous_index]--;
          count_current_index++;
          *it_indices = current_index;
          *it_angles = abs(acos_function_object(normed_inner));

          bool sign = normed_inner >= 0;
          if(sign != *itb)
          {
            update = true;

            // update the value of the accumulator according to sign change
            *itb = sign;
            
            scalar_t *p(p_acc);
            scalar_t const * current_line_copy= current_line;
            scalar_t const * const current_line_end = current_line + data_dimension;
            
            if(sign)
            {
              for(; current_line_copy < current_line_end; ++current_line_copy, ++p)
              {
                *p += *current_line_copy;
              }
            }
            else 
            {
              for(; current_line_copy < current_line_end; ++current_line_copy, ++p)
              {
                *p -= *current_line_copy;
              }
            }
          }
        }

        // posts the new value to the listeners
        if(update)
        {
          scalar_t *p(p_acc);
          for(size_t d = 0; d < data_dimension; d++, ++p)
          {
            *p *= 2;
          }          
          signal_acc(&accumulator);
        }

        // it could be that there is no update in the accumulator, still there might be one
        // in the references
        map_update[current_index] = static_cast<int>(count_current_index);
        signal_reference_count(&map_update);

        signal_counter();
      }

      //! Project the data onto the orthogonal subspace of the provided vector
      void project_onto_orthogonal_subspace(data_t const &mu)
      {
        // update of vectors in the orthogonal space, and update of the norms at the same time. 
        
        scalar_t const * const p_mu = &mu.data()[0];
        scalar_t * current_line = p_c_matrix;
        
        for(size_t line = 0; line < nb_elements; line ++, current_line += data_padding)
        {
          scalar_t const inner_prod = inner_product(p_mu, current_line);
          for(int column = 0; column < data_dimension; column++)
          {
            current_line[column] -= inner_prod * p_mu[column];
          }               
        }
  
        signal_counter();
      }

    };

}

#endif /* GRASSMANN_AVERAGES_PCA_CHUNK_PROCESSORS_HPP__ */
