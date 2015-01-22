// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)


/*!@file
 * This file contains an application of the EM-PCA to the frames of a movie, in order to compare the results
 * with the Grassmann average or the trimmed grassmann average. You may adapt the code to your needs by modifying
 * the function @c number2filename, which from the index of the frame returns a full path of the file of this frame.
 * To limit the memory footprint, the data is stored directly in the temporary memory of the algorithm as it is loaded
 * (see @c iterator_on_image_files). An observer flushes the results to the disk as they arrive from the algorithm (see 
 * @c grassmann_pca_observer).
 */

#include <cstdio>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>


#include <include/private/em_pca.hpp>
#include <include/private/boost_ublas_external_storage.hpp>
#include <include/private/boost_ublas_row_iterator.hpp>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include <boost/iterator/iterator_facade.hpp>

#include <string>
#include <fstream>

#ifndef MAX_PATH
  // "funny" differences win32/posix
  #define MAX_PATH PATH_MAX
#endif

static std::string movielocation = "/is/ps/shared/users/jonas/movies/";
static std::string eigenvectorslocation = "./";


namespace grassmann_averages_pca
{
  namespace applications
  {
    //! Transforms a frame index to a full path name.
    void number2filename(size_t file_number, char *filename)
    {
      const size_t dir_num = file_number / 10000;
      sprintf(filename, 
              (movielocation + "/starwars_%.3d/frame%.7d.png").c_str(), 
              dir_num, 
              file_number);
    }

    //! Iterator loading the images on demand instead of storing everything on memory
    template <class T>
    class iterator_on_image_files : 
      public boost::iterator_facade<
            iterator_on_image_files<T>
          , boost::numeric::ublas::vector<T>
          , std::random_access_iterator_tag
          , boost::numeric::ublas::vector<T> const& // const reference
        >
    {
    public:
      typedef iterator_on_image_files<T> this_type;
      typedef boost::numeric::ublas::vector<T> image_vector_type;
      iterator_on_image_files() : m_index(std::numeric_limits<size_t>::max()) 
      {}

      explicit iterator_on_image_files(size_t index)
        : m_index(index) 
      {}

    private:
      friend class boost::iterator_core_access;

      typename this_type::difference_type distance_to(this_type const& r) const
      {
        return typename this_type::difference_type(r.m_index) - typename this_type::difference_type(m_index); // sign promotion
      }

      void increment() 
      { 
        m_index++; 
        image_vector.resize(0, false);
      }

      bool equal(this_type const& other) const
      {
        return this->m_index == other.m_index;
      }

      image_vector_type const& dereference() const 
      { 
        if(image_vector.empty())
        {
          read_image();
        }
        return image_vector; 
      }
  
  
      void read_image() const
      {
        char filename[MAX_PATH];
        number2filename(m_index, filename);
    
        if((m_index % 1000) == 0)
        {
          lock_t guard(internal_mutex);
          std::cout << "[THREAD " << boost::this_thread::get_id() << "] Reading " << filename << std::endl;
        }
    
        cv::Mat image = cv::imread(filename, CV_LOAD_IMAGE_COLOR);
        if(!image.data)
        {
          std::ostringstream o;
          o << "error: could not load image '" << filename << "'";
          std::cerr << o.str() << std::endl;
          throw std::runtime_error(o.str());
        }

        const int w = image.size().width;
        const int h = image.size().height;

        image_vector.resize(w * h * 3);
        typename boost::numeric::ublas::vector<T>::iterator it = image_vector.begin();

        for(int y = 0; y < h; y++)
        {
          for(int x = 0; x < w; x++)
          {
            cv::Vec3b pixel = image.at<cv::Vec3b>(y, x);
            *it++ = pixel[0];
            *it++ = pixel[1];
            *it++ = pixel[2];
          }
        }
      }

      void advance(typename this_type::difference_type n)
      {
        if(n < 0)
        {
          assert((-n) <= static_cast<typename this_type::difference_type>(m_index));
          m_index -= n;
        }
        else
        {
          //assert(n + index <= matrix->size1());
          m_index += n;
        }
        if(n != 0)
        {
          image_vector.resize(0, false);
        }
      }  

      size_t m_index;
      mutable boost::numeric::ublas::vector<T> image_vector;

      // for being able to log in a thread safe manner
      typedef boost::recursive_mutex mutex_t;
      typedef boost::lock_guard<mutex_t> lock_t;

      static mutex_t internal_mutex;


    };

    template <class T>
    typename iterator_on_image_files<T>::mutex_t iterator_on_image_files<T>::internal_mutex;

    //! Observer for saving the results of the algorithm as they arrive, and monitor the progress
    //! for very long runs.
    template <class data_t>
    struct grassmann_pca_observer
    {
    private:
      size_t element_per_line_during_save;
      size_t last_nb_iteration;
      std::string filename_from_template(std::string template_, size_t i) const
      {
        char filename[MAX_PATH];
        sprintf(filename, template_.c_str(), i);
        return filename;
      }

      void save_vector(const data_t& v, std::string filename) const
      {
        std::ofstream f(filename);
        if(!f.is_open())
        {
          std::cerr << "[ERROR] Cannot open the file " << filename << " for writing" << std::endl;
          return;
        }
    
        std::cout << "-\tWriting file " << filename;
    
        typedef typename data_t::const_iterator element_iterator;
    
        element_iterator itelement(v.begin());
        for(int i = 0; i < v.size(); i++, ++itelement)
        {
          if((i + 1) % element_per_line_during_save == 0)
          {
            f << std::endl;
          }
      
          f << *itelement << " ";
        }
    
        f.close();
        std::cout << " -- done" << std::endl;
      }

    public:
      grassmann_pca_observer(size_t element_per_line_during_save_) : 
        element_per_line_during_save(element_per_line_during_save_)
      {}


      void log_error_message(const char* message) const
      {
        std::cout << message << std::endl;
      }


      //! This is called after centering the data in order to keep track 
      //! of the mean of the dataset
      void signal_mean(const data_t& mean) const
      {
        std::cout << "* Mean computed" << std::endl;
        save_vector(mean, "./mean_vector.txt");
      }

      //! Called after the computation of the PCA
      void signal_pca(const data_t& pca,
                      size_t current_eigenvector_dimension) const
      {
        std::cerr << "[ERROR] PCA callback called: this is inconsistent with the algorithm" << std::endl;
        throw std::runtime_error("[ERROR] PCA callback called: inconsistent with the requested algorithm");
      }

      //! Called each time a new eigenvector is computed
      void signal_eigenvector(const data_t& current_eigenvector, 
                              size_t current_eigenvector_dimension) const
      {
        std::cout << "* Eigenvector subspace " << current_eigenvector_dimension << " computed in # " << last_nb_iteration << " iterations " << std::endl;
        save_vector(current_eigenvector, filename_from_template("./vector_empca_subspace_%.7d.txt", current_eigenvector_dimension));
      }

      //! Called at every step of the algorithm, at the end of the step
      void signal_intermediate_result(
        const data_t& current_eigenvector_state, 
        size_t current_eigenvector_dimension,
        size_t current_iteration_step) 
      {
        last_nb_iteration = current_iteration_step;
        if((current_iteration_step % 100) == 0)
        {
          std::cout << "* PCA subspace " << current_eigenvector_dimension << " @ iteration " << current_iteration_step << std::endl;
        }
      }

    };


  }
}





int main(int argc, char *argv[])
{
  namespace po = boost::program_options;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("nb-frames,f",       po::value<int>(),           "number of frames in the movie")
    ("movie,m",           po::value<std::string>(),   "full path to the movie location")
    ("max-dimensions,d",  po::value<int>(),           "requested number of components for the computation of the trimmed grassmann average")
    ("max-iterations",    po::value<int>(),           "maximum number of iterations (defaults to the number of frames)")
    ("nb-processors",      po::value<int>(),          "number of processors used (defaults to 1)")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if(vm.count("help")) 
  {
    std::cout << desc << "\n";
    return 1;
  }


  size_t num_frames = 100; //179415;
  if(vm.count("nb-frames")) 
  {
    std::cout << "The movie contains #" << vm["nb-frames"].as<int>() << " frames.\n";
    num_frames = vm["nb-frames"].as<int>();
  } 
  else 
  {
    std::cout << "Nb of frames not set, setting it to" << num_frames << "\n";
  }
  
  if(vm.count("movie")) 
  {
    std::cout << "Directory of the movie #" << vm["movie"].as<std::string>() << "\n";
    movielocation = vm["movie"].as<std::string>();
  }   

  size_t max_dimension = 30; 
  if(vm.count("max-dimensions")) 
  {
    std::cout << "Number of components requested #" << vm["max-dimensions"].as<int>() << "\n";
    max_dimension = vm["max-dimensions"].as<int>();
  } 
  else 
  {
    std::cout << "Nb of components not set, setting it to" << max_dimension << "\n";
  }

  size_t max_iterations = num_frames;
  if(vm.count("max-iterations")) 
  {
    std::cout << "Maximum number of iterations #" << vm["max-iterations"].as<int>() << "\n";
    max_iterations = vm["max-iterations"].as<int>();
  } 

  int nb_processors = 0;
  if(vm.count("nb-processors")) 
  {
    std::cout << "Number of processors #" << vm["nb-processors"].as<int>() << "\n";
    nb_processors = vm["nb-processors"].as<int>();
  } 




  namespace ub = boost::numeric::ublas;
  using namespace grassmann_averages_pca;
  using namespace grassmann_averages_pca::applications;
  using namespace grassmann_averages_pca::details::ublas_helpers;

  // Reading the first image to have the dimensions
  size_t rows(0);
  size_t cols(0);
   
  {
    char filename[MAX_PATH];
    number2filename(1, filename);

    cv::Mat image = cv::imread(filename, CV_LOAD_IMAGE_COLOR);
    cv::Size image_size = image.size();
    rows = image_size.height;
    cols = image_size.width;
  }

  
  
  // type of the scalars manipulated by the algorithm
  typedef float input_array_type;


  // Allocate data
  std::cout << "=== Allocate ===" << std::endl;
  std::cout << "  Number of images: " << num_frames << std::endl;
  std::cout << "  Image size:       " << cols << "x" << rows << " (RGB)" << std::endl;
  std::cout << "  Data size:        " << (num_frames*rows*cols * 3 * sizeof(input_array_type)) / (1024 * 1024) << " MB" << std::endl;
    
  iterator_on_image_files<float> iterator_file_begin(1);
  iterator_on_image_files<float> iterator_file_end(num_frames + 1);


  // type of the data extracted from the input iterators
  typedef ub::vector<input_array_type> data_t;
  // type of the observer
  typedef grassmann_pca_observer<data_t> observer_t;
  // type of the em-pca algorithm
  typedef em_pca< data_t, observer_t > em_pca_t;

  // main instance
  em_pca_t instance;
  
  // storage for computed subspaces
  typedef std::vector<data_t> output_eigenvector_collection_t;
  output_eigenvector_collection_t v_output_eigenvectors(max_dimension);


  if(nb_processors > 0)
  {
    if(!instance.set_nb_processors(nb_processors))
    {
      std::cerr << "[configuration]" << "Incorrect number of processors. Please consult the documentation (was " << nb_processors << ")" << std::endl;
      return 1;
    }
  }

  // setting the observer
  observer_t my_simple_observer(cols);
  if(!instance.set_observer(&my_simple_observer))
  {
    std::cerr << "[configuration]" << "Error while setting the observer" << std::endl;
    return 1;
  }

  // requesting the centering of the data
  if(!instance.set_centering(true))
  {
    std::cerr << "[configuration]" << "Error while configuring the centering" << std::endl;
    return 1;
  }

  // running the computation
  bool ret = instance.batch_process(
    max_iterations,
    max_dimension,
    iterator_file_begin,
    iterator_file_end,
    v_output_eigenvectors.begin());

  // Results are saved in the observer instance

  if(!ret)
  {
    std::cerr << "The process returned an error" << std::endl;
    return 1;
  }


  return 0;
}

