#include <cstdio>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>


#include <include/grassmann_pca_with_trimming.hpp>
#include <include/private/boost_ublas_external_storage.hpp>
#include <include/private/boost_ublas_row_iterator.hpp>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include <boost/iterator/iterator_facade.hpp>

#include <string>
#include <fstream>

static std::string movielocation = "/is/ps/shared/users/jonas/movies/";
static std::string eigenvectorslocation = "./";

void number2filename(size_t file_number, char *filename)
{
  const char *fn_template = (movielocation + "starwars_%.3u/frame%.7u.png").c_str();
  const size_t dir_num = file_number / 10000;
  sprintf(filename, fn_template, dir_num, file_number);
}


//! An container over the content of a cvMat, addressable
//! in a vector fashion.
template <class T>
struct cvMatAsVector
{
  cvMatAsVector(const cv::Mat& image_) : 
    image(image_),
    width(image_.size().width),
    height(image_.size().height)
  {}
  
  
  T operator()(std::size_t index) const
  {
    int colorchannel = index % 3;
    index /= 3;
    int column = index % width;
    index /= width;
    return image.at<cv::Vec3b>(index, column)[colorchannel];
  }

  struct internal_iterator
  {
    cvMatAsVector<T> &outer_instance;
    
    internal_iterator() : current_index(0){}
    internal_iterator(size_t index) : current_index(index) {}

    T const& operator*() const
    {
      return outer_instance[current_index];
    }

    internal_iterator& operator++()
    {
      current_index++;
      return *this;
    }

    size_t current_index;
  };

  internal_iterator begin() const
  {
    return internal_iterator(0);
  }

  internal_iterator end() const
  {
    return internal_iterator(size());
  }
  
  size_t size() const
  {
    return width * height * 3;
  }
  
  const cv::Mat& image;
  const size_t width;
  const size_t height;
  
};


//! An iterator that will load the images on demand instead of storing everything on memory
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
    image_vector.clear();
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
    std::cout << "Reading " << filename << std::endl;
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
      image_vector.clear();
    }
  }  

  size_t m_index;
  mutable boost::numeric::ublas::vector<T> image_vector;
};

//! An output vector that also writes to a file
template <class data_iterator_t>
class iterator_on_output_data : 
  public boost::iterator_facade<
        iterator_on_output_data<data_iterator_t>
      , typename data_iterator_t::value_type
      , std::random_access_iterator_tag
      , typename data_iterator_t::reference
    >
{

public:
  typedef typename data_iterator_t::reference reference;
  typedef iterator_on_output_data<data_iterator_t> this_type;
  
  
  iterator_on_output_data() : 
    m_index(std::numeric_limits<size_t>::max()), 
    m_max_element_per_line(100),
    internal_it()
  {}

  explicit iterator_on_output_data(size_t index, data_iterator_t it_) : 
    m_index(index),
    m_max_element_per_line(100),
    internal_it(it_)
  {}

  void save_eigenvector()
  {
    char filename[MAX_PATH];
    const char * eigen_vector_template = "eigenvector_%.7d.txt";
    sprintf(filename, eigen_vector_template, m_index);
    
    std::ofstream f(filename);
    if(!f.is_open())
    {
      std::cerr << "Cannot open the file " << filename << " for writing" << std::endl;
      return;
    }
    
    std::cout << "Writing eigenvector file " << filename << std::endl;
    
    typedef typename data_iterator_t::value_type::iterator element_iterator;
    
    element_iterator itelement(internal_it->begin());
    for(int i = 0; i < internal_it->size(); i++, ++itelement)
    {
      if((i + 1) % m_max_element_per_line == 0)
      {
        f << std::endl;
      }
      
      f << *itelement;
    }
    
    f.close();
    std::cout << "Writing eigenvector file " << filename << " -- ok " << std::endl;
  }

private:
  friend class boost::iterator_core_access;

  typename this_type::difference_type distance_to(this_type const& r) const
  {
    return typename this_type::difference_type(r.m_index) - typename this_type::difference_type(m_index); // sign promotion
  }

  void increment() {
    // here we save before going further, except if it has not been changed
    save_eigenvector();
    m_index++; 
    ++internal_it;
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
    
  }  
  

  bool equal(this_type const& other) const
  {
    return this->internal_it == other.internal_it;
  }

  reference dereference() const
  { 
    return *internal_it;
  }
  
  


  size_t m_index;
  // max number of element on one line during the save
  size_t m_max_element_per_line;
  data_iterator_t internal_it;
};

#if 0
template <class data_t>
class data_and_save : data_t
{
  size_t index;
  
  //! sets the index that will determine the filename to which this eigenvector will 
  //! be saved.
  void set_index(size_t index_)
  {
    index = index_;
  }
  
  data_t& operator*()
};
#endif

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
    ("nb-pca-steps",      po::value<int>(),           "number of pca steps")
    ("trimming-percentage",      po::value<float>(),  "percentage of trimming")
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

  if(!vm.count("trimming-percentage")) 
  {
    std::cerr << "the parameter 'trimming-percentage' should be set, now exiting...\n";
    return 1;
  }
  const float trimming_percentage = vm["trimming-percentage"].as<float>();

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

  size_t nb_pca_steps = 3;
  if(vm.count("nb-pca-steps")) 
  {
    std::cout << "Number of PCA steps #" << vm["nb-pca-steps"].as<int>() << "\n";
    nb_pca_steps = vm["nb-pca-steps"].as<int>();
  } 
  else 
  {
    std::cout << "Number of PCA steps not set, setting it to" << nb_pca_steps << "\n";
  }


  int nb_processors = 0;
  if(vm.count("nb-processors")) 
  {
    std::cout << "Number of processors #" << vm["nb-processors"].as<int>() << "\n";
    nb_processors = vm["nb-processors"].as<int>();
  } 





  
  size_t rows(0);
  size_t cols(0);
   
  {
    char filename[MAX_PATH];
    // Read first image to get image size
    number2filename(1, filename);

    cv::Mat image = cv::imread(filename, CV_LOAD_IMAGE_COLOR);
    cv::Size image_size = image.size();
    rows = image_size.height;
    cols = image_size.width;
  }
  
  // Allocate data
  std::cout << "=== Allocate ===" << std::endl;
  std::cout << "  Number of images: " << num_frames << std::endl;
  std::cout << "  Image size:       " << cols << "x" << rows << " (RGB)" << std::endl;
  std::cout << "  Data size:        " << (num_frames*rows*cols * 3) / (1024 * 1024) << " MB" << std::endl;
    
  iterator_on_image_files<float> iterator_file_begin(1);
  iterator_on_image_files<float> iterator_file_end(num_frames + 1);

  // Compute trimmed Grassmann Average
  // XXX: Ask Raffi for help here
  // XXX: We need to time this computation


  namespace ub = boost::numeric::ublas;

  // type for the final eigen vectors
  typedef ub::vector<double> data_t;  


  using namespace grassmann_averages_pca;
  using namespace grassmann_averages_pca::details::ublas_helpers;


  typedef double input_array_type;

  const size_t dimension = rows * cols;
  const size_t nb_elements = num_frames;


  // this is the form of the data extracted from the storage
  typedef ub::vector<input_array_type> data_t;
  typedef grassmann_pca_with_trimming< data_t > grassmann_pca_with_trimming_t;


  // main instance
  grassmann_pca_with_trimming_t instance(trimming_percentage / 100);
  
  
  typedef std::vector<data_t> output_eigenvector_collection_t;
  output_eigenvector_collection_t v_output_eigenvectors(nb_elements);


  if(nb_processors > 0)
  {
    if(!instance.set_nb_processors(nb_processors))
    {
      std::cerr << "[configuration]" << "Incorrect number of processors. Please consult the documentation (was " << nb_processors << ")" << std::endl;
      return 1;
    }
  }



  if(!instance.set_nb_steps_pca(nb_pca_steps))
  {
    std::cerr << "[configuration]" << "Incorrect number of regular PCA steps. Please consult the documentation (was " << nb_pca_steps << ")" << std::endl;
    return false;
  }


  bool ret = instance.batch_process(
    max_iterations,
    max_dimension,
    iterator_file_begin,
    iterator_file_end,
    iterator_on_output_data<output_eigenvector_collection_t::iterator>(0, v_output_eigenvectors.begin()));


  // Save resulting components
  // XXX: Ask Raffi for help here


  return 0;
}

