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

void number2filename(int file_number, char *filename)
{
  const char *fn_template = "/is/ps/shared/users/jonas/movies/starwars_%.3d/frame%.7d.png";
  const int dir_num = file_number / 10000;
  sprintf(filename, fn_template, dir_num, file_number);
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




  char filename[100];
  cv::Mat image;
  

  // Read first image to get image size
  number2filename(1, filename);
  image = cv::imread(filename, CV_LOAD_IMAGE_COLOR);
  cv::Size image_size = image.size();
  const size_t rows = image_size.height;
  const size_t cols = image_size.width;

  // Allocate data
  std::cout << "=== Allocate ===" << std::endl;
  std::cout << "  Number of images: " << num_frames << std::endl;
  std::cout << "  Image size:       " << cols << "x" << rows << " (RGB)" << std::endl;
  std::cout << "  Data size:        " << (num_frames*rows*cols * 3) / (1024 * 1024) << " MB" << std::endl;
  float *data = new float[num_frames*rows*cols * 3];

  size_t idx = 0;
  for(size_t file_number = 1; file_number <= num_frames; file_number++)
  {
    number2filename(file_number, filename);
    image = cv::imread(filename, CV_LOAD_IMAGE_COLOR);
    if(!image.data)
    {
      std::cerr << "error: could not load image '" << filename << "'" << std::endl;
      break;
    }
    else
    {
      for(int r = 0; r < rows; r++)
      {
        for(int c = 0; c < cols; c++)
        {
          // XXX: Is this the right way to store the data for TGA?
          const cv::Vec3b intensity = image.at<cv::Vec3b>(r, c);
          data[idx] = intensity.val[0]; // blue
          data[idx + 1] = intensity.val[1]; // green
          data[idx + 2] = intensity.val[2]; // red
          idx += 3;
        }
      }
    }
  }

  // Compute trimmed Grassmann Average
  // XXX: Ask Raffi for help here
  // XXX: We need to time this computation


  namespace ub = boost::numeric::ublas;

  using namespace grassmann_averages_pca;
  using namespace grassmann_averages_pca::details::ublas_helpers;


  typedef float input_array_type;

  typedef external_storage_adaptor<input_array_type> input_storage_t;
  typedef ub::matrix<input_array_type, ub::column_major, input_storage_t> input_matrix_t;

  typedef external_storage_adaptor<input_array_type> output_storage_t;
  typedef ub::matrix<input_array_type, ub::row_major, output_storage_t> output_matrix_t; // this is in fact column_major, it should be in accordance with the
                                                                               // dimension of the matrix output_basis_vectors (we take the transpose of it)

  const size_t dimension = image_size.width * image_size.height;
  const size_t nb_elements = num_frames;


  // input data matrix, external storage.
  input_storage_t input_storage(nb_elements * dimension, static_cast<input_array_type *>(mxGetData(X)));
  input_matrix_t input_data(nb_elements, dimension, input_storage);

  // output data matrix, also external storage for uBlas
  output_storage_t storageOutput(dimension * max_dimension, static_cast<input_array_type *>(mxGetData(outputMatrix)));
  output_matrix_t output_basis_vectors(max_dimension, dimension, storageOutput);


  // input data matrix, external storage.
  input_storage_t input_storage(nb_elements*dimension, static_cast<input_array_type*>(mxGetData(X)));
  input_matrix_t input_data(nb_elements, dimension, input_storage);

  // output data matrix, also external storage for uBlas
  output_storage_t storageOutput(dimension * max_dimension, static_cast<input_array_type *>(mxGetData(outputMatrix)));
  output_matrix_t output_basis_vectors(max_dimension, dimension, storageOutput);





  // this is the form of the data extracted from the storage
  typedef ub::vector<input_array_type> data_t;
  typedef grassmann_pca_with_trimming< data_t > grassmann_pca_with_trimming_t;

  typedef row_iter<const input_matrix_t> const_input_row_iter_t;
  typedef row_iter<output_matrix_t> output_row_iter_t;


  // main instance
  grassmann_pca_with_trimming_t instance(trimming_percentage / 100);


  if(nb_processors > 0)
  {
    if(!instance.set_nb_processors(nb_processors))
    {
      std::cerr << "[configuration]" << "Incorrect number of processors. Please consult the documentation (was " << nb_processors << ")" << std::endl;
      return 1;
    }
  }


#if 0

  // initialisation vector if given
  std::vector<data_t> init_vectors;
  if(algorithm_configuration.initial_vectors != 0)
  { 
    init_vectors.resize(max_dimension);
    input_storage_t input_init_vector_storage(max_dimension*dimension, static_cast<input_array_type*>(mxGetData(algorithm_configuration.initial_vectors)));
    input_matrix_t input_init_vector_data(dimension, max_dimension, input_init_vector_storage);
    for(size_t index = 0;
        index < max_dimension;
        index++)
    {
      init_vectors[index] = ub::column(input_init_vector_data, index);
    }

    // if the initial vectors are set, we avoid the computation of the regular PCA.
    nb_pca_steps = 0;
  }
#endif

  if(!instance.set_nb_steps_pca(nb_pca_steps))
  {
    std::cerr << "[configuration]" << "Incorrect number of regular PCA steps. Please consult the documentation (was " << nb_pca_steps << ")" << std::endl;
    return false;
  }


  bool ret = instance.batch_process(
    max_iterations,
    max_dimension,
    const_input_row_iter_t(input_data, 0),
    const_input_row_iter_t(input_data, input_data.size1()),
    output_row_iter_t(output_basis_vectors, 0),
    0);












  // Save resulting components
  // XXX: Ask Raffi for help here

  // Delete data
  delete[] data;

  return 0;
}

