// Copyright 2014, Max Planck Society.
// Distributed under the BSD 3-Clause license.
// (See accompanying file LICENSE.txt or copy at
// http://opensource.org/licenses/BSD-3-Clause)

/*!@file
 * This file is dedicated to a faster implementation of the acos function
 */



#include <boost/test/unit_test.hpp>

#include <test/test_main.hpp>
#include <include/private/utilities.hpp>


// boost chrono
#include <boost/chrono/include.hpp>


//! An approximation of the acos function with Chebyshev polynomials
template <class precision = double>
struct safe_acos_chebyshev_approximation
{
  typedef precision result_type;

  std::vector<precision> v_coefficients;
  precision const* p_coefficients;

  //! The number of coefficients used for evaluation, which could be lower than the 
  //! number of coefficients used for computing the Chebyshev coefficients.
  size_t truncation;

  safe_acos_chebyshev_approximation(size_t nb_coefficients = 25, size_t truncation_ = std::numeric_limits<size_t>::max()) : 
     v_coefficients(nb_coefficients), 
     truncation(std::min(truncation_, nb_coefficients)),
     p_coefficients(0)
  {
    using namespace std; // bringing acos

    // evaluating the first terms, which siplifies because of the acos(cos(x)) = x on this interval
    std::vector<double> v_fk(nb_coefficients, M_PI/nb_coefficients);
    for(size_t s = 0; s < nb_coefficients; s++)
    {
      v_fk[s] *= s + 0.5;
    }

    for(size_t s = 0; s < nb_coefficients; s++)
    {
      double acc(0);
      for(size_t j = 0; j < nb_coefficients; j++)
      {
        acc += v_fk[j] * cos(M_PI/nb_coefficients * s * (j + 0.5));
      }

      v_coefficients[s] = static_cast<precision>(2 * acc / nb_coefficients);
    }

    v_coefficients[0] /= 2; // this one is always * 0.5
    p_coefficients = &v_coefficients[0];
  }

  template <class T>
  result_type operator()(T v) const
  {
    // this branch is killing everything
    if(v >= 1) 
    {
      return 0;
    }
    else if(v <= -1)
    {
      return M_PI;
    }

#if 1
    precision const* p_coeff = p_coefficients;
    precision a = 1;
    precision b = v;
    v *= 2;
    precision acc = *(p_coeff++);  // p_coefficients[0] already multiplied by 0.5
    precision next;

    for(size_t s = 1; s < truncation - 1; s++, p_coeff++)
    {
      acc += (*p_coeff) * b;
      next = v * b - a;
      a = b;
      b = next;
    }

    return acc + (*p_coeff) * b;
#else
    // Clenshaw recurrence, see numerical recipes
    // actually slower
    precision const* p_coeff = p_coefficients + truncation - 1;
    precision d = 0;
    precision dd = 0;
    precision sv;
    precision v2 = 2 * v;
    

    for(; p_coeff > p_coefficients; p_coeff--)
    {
      sv = d;
      d = v2 * d - dd + *p_coeff;
      dd = sv;
    }

    return v * d -  dd + (*p_coeff); // p_coefficients[0] already multiplied by 0.5



#endif
  }
};



//! An approximation of the acos function with Chebyshev polynomials
template <class precision = double>
struct safe_acos_linear_interpolation
{
  typedef precision result_type;

  std::vector<precision> v_table;
  precision step;
  precision step_inv;

  safe_acos_linear_interpolation(size_t nb_coefficients = 100) : 
     v_table(nb_coefficients + 1)
  {
    using namespace std; // bringing acos

    step = static_cast<precision>(2./nb_coefficients);
    step_inv = static_cast<precision>(nb_coefficients / 2.);

    for(size_t s = 0; s < nb_coefficients; s++)
    {
      v_table[s] = acos(-1 + s * step);
    }
    v_table[nb_coefficients] = 0;
  }

  template <class T>
  result_type operator()(T v) const
  {
    // this branch is killing everything
    if(v >= 1) 
    {
      return 0;
    }
    else if(v <= -1)
    {
      return M_PI;
    }

    v = (v + 1) * step_inv;
    size_t index = static_cast<size_t>(v);

    precision a = v_table[index];
    precision b = v_table[index + 1];

    return a + (v - index) * (b - a);

  }
};


//! An approximation of the acos function with Chebyshev polynomials
template <class precision = double>
struct safe_acos_quadratic_interpolation
{
  typedef precision result_type;

  std::vector<precision> v_table;
  precision step;
  precision step_inv;

  safe_acos_quadratic_interpolation(size_t nb_coefficients = 50) : 
     v_table(nb_coefficients + 1)
  {
    using namespace std; // bringing acos

    step = static_cast<precision>(2./nb_coefficients);
    step_inv = static_cast<precision>(nb_coefficients / 2.);

    for(size_t s = 0; s < nb_coefficients; s++)
    {
      v_table[s] = acos(-1 + s * step);
    }
    v_table[nb_coefficients] = 0;
  }

  template <class T>
  result_type operator()(T v) const
  {
    // this branch is killing everything
    if(v >= 1) 
    {
      return 0;
    }
    else if(v <= -1)
    {
      return M_PI;
    }

    v = (v + 1) * step_inv;

    size_t index = static_cast<size_t>(v);
    if(index == 0)
    {
      return M_PI + v * (v_table[1] - M_PI);
    }

    precision y = v_table[index];
    precision yp1 = v_table[index + 1];
    precision ym1 = v_table[index - 1];

    v -= index;
    return (0.5*(yp1 + ym1) - y) * v*v + 0.5*(yp1 - ym1) * v + y;

  }
};


//! An approximation of the acos function with Chebyshev polynomials
template <class precision = double>
struct safe_acos_taylor_approximation_order1
{
  typedef precision result_type;

  std::vector< std::pair<precision, precision> > v_table;
  precision step;
  precision step_inv;

  safe_acos_taylor_approximation_order1(size_t nb_coefficients = 100) : 
     v_table(nb_coefficients + 1)
  {
    using namespace std; // bringing acos

    step = static_cast<precision>(2./nb_coefficients);
    step_inv = static_cast<precision>(nb_coefficients / 2.);

    for(size_t s = 0; s < nb_coefficients; s++)
    {
      precision v = -1 + s * step;
      v_table[s].first = acos(v);
      v_table[s].second = -step/sqrt(1-v*v);
    }
    v_table[nb_coefficients] = std::make_pair(0, 0);
  }

  template <class T>
  result_type operator()(T v) const
  {
    // this branch is killing everything
    if(v >= 1 - step) 
    {
      return v >= 1 ? 0 : acos(v);
    }
    else if(v <= -1 + step)
    {
      return v <= -1 ? M_PI : acos(v);
    }

    v = (v + 1) * step_inv;
    size_t index = static_cast<size_t>(v);

    v -= index;

    // little wart to cope for some discontinuities on junction points, not super
    return 0.5*(v_table[index].first + v_table[index].second * v + v_table[index+1].first + v_table[index+1].second * (v-1));

  }
};


//! An approximation of the acos function with Cublic splines
template <class precision = double>
struct safe_acos_cublic_spline_approximation
{
  typedef precision result_type;

  std::vector< std::pair<precision, precision> > v_table;
  precision step;
  precision step_inv;

  safe_acos_cublic_spline_approximation(size_t nb_coefficients = 30) : 
     v_table(nb_coefficients + 1)
  {
    using namespace std; // bringing acos

    step = static_cast<precision>(2./nb_coefficients);
    step_inv = static_cast<precision>(nb_coefficients / 2.);

    for(size_t s = 0; s < nb_coefficients; s++)
    {
      precision v = -1 + s * step;
      v_table[s].first = acos(v);
      v_table[s].second = -1/sqrt(1-v*v);
    }
    v_table[nb_coefficients] = std::make_pair(0, 0);
  }

  template <class T>
  result_type operator()(T v) const
  {
    // this branch is killing everything
    if(v >= 1 - step) 
    {
      return v >= 1 ? 0 : acos(v);
    }
    else if(v <= -1 + step)
    {
      return v <= -1 ? M_PI : acos(v);
    }

    v = (v + 1) * step_inv;
    size_t index = static_cast<size_t>(v);

    precision p0 = v_table[index].first;
    precision m0 = v_table[index].second;
    precision p1 = v_table[index+1].first;
    precision m1 = v_table[index+1].second;


    v -= index;
    assert(v >= 0 && v < 1);


    return (2*v*v*v -3*v*v + 1) * p0 + (v*v*v -2 * v*v + v) * m0 + (-2*v*v + 3 * v * v) * p1 + (v*v*v - v*v) * m1 ;

  }
};




BOOST_AUTO_TEST_CASE(acos_approximation_chebyshev_method)
{
  safe_acos_chebyshev_approximation<> acos_obj;

  {
    const int nb_values = 97;
    const double increment = 2./nb_values;
    double current = -1;
    for(size_t s = 0; s < nb_values; s++, current += increment)
    {
      BOOST_CHECK_CLOSE(acos_obj(current), acos(current), 1); // 1% radian
    }
  }

  const std::string filename = "./toto_chebyshev.txt";
  std::ofstream ff(filename.c_str());

  BOOST_REQUIRE(ff.is_open());

  for(int i = 0; i < 10000; i++)
  {
    ff << acos_obj(-1 + 2.*i/ 10000) << std::endl;
  }

}

BOOST_AUTO_TEST_CASE(acos_approximation_linear_method)
{
  safe_acos_linear_interpolation<> acos_obj;

  {
    const int nb_values = 97;
    const double increment = 2./nb_values;
    double current = -1;
    for(size_t s = 0; s < nb_values; s++, current += increment)
    {
      BOOST_CHECK_CLOSE(acos_obj(current), acos(current), 1); // 1% radian
    }
  }

  const std::string filename = "./toto_linear.txt";
  std::ofstream ff(filename.c_str());

  BOOST_REQUIRE(ff.is_open());

  for(int i = 0; i < 10000; i++)
  {
    ff << acos_obj(-1 + 2.*i/ 10000) << std::endl;
  }

}


BOOST_AUTO_TEST_CASE(acos_approximation_quadratic_method)
{
  safe_acos_quadratic_interpolation<> acos_obj;

  {
    const int nb_values = 97;
    const double increment = 2./nb_values;
    double current = -1;
    for(size_t s = 0; s < nb_values; s++, current += increment)
    {
      BOOST_CHECK_CLOSE(acos_obj(current), acos(current), 1); // 1% radian
    }
  }


  const std::string filename = "./toto_quad.txt";
  std::ofstream ff(filename.c_str());

  BOOST_REQUIRE(ff.is_open());

  for(int i = 0; i < 10000; i++)
  {
    ff << acos_obj(-1 + 2.*i/ 10000) << std::endl;
  }


}

BOOST_AUTO_TEST_CASE(acos_approximation_taylor1_method)
{
  safe_acos_taylor_approximation_order1<> acos_obj;

  {
    const int nb_values = 97;
    const double increment = 2./nb_values;
    double current = -1;
    for(size_t s = 0; s < nb_values; s++, current += increment)
    {
      BOOST_CHECK_CLOSE(acos_obj(current), acos(current), 1); // 1% radian
    }
  }

  const std::string filename = "./toto_taylor1.txt";
  std::ofstream ff(filename.c_str());

  BOOST_REQUIRE(ff.is_open());

  for(int i = 0; i < 10000; i++)
  {
    ff << acos_obj(-1 + 2.*i/ 10000) << std::endl;
  }

}


BOOST_AUTO_TEST_CASE(acos_approximation_cubic_spline_method)
{
  safe_acos_cublic_spline_approximation<> acos_obj;

  {
    const int nb_values = 97;
    const double increment = 2./nb_values;
    double current = -1;
    for(size_t s = 0; s < nb_values; s++, current += increment)
    {
      BOOST_CHECK_CLOSE(acos_obj(current), acos(current), 1); // 1% radian
    }
  }

  const std::string filename = "./toto_cubic_spline.txt";
  std::ofstream ff(filename.c_str());

  BOOST_REQUIRE(ff.is_open());

  for(int i = 0; i < 10000; i++)
  {
    ff << acos_obj(-1 + 2.*i/ 10000) << std::endl;
  }

}
BOOST_AUTO_TEST_CASE(acos_approximation_chebyshev_method_performances)
{
  typedef float precision;
  safe_acos_chebyshev_approximation<precision> acos_obj(10);
  grassmann_averages_pca::details::safe_acos<precision> acos_obj_safe;

  const size_t repetition = 10000;
  const size_t nb_values = 1000;
  const precision increment = 2./nb_values;


  typedef boost::chrono::steady_clock clock_type;
  clock_type::duration elapsed1, elapsed2;

  precision accg = 0;

  // acos
  {
    clock_type::time_point start = clock_type::now();


    precision current = -1;
    for(size_t s = 0; s < nb_values; s++, current += increment)
    {
      float acc = 0;
      for(int i = 0; i < repetition; i++)
      {
        acc += acos_obj_safe(current);
      }
      accg += acc;
    }
    elapsed1 = clock_type::now() - start;

  }


  // acos
  {
    clock_type::time_point start = clock_type::now();


    precision current = -1;
    for(size_t s = 0; s < nb_values; s++, current += increment)
    {
      precision acc = 0;
      for(int i = 0; i < repetition; i++)
      {
        acc += acos_obj(current);
      }
      accg += acc;
    }

    elapsed2 = clock_type::now() - start;
  }

  std::cout << "processing " << nb_values * repetition << " elements " << std::endl
            << " - acos " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed1) << std::endl
            << " - Chebyshev approximation of acos " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed2) << std::endl
            << " - dummy value " << (accg > 0 ? 1: -1) << std::endl;



}



BOOST_AUTO_TEST_CASE(acos_approximation_linear_method_performances)
{
  typedef float precision;
  safe_acos_linear_interpolation<precision> acos_obj(10);
  grassmann_averages_pca::details::safe_acos<precision> acos_obj_safe;

  const size_t repetition = 10000;
  const size_t nb_values = 1000;
  const precision increment = 2./nb_values;


  typedef boost::chrono::steady_clock clock_type;
  clock_type::duration elapsed1, elapsed2;

  precision accg = 0;

  // acos
  {
    clock_type::time_point start = clock_type::now();


    precision current = -1;
    for(size_t s = 0; s < nb_values; s++, current += increment)
    {
      float acc = 0;
      for(int i = 0; i < repetition; i++)
      {
        acc += acos_obj_safe(current);
      }
      accg += acc;
    }
    elapsed1 = clock_type::now() - start;

  }


  // acos
  {
    clock_type::time_point start = clock_type::now();


    precision current = -1;
    for(size_t s = 0; s < nb_values; s++, current += increment)
    {
      precision acc = 0;
      for(int i = 0; i < repetition; i++)
      {
        acc += acos_obj(current);
      }
      accg += acc;
    }

    elapsed2 = clock_type::now() - start;
  }

  std::cout << "processing " << nb_values * repetition << " elements " << std::endl
            << " - acos " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed1) << std::endl
            << " - Linear approximation of acos " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed2) << std::endl
            << " - dummy value " << (accg > 0 ? 1: -1) << std::endl;



}




BOOST_AUTO_TEST_CASE(acos_approximation_quadratic_method_performances)
{
  typedef float precision;
  safe_acos_quadratic_interpolation<precision> acos_obj(10);
  grassmann_averages_pca::details::safe_acos<precision> acos_obj_safe;

  const size_t repetition = 10000;
  const size_t nb_values = 1000;
  const precision increment = 2./nb_values;


  typedef boost::chrono::steady_clock clock_type;
  clock_type::duration elapsed1, elapsed2;

  precision accg = 0;

  // acos
  {
    clock_type::time_point start = clock_type::now();


    precision current = -1;
    for(size_t s = 0; s < nb_values; s++, current += increment)
    {
      float acc = 0;
      for(int i = 0; i < repetition; i++)
      {
        acc += acos_obj_safe(current);
      }
      accg += acc;
    }
    elapsed1 = clock_type::now() - start;

  }


  // acos
  {
    clock_type::time_point start = clock_type::now();


    precision current = -1;
    for(size_t s = 0; s < nb_values; s++, current += increment)
    {
      precision acc = 0;
      for(int i = 0; i < repetition; i++)
      {
        acc += acos_obj(current);
      }
      accg += acc;
    }

    elapsed2 = clock_type::now() - start;
  }

  std::cout << "processing " << nb_values * repetition << " elements " << std::endl
            << " - acos " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed1) << std::endl
            << " - Quadratic approximation of acos " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed2) << std::endl
            << " - dummy value " << (accg > 0 ? 1: -1) << std::endl;



}



BOOST_AUTO_TEST_CASE(acos_approximation_taylor1_method_performances)
{
  typedef float precision;
  safe_acos_taylor_approximation_order1<precision> acos_obj(100);
  grassmann_averages_pca::details::safe_acos<precision> acos_obj_safe;

  const size_t repetition = 10000;
  const size_t nb_values = 1000;
  const precision increment = 2./nb_values;


  typedef boost::chrono::steady_clock clock_type;
  clock_type::duration elapsed1, elapsed2;

  precision accg = 0;

  // acos
  {
    clock_type::time_point start = clock_type::now();


    precision current = -1;
    for(size_t s = 0; s < nb_values; s++, current += increment)
    {
      float acc = 0;
      for(int i = 0; i < repetition; i++)
      {
        acc += acos_obj_safe(current);
      }
      accg += acc;
    }
    elapsed1 = clock_type::now() - start;

  }


  // acos
  {
    clock_type::time_point start = clock_type::now();


    precision current = -1;
    for(size_t s = 0; s < nb_values; s++, current += increment)
    {
      precision acc = 0;
      for(int i = 0; i < repetition; i++)
      {
        acc += acos_obj(current);
      }
      accg += acc;
    }

    elapsed2 = clock_type::now() - start;
  }

  std::cout << "processing " << nb_values * repetition << " elements " << std::endl
            << " - acos " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed1) << std::endl
            << " - Taylor approximation of acos " << boost::chrono::duration_cast<boost::chrono::microseconds>(elapsed2) << std::endl
            << " - dummy value " << (accg > 0 ? 1: -1) << std::endl;



}