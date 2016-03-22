// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifdef ENABLE_MPI
# undef ENABLE_MPI
#endif
#define ENABLE_MPI 0

#include "config.h"

#if HAVE_ALUGRID

# include <dune/grid/alugrid.hh>

# include "problems/AO2013.hh"
# include "eocexpectations.hh"

namespace Dune {
namespace GDT {
namespace Test {

// polorder 1, conforming

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     1,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {}; // {5.33e-02, 1.69e-02};
      else
        return {7.42e-03, 4.32e-03, 1.19e-03, 2.53e-04}; // {5.57e-02, 1.99e-02, 5.54e-03, 1.29e-03};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {}; // {3.82e-01, 2.29e-01};
      else
        return {3.93e-01, 3.22e-01, 1.71e-01, 7.75e-02}; // {4.32e-01, 2.93e-01, 1.50e-01, 6.54e-02};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {};
      else
        return {8.47e-02, 5.12e-02, 2.66e-02, 1.21e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, conforming

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     2,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {}; // {1.18e-02, 2.12e-03};
      else
        return {}; // {1.18e-02, 2.11e-03, 3.89e-04, 7.76e-05};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {}; // {1.67e-01, 5.58e-02};
      else
        return {}; // {1.69e-01, 5.96e-02, 1.94e-02, 6.04e-03};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {};
      else
        return {};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 1, nonconforming

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     1,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {};
      else
        return {1.18e-02, 3.50e-03, 8.48e-04, 1.72e-04};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {};
      else
        return {5.40e-01, 2.86e-01, 1.37e-01, 5.96e-02};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {};
      else
        return {7.82e-02, 4.12e-02, 2.06e-02, 9.28e-03};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations



template class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              1 >;
template class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              2 >;
template class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              1 >;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
