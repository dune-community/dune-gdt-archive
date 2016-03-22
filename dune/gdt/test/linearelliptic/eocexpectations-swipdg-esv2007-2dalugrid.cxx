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

# include "problems/ESV2007.hh"
# include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {

// polorder 1, conforming alugrid

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     1,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {2.32e-02, 4.53e-03, 1.12e-03, 2.78e-04}; // {1.15e-01, 3.04e-02, 7.51e-03, 1.86e-03};
    else if (type == "H1_semi" || type == "energy")
      return {3.32e-01, 1.62e-01, 8.04e-02, 4.01e-02}; // {3.79e-01, 1.90e-01, 9.38e-02, 4.67e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, conforming alugrid

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     2,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 2 >
{
  typedef LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {}; // {1.25e-02, 1.42e-03, 1.69e-04, 2.08e-05};
    else if (type == "H1_semi" || type == "energy")
      return {}; // {7.84e-02, 2.01e-02, 5.02e-03, 1.26e-03};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 1, nonconforming alugrid

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     1,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {2.32e-02, 5.97e-03, 1.50e-03, 3.76e-04};
    else if (type == "H1_semi" || type == "energy")
      return {3.32e-01, 1.63e-01, 8.07e-02, 4.01e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations



template class LinearEllipticEocExpectations< LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              1 >;
template class LinearEllipticEocExpectations< LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              2 >;
template class LinearEllipticEocExpectations< LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              1 >;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
