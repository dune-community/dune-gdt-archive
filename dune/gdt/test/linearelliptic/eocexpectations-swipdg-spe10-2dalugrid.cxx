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

# include "problems/spe10.hh"
# include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {

// polorder 1, conforming

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     1,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {7.22e-02, 2.59e-02};
    else if (type == "H1_semi")
      return {5.28e-01, 3.48e-01};
//    else if (type == "energy")
//      return {};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, conforming

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     2,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {2.08e-02, 3.74e-03};
    else if (type == "H1_semi")
      return {2.56e-01, 8.47e-02};
//    else if (type == "energy")
//      return {};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

//template< bool anything >
//class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
//                                     LinearElliptic::ChooseDiscretizer::swipdg,
//                                     1,
//                                     anything >
//  : public internal::LinearEllipticEocExpectationsBase< 1 >
//{
//  typedef LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 > TestCaseType;
//public:
//  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
//  {
//    if (type == "L2")
//      return {};
//    else if (type == "H1_semi")
//      return {};
//    else if (type == "energy")
//      return {};
//    else
//      EXPECT_TRUE(false) << "test results missing for type: " << type;
//    return {};
//  } // ... results(...)
//}; // LinearEllipticEocExpectations



template class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              1 >;
template class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              2 >;
//template class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
//                                              LinearElliptic::ChooseDiscretizer::swipdg,
//                                              1 >;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
