// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/sgrid.hh>

#include "problems/ER2007.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {


template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::ER2007TestCase< SGrid< 2, 2 >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     1,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::ER2007TestCase< SGrid< 2, 2 >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {9.46e-01, 9.67e-01};
    else if (type == "H1_semi" || type == "energy")
      return {1.65e+01, 2.22e+01};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations< LinearElliptic::ER2007TestCase< SGrid< 2, 2 >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              1 >;


} // namespace Test
} // namespace GDT
} // namespace Dune
