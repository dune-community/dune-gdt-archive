// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_SWIPDG_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_SWIPDG_HH

#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/la/container.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretizations/default.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localfunctional/codim0.hh>
#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/playground/localevaluation/swipdg.hh>
#include <dune/gdt/spaces/dg.hh>

#include "../problems/interface.hh"
#include "base.hh"

namespace Dune {
namespace GDT {
namespace LinearElliptic {


/**
 * \brief Discretizes a linear elliptic PDE using a continuous Galerkin Finite Element method.
 */
template< class GridType,
          Stuff::Grid::ChooseLayer layer = Stuff::Grid::ChooseLayer::leaf,
          ChooseSpaceBackend space_backend = Spaces::default_dg_backend,
          Stuff::LA::ChooseBackend la_backend = Stuff::LA::default_sparse_backend,
          int pol = 1,
          class RangeFieldType = double,
          size_t dimRange = 1 >
class SwipdgDiscretizer
{
public:
  typedef ProblemInterface< typename GridType::template Codim< 0 >::Entity,
                            typename GridType::ctype,
                            GridType::dimension,
                            RangeFieldType,
                            dimRange > ProblemType;
  typedef Spaces::DGProvider< GridType, layer, space_backend, pol, RangeFieldType, dimRange > SpaceProvider;
  typedef typename SpaceProvider::Type                                                        SpaceType;
  typedef typename Stuff::LA::Container< RangeFieldType, la_backend >::MatrixType             MatrixType;
  typedef typename Stuff::LA::Container< RangeFieldType, la_backend >::VectorType             VectorType;
  typedef Discretizations::StationaryContainerBasedDefault
      < ProblemType, SpaceType, MatrixType, VectorType, SpaceType >                           DiscretizationType;
  static const constexpr ChooseDiscretizer type = ChooseDiscretizer::swipdg;
  static const int polOrder = pol;

  static std::string static_id() //                                                        int() needed, otherwise we
  { //                                                                                     get a linker error
    return std::string("gdt.linearelliptic.discretization.swipdg.order_") + DSC::to_string(int(polOrder));
  }

  static DiscretizationType discretize(Stuff::Grid::ProviderInterface< GridType >& grid_provider,
                                       const ProblemType& problem,
                                       const int level = 0)
  {
    auto logger = Stuff::Common::TimedLogger().get(static_id());
    logger.info() << "Creating space... " << std::endl;
    auto space = SpaceProvider::create(grid_provider, level);
    logger.debug() << "grid has " << space.grid_view().indexSet().size(0) << " elements" << std::endl;
    typedef typename SpaceType::GridViewType    GridViewType;
    typedef typename GridViewType::Intersection IntersectionType;
    auto boundary_info = Stuff::Grid::BoundaryInfoProvider< IntersectionType >::create(problem.boundary_info_cfg());
    logger.info() << "Assembling... " << std::endl;
    VectorType rhs_vector(space.mapper().size(), 0.0);
    MatrixType system_matrix(space.mapper().size(),
                             space.mapper().size(),
                             space.compute_face_and_volume_pattern());
    typedef typename ProblemType::DiffusionFactorType DiffusionFactorType;
    typedef typename ProblemType::DiffusionTensorType DiffusionTensorType;
    typedef typename ProblemType::FunctionType FunctionType;
    // volume terms
    // * lhs
    typedef LocalOperator::Codim0Integral
        < LocalEvaluation::Elliptic< DiffusionFactorType, DiffusionTensorType > > EllipticOperatorType;
    const EllipticOperatorType                                  elliptic_operator(problem.diffusion_factor(),
                                                                                  problem.diffusion_tensor());
    const LocalAssembler::Codim0Matrix< EllipticOperatorType >  diffusion_matrix_assembler(elliptic_operator);
    // * rhs
    typedef LocalFunctional::Codim0Integral< LocalEvaluation::Product< FunctionType > > ForceFunctionalType;
    const ForceFunctionalType                                 force_functional(problem.force());
    const LocalAssembler::Codim0Vector< ForceFunctionalType > force_vector_assembler(force_functional);
    // inner face terms
    typedef LocalOperator::Codim1CouplingIntegral
        < LocalEvaluation::SWIPDG::Inner< DiffusionFactorType, DiffusionTensorType > > CouplingOperatorType;
    const CouplingOperatorType                                          coupling_operator(problem.diffusion_factor(),
                                                                                          problem.diffusion_tensor());
    const LocalAssembler::Codim1CouplingMatrix< CouplingOperatorType >  coupling_matrix_assembler(coupling_operator);
    // dirichlet boundary face terms
    // * lhs
    typedef LocalOperator::Codim1BoundaryIntegral
        < LocalEvaluation::SWIPDG::BoundaryLHS< DiffusionFactorType, DiffusionTensorType > > DirichletOperatorType;
    const DirichletOperatorType                                         dirichlet_operator(problem.diffusion_factor(),
                                                                                           problem.diffusion_tensor());
    const LocalAssembler::Codim1BoundaryMatrix< DirichletOperatorType > dirichlet_matrix_assembler(dirichlet_operator);
    // * rhs
    typedef LocalFunctional::Codim1Integral
        < LocalEvaluation::SWIPDG::BoundaryRHS< DiffusionFactorType, FunctionType, DiffusionTensorType > >
            DirichletFunctionalType;
    const DirichletFunctionalType                                 dirichlet_functional(problem.diffusion_factor(),
                                                                                       problem.diffusion_tensor(),
                                                                                       problem.dirichlet());
    const LocalAssembler::Codim1Vector< DirichletFunctionalType > dirichlet_vector_assembler(dirichlet_functional);
    // neumann boundary face terms
    // * rhs
    typedef LocalFunctional::Codim1Integral< LocalEvaluation::Product< FunctionType > > NeumannFunctionalType;
    const NeumannFunctionalType                                 neumann_functional(problem.neumann());
    const LocalAssembler::Codim1Vector< NeumannFunctionalType > neumann_vector_assembler(neumann_functional);
    // register everything for assembly in one grid walk
    SystemAssembler< SpaceType > assembler(space);
    assembler.add(diffusion_matrix_assembler, system_matrix);
    assembler.add(force_vector_assembler, rhs_vector);
    assembler.add(coupling_matrix_assembler,
                  system_matrix,
                  new Stuff::Grid::ApplyOn::InnerIntersectionsPrimally< GridViewType >());
    assembler.add(dirichlet_matrix_assembler,
                  system_matrix,
                  new Stuff::Grid::ApplyOn::DirichletIntersections< GridViewType >(*boundary_info));
    assembler.add(dirichlet_vector_assembler,
                  rhs_vector,
                  new Stuff::Grid::ApplyOn::DirichletIntersections< GridViewType >(*boundary_info));
    assembler.add(neumann_vector_assembler,
                  rhs_vector,
                  new Stuff::Grid::ApplyOn::NeumannIntersections< GridViewType >(*boundary_info));
    assembler.assemble();
    // create the discretization (no copy of the containers done here, bc. of cow)
    return DiscretizationType(problem, space, system_matrix, rhs_vector);
  } // ... discretize(...)
}; // class SwipdgDiscretizer


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_SWIPDG_HH
