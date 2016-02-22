// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_LOCALFLUXES_GODUNOV_HH
#define DUNE_GDT_LOCALFLUXES_GODUNOV_HH

#include <tuple>
#include <memory>

# if HAVE_EIGEN
#   include <Eigen/Eigenvalues>
# endif

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/affine.hh>
#include <dune/stuff/la/container/eigen.hh>
#include <dune/stuff/common/parallel/threadstorage.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forwards
template< class AnalyticalFluxImp, size_t domainDim >
class GodunovNumericalCouplingFlux;

template< class AnalyticalFluxImp, class BoundaryValueFunctionType, size_t domainDim >
class GodunovNumericalBoundaryFlux;


namespace internal {

/**
 *  \brief  Traits for the Lax-Friedrichs flux evaluation.
 */
template< class AnalyticalCouplingFluxImp, size_t domainDim = AnalyticalCouplingFluxImp::dimDomain >
class GodunovNumericalCouplingFluxTraits
{
  static_assert(std::is_base_of< Dune::GDT::IsAnalyticalFlux, AnalyticalFluxImp >::value,
                "AnalyticalFluxImp has to be derived from GDT::IsAnalyticalFlux.");
public:
  typedef AnalyticalCouplingFluxImp                            AnalyticalFluxType;
  typedef GodunovNumericalCouplingFlux< AnalyticalFluxType, domainDim >               derived_type;
  typedef typename AnalyticalFluxType::EntityType              EntityType;
  typedef typename AnalyticalFluxType::DomainFieldType         DomainFieldType;
  typedef typename AnalyticalFluxType::RangeFieldType          RangeFieldType;
  typedef std::tuple< double >                                 LocalfunctionTupleType;
  static const size_t dimDomain = domainDim;
  static const size_t dimRange = AnalyticalFluxType::dimRange;
  static_assert(AnalyticalFluxType::dimRangeCols == 1, "Not implemented for dimRangeCols > 1!");

  typedef typename AnalyticalFluxType::RangeType                                        RangeType;
  typedef typename AnalyticalFluxType::FluxRangeType                                    FluxRangeType;
  typedef typename AnalyticalFluxType::FluxJacobianRangeType                            FluxJacobianRangeType;
  typedef typename Dune::Stuff::LA::EigenDenseMatrix< RangeFieldType >                  EigenMatrixType;
}; // class GodunovNumericalCouplingFluxTraits

template< class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp, size_t domainDim = AnalyticalBoundaryFluxImp::dimDomain >
class GodunovNumericalBoundaryFluxTraits
   : public GodunovNumericalCouplingFluxTraits< AnalyticalBoundaryFluxImp, domainDim >
{
  typedef GodunovNumericalCouplingFluxTraits< AnalyticalBoundaryFluxImp, domainDim > BaseType;
public:
  typedef AnalyticalBoundaryFluxImp                                         AnalyticalFluxType;
  typedef BoundaryValueFunctionImp                                          BoundaryValueFunctionType;
  typedef typename BoundaryValueFunctionType::LocalfunctionType             LocalfunctionType;
  typedef GodunovNumericalBoundaryFlux< AnalyticalFluxType, BoundaryValueFunctionType, domainDim >   derived_type;
  typedef std::tuple< std::shared_ptr< LocalfunctionType > >                LocalfunctionTupleType;
}; // class GodunovNumericalBoundaryFluxTraits


} // namespace internal



#define PAPERFLUX 0

template< class AnalyticalCouplingFluxImp, size_t domainDim = AnalyticalCouplingFluxImp::dimDomain >
class GodunovNumericalCouplingFlux
  : public NumericalCouplingFluxInterface< internal::GodunovNumericalCouplingFluxTraits< AnalyticalCouplingFluxImp > >
{
public:
  typedef internal::GodunovNumericalCouplingFluxTraits< AnalyticalCouplingFluxImp, domainDim >           Traits;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType                    FluxJacobianRangeType;
  typedef typename Traits::RangeType                                RangeType;
  typedef typename Traits::EigenMatrixType                          EigenMatrixType;
  static constexpr size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit GodunovNumericalCouplingFlux(const AnalyticalFluxType& analytical_flux,
                 const bool is_linear = false)
    : analytical_flux_(analytical_flux)
    , is_linear_(is_linear)
  {
    if (!jacobians_constructed_ && is_linear_)
      initialize_jacobians();
  }

  LocalfunctionTupleType local_functions(const EntityType& /*entity*/) const
  {
    return LocalfunctionTupleType();
  }

  template< class IntersectionType >
  ResultType evaluate(const LocalfunctionTupleType& /*local_functions_tuple_entity*/,
                      const LocalfunctionTupleType& /*local_functions_tuple_neighbor*/,
                      const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& local_source_entity,
                      const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& local_source_neighbor,
                      const IntersectionType& intersection,
                      const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& x_intersection) const
  {
    // get function values
    const RangeType u_i = local_source_entity.evaluate(intersection.geometryInInside().global(x_intersection));
    const RangeType u_j = local_source_neighbor.evaluate(intersection.geometryInOutside().global(x_intersection));
    // get flux value
    const FluxRangeType f_u_i = analytical_flux_.evaluate(u_i);

    if (!is_linear_) { // use simple linearized Riemann solver, LeVeque p.316
      reinitialize_jacobians(u_i, u_j);
    }
    // get jump at the intersection
    const RangeType delta_u = u_i - u_j;
    // get unit outer normal
    const auto n_ij = intersection.unitOuterNormal(x_intersection);
    // find direction of unit outer normal
    size_t coord = 0;
#ifndef NDEBUG
    size_t num_zeros = 0;
#endif //NDEBUG
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(1)) || DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(-1)))
        coord = ii;
      else if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(0))) {
#ifndef NDEBUG
        ++num_zeros;
#endif //NDEBUG
      }
      else
        DUNE_THROW(Dune::NotImplemented, "Godunov flux is only implemented for axis parallel cube grids");
    }
    assert(num_zeros == dimDomain - 1);
    // calculate return vector
    ResultType ret;
    RangeFieldType vol_intersection = intersection.geometry().volume();
    if (n_ij[coord] > 0) {
      RangeType negative_waves(RangeFieldType(0));
      jacobian_neg_[coord].mv(delta_u, negative_waves);
      for (size_t kk = 0; kk < dimRange; ++kk)
        ret[0][kk] = (f_u_i[kk][coord] - negative_waves[kk]*n_ij[coord])*vol_intersection;
    } else {
      RangeType positive_waves(RangeFieldType(0));
      jacobian_pos_[coord].mv(delta_u, positive_waves);
      for (size_t kk = 0; kk < dimRange; ++kk)
        ret[0][kk] = (-f_u_i[kk][coord] - positive_waves[kk]*n_ij[coord])*vol_intersection;
    }
    return ret;
  } // ResultType evaluate(...) const

private:
  void initialize_jacobians()
  {
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(RangeType(0)));
    std::vector< EigenMatrixType > jacobian_eigen;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      jacobian_eigen.emplace_back(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian[ii], 15)));
    calculate_jacobians(std::move(jacobian_eigen));
    jacobians_constructed_ = true;
  } // void initialize_jacobians()

  void reinitialize_jacobians(const RangeType& u_i,
                              const RangeType& u_j) const
  {
    // calculate jacobian as jacobian(0.5*(u_i+u_j)
    RangeType u_mean = u_i + u_j;
    u_mean *= RangeFieldType(0.5);
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(u_mean));
    std::vector< EigenMatrixType > jacobian_eigen;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      jacobian_eigen.emplace_back(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian[ii], 15)));
    // calculate jacobian_neg and jacobian_pos
    calculate_jacobians(std::move(jacobian_eigen));
  } // void reinitialize_jacobians(...)

  void calculate_jacobians(std::vector< EigenMatrixType >&& jacobian) const
  {
#if HAVE_EIGEN
    EigenMatrixType diag_jacobian_pos_tmp(dimRange, dimRange, RangeFieldType(0));
    EigenMatrixType diag_jacobian_neg_tmp(dimRange, dimRange, RangeFieldType(0));
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      diag_jacobian_pos_tmp.scal(RangeFieldType(0));
      diag_jacobian_neg_tmp.scal(RangeFieldType(0));
      // create EigenSolver
      ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
          eigen_solver(jacobian[ii].backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto eigenvalues = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
      const auto eigenvectors = eigen_solver.eigenvectors(); // <- this should be an Eigen matrix of std::complex
      assert(boost::numeric_cast< size_t >(eigenvalues.size()) == dimRange);
      for (size_t jj = 0; jj < dimRange; ++jj) {
        // assert this is real
        assert(std::abs(eigenvalues[jj].imag()) < 1e-15);
        const RangeFieldType eigenvalue = eigenvalues[jj].real();
        if (eigenvalue < 0)
          diag_jacobian_neg_tmp.set_entry(jj, jj, eigenvalue);
        else
          diag_jacobian_pos_tmp.set_entry(jj, jj, eigenvalue);
      }
      const auto eigenvectors_inverse = eigenvectors.inverse();
      EigenMatrixType jacobian_neg_eigen(eigenvectors.real()*diag_jacobian_neg_tmp.backend()*eigenvectors_inverse.real());
      EigenMatrixType jacobian_pos_eigen(eigenvectors.real()*diag_jacobian_pos_tmp.backend()*eigenvectors_inverse.real());
      jacobian_neg_[ii] = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_neg_eigen, 15));
      jacobian_pos_[ii] = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_pos_eigen, 15));
    }
#else
    static_assert(AlwaysFalse< FluxJacobianRangeType >::value, "You are missing eigen!");
#endif
  } // void calculate_jacobians(...)

  const AnalyticalFluxType& analytical_flux_;
  static FluxJacobianRangeType jacobian_neg_;
  static FluxJacobianRangeType jacobian_pos_;
  static bool jacobians_constructed_;
  const bool is_linear_;
}; // class GodunovNumericalCouplingFlux

template < class AnalyticalCouplingFluxImp, size_t dimDomain >
typename internal::GodunovNumericalCouplingFluxTraits< AnalyticalCouplingFluxImp, dimDomain >::FluxJacobianRangeType
GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, dimDomain >::jacobian_neg_(0);

template < class AnalyticalCouplingFluxImp, size_t dimDomain >
typename internal::GodunovNumericalCouplingFluxTraits< AnalyticalCouplingFluxImp, dimDomain >::FluxJacobianRangeType
GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, dimDomain >::jacobian_pos_(0);

template < class AnalyticalCouplingFluxImp, size_t dimDomain >
bool
GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, dimDomain >::jacobians_constructed_(false);

template< class AnalyticalCouplingFluxImp >
class GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >
  : public NumericalCouplingFluxInterface< internal::GodunovNumericalCouplingFluxTraits< AnalyticalCouplingFluxImp > >
{
public:
  typedef internal::GodunovNumericalCouplingFluxTraits< AnalyticalCouplingFluxImp, 1 >        Traits;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType                    FluxJacobianRangeType;
  typedef typename Traits::RangeType                                RangeType;
  typedef typename Traits::EigenMatrixType                          EigenMatrixType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  typedef typename DS::Functions::Affine< typename AnalyticalFluxType::EntityType,
                                          RangeFieldType, dimRange,
                                          RangeFieldType, dimRange, 1 >         AffineFunctionType;

  explicit GodunovNumericalCouplingFlux(const AnalyticalFluxType& analytical_flux,
                 const bool is_linear = false)
    : analytical_flux_(analytical_flux)
    , is_linear_(is_linear)
  {
    if (!jacobians_constructed_ && is_linear_)
      initialize_jacobians();
  }

  LocalfunctionTupleType local_functions(const EntityType& /*entity*/) const
  {
    return LocalfunctionTupleType();
  }

  template< class IntersectionType >
  ResultType evaluate(const LocalfunctionTupleType& /*local_functions_tuple_entity*/,
                      const LocalfunctionTupleType& /*local_functions_tuple_neighbor*/,
                      const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& local_source_entity,
                      const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& local_source_neighbor,
                      const IntersectionType& intersection,
                      const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& x_intersection) const
  {
    // get function values
    const RangeType u_i = local_source_entity.evaluate(intersection.geometryInInside().global(x_intersection));
    const RangeType u_j = local_source_neighbor.evaluate(intersection.geometryInOutside().global(x_intersection));

    if (!jacobians_constructed_ && is_linear_)
      initialize_jacobians();

    if (!is_linear_) { // use simple linearized Riemann solver, LeVeque p.316
      reinitialize_jacobians(u_i, u_j);
    }
    // get unit outer normal
    const RangeFieldType n_ij = intersection.unitOuterNormal(x_intersection)[0];
    // calculate return vector
    ResultType ret;
#if PAPERFLUX
    //get flux values, FluxRangeType should be FieldVector< ..., dimRange >
    const FluxRangeType f_u_i_plus_f_u_j = is_linear_
                                           ? analytical_flux_.evaluate(u_i + u_j)
                                           : analytical_flux_.evaluate(u_i) + analytical_flux_.evaluate(u_j);
    RangeType waves(RangeFieldType(0));
    if (n_ij > 0) {
      // flux = 0.5*(f_u_i + f_u_j + |A|*(u_i-u_j))*n_ij
      jacobian_abs_function_.evaluate(u_i - u_j, waves);
      ret.axpy(RangeFieldType(0.5), f_u_i_plus_f_u_j + waves);
    } else {
      // flux = 0.5*(f_u_i + f_u_j - |A|*(u_i-u_j))*n_ij
      jacobian_abs_function_.evaluate(u_j - u_i, waves);
      ret.axpy(RangeFieldType(-0.5), f_u_i_plus_f_u_j + waves);
    }
#else
    const FluxRangeType f_u_i = analytical_flux_.evaluate(u_i);
    if (n_ij > 0) {
      RangeType negative_waves(RangeFieldType(0));
      jacobian_neg_function_.evaluate(u_j - u_i, negative_waves);
      ret = Dune::DynamicVector< RangeFieldType >(f_u_i + negative_waves);
    } else {
      RangeType positive_waves(RangeFieldType(0));
      jacobian_pos_function_.evaluate(u_i - u_j, positive_waves);
      ret = Dune::DynamicVector< RangeFieldType >(positive_waves - f_u_i);
    }
#endif
  } // void evaluate(...) const

private:
  void initialize_jacobians() const
  {
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(RangeType(0)));
    EigenMatrixType jacobian_eigen(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian, 15)));
    calculate_jacobians(std::move(jacobian_eigen));
    jacobians_constructed_ = true;
  } // void initialize_jacobians()

  void reinitialize_jacobians(const RangeType& u_i,
                              const RangeType& u_j) const
  {
    // calculate jacobian as jacobian(0.5*(u_i+u_j)
    RangeType u_mean = u_i + u_j;
    u_mean *= RangeFieldType(0.5);
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(u_mean));
    EigenMatrixType jacobian_eigen = DSC::fromString< EigenMatrixType >(DSC::toString(jacobian, 15));

    // calculate jacobian_neg and jacobian_pos
    calculate_jacobians(std::move(jacobian_eigen));
  } // void reinitialize_jacobians(...)

  void calculate_jacobians(EigenMatrixType&& jacobian) const
  {
#if HAVE_EIGEN
    EigenMatrixType diag_jacobian_pos_tmp(dimRange, dimRange, RangeFieldType(0));
    EigenMatrixType diag_jacobian_neg_tmp(dimRange, dimRange, RangeFieldType(0));
    diag_jacobian_pos_tmp.scal(RangeFieldType(0));
    diag_jacobian_neg_tmp.scal(RangeFieldType(0));
    // create EigenSolver
    ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
        eigen_solver(jacobian.backend());
    assert(eigen_solver.info() == ::Eigen::Success);
    const auto eigenvalues = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
    const auto eigenvectors = eigen_solver.eigenvectors(); // <- this should be an Eigen matrix of std::complex
    assert(boost::numeric_cast< size_t >(eigenvalues.size()) == dimRange);
    for (size_t jj = 0; jj < dimRange; ++jj) {
      // assert this is real
      assert(std::abs(eigenvalues[jj].imag()) < 1e-15);
      const RangeFieldType eigenvalue = eigenvalues[jj].real();
      if (eigenvalue < 0)
        diag_jacobian_neg_tmp.set_entry(jj, jj, eigenvalue);
      else
        diag_jacobian_pos_tmp.set_entry(jj, jj, eigenvalue);
      const auto eigenvectors_inverse = eigenvectors.inverse();
      EigenMatrixType jacobian_neg_eigen(eigenvectors.real()*diag_jacobian_neg_tmp.backend()*eigenvectors_inverse.real());
      EigenMatrixType jacobian_pos_eigen(eigenvectors.real()*diag_jacobian_pos_tmp.backend()*eigenvectors_inverse.real());
      // set jacobian_neg_ and jacobian_pos_
      jacobian_neg_ = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_neg_eigen, 15));
      jacobian_pos_ = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_pos_eigen, 15));
      // jacobian_abs_ = jacobian_pos_ - jacobian_neg_;
      jacobian_abs_ = jacobian_neg_;
      jacobian_abs_ *= RangeFieldType(-1.0);
      jacobian_abs_ += jacobian_pos_;
# if PAPERFLUX
      jacobian_abs_function_ = AffineFunctionType(jacobian_abs_, RangeType(0), true);
# else
      jacobian_neg_function_ = AffineFunctionType(jacobian_neg_, RangeType(0), true);
      jacobian_pos_function_ = AffineFunctionType(jacobian_pos_, RangeType(0), true);
# endif
    }
#else
    static_assert(AlwaysFalse< FluxJacobianRangeType >::value, "You are missing eigen!");
#endif
  } // void calculate_jacobians(...)

  const AnalyticalFluxType& analytical_flux_;
  thread_local static FluxJacobianRangeType jacobian_neg_;
  thread_local static FluxJacobianRangeType jacobian_pos_;
  thread_local static FluxJacobianRangeType jacobian_abs_;
  thread_local static AffineFunctionType jacobian_neg_function_;
  thread_local static AffineFunctionType jacobian_pos_function_;
  thread_local static AffineFunctionType jacobian_abs_function_;
  thread_local static bool jacobians_constructed_;
  const bool is_linear_;
}; // class GodunovNumericalCouplingFlux< ..., 1 >

template < class AnalyticalCouplingFluxImp >
thread_local typename GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::FluxJacobianRangeType
GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::jacobian_neg_{typename GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::FluxJacobianRangeType()};

template < class AnalyticalCouplingFluxImp >
thread_local typename GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::FluxJacobianRangeType
GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::jacobian_pos_{typename GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::FluxJacobianRangeType()};

template < class AnalyticalCouplingFluxImp >
thread_local typename GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::FluxJacobianRangeType
GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::jacobian_abs_{typename GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::FluxJacobianRangeType()};

template < class AnalyticalCouplingFluxImp >
thread_local typename GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::AffineFunctionType
GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::jacobian_neg_function_(GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::AffineFunctionType(GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::FluxJacobianRangeType(0)));

template < class AnalyticalCouplingFluxImp >
thread_local typename GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::AffineFunctionType
GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::jacobian_pos_function_(GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::AffineFunctionType(GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::FluxJacobianRangeType(0)));

template < class AnalyticalCouplingFluxImp >
thread_local typename GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::AffineFunctionType
GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::jacobian_abs_function_(GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::AffineFunctionType(GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::FluxJacobianRangeType(0)));

template < class AnalyticalCouplingFluxImp >
thread_local bool
GodunovNumericalCouplingFlux< AnalyticalCouplingFluxImp, 1 >::jacobians_constructed_(false);

template< class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp, size_t domainDim = AnalyticalBoundaryFluxImp::dimDomain >
class GodunovNumericalBoundaryFlux
  : public NumericalBoundaryFluxInterface
                          < internal::GodunovNumericalBoundaryFluxTraits< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp > >
{
public:
  typedef internal::GodunovNumericalBoundaryFluxTraits< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, domainDim >  Traits;
  typedef typename Traits::BoundaryValueFunctionType                BoundaryValueFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType                    FluxJacobianRangeType;
  typedef typename Traits::RangeType                                RangeType;
  typedef typename Traits::EigenMatrixType                          EigenMatrixType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit GodunovNumericalBoundaryFlux(const AnalyticalFluxType& analytical_flux,
                     const BoundaryValueFunctionType& boundary_values,
                     const bool is_linear = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , is_linear_(is_linear)
  {
    if (!jacobians_constructed_ && is_linear_)
      initialize_jacobians();
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(boundary_values_.local_function(entity));
  }


  template< class IntersectionType >
  ResultType evaluate(
      const LocalfunctionTupleType& local_functions_tuple,
      const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& local_source_entity,
      const IntersectionType& intersection,
      const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& x_intersection) const
  {
      const auto x_intersection_entity_coords = intersection.geometryInInside().global(x_intersection);
      const RangeType u_i = local_source_entity.evaluate(intersection.geometryInInside().global(x_intersection_entity_coords));
      const RangeType u_j = std::get< 0 >(local_functions_tuple)->evaluate(x_intersection_entity_coords);
      // get flux values
      const FluxRangeType f_u_i = analytical_flux_.evaluate(u_i);

      if (!is_linear_) { // use simple linearized Riemann solver, LeVeque p.316
        reinitialize_jacobians(u_i, u_j);
      }
      // get jump at the intersection
      const RangeType delta_u = u_i - u_j;
      // get unit outer normal
      const auto n_ij = intersection.unitOuterNormal(localPoint);
      // find direction of unit outer normal
      size_t coord = 0;
  #ifndef NDEBUG
      size_t num_zeros = 0;
  #endif // NDEBUG
      for (size_t ii = 0; ii < dimDomain; ++ii) {
        if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(1)) || DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(-1)))
          coord = ii;
        else if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(0))) {
  #ifndef NDEBUG
          ++num_zeros;
  #endif // NDEBUG
        }
        else
          DUNE_THROW(Dune::NotImplemented, "Godunov flux is only implemented for axis parallel cube grids");
      }
      assert(num_zeros == dimDomain - 1);
      // calculate return vector
      ReturnType ret;
      RangeFieldType vol_intersection = intersection.geometry().volume();
      if (n_ij[coord] > 0) {
        RangeType negative_waves(RangeFieldType(0));
        jacobian_neg_[coord].mv(delta_u, negative_waves);
        for (size_t kk = 0; kk < dimRange; ++kk)
          ret[kk] = (f_u_i[kk][coord] - negative_waves[kk]*n_ij[coord])*vol_intersection;
      } else {
        RangeType positive_waves(RangeFieldType(0));
        jacobian_pos_[coord].mv(delta_u, positive_waves);
        for (size_t kk = 0; kk < dimRange; ++kk)
          ret[kk] = (-f_u_i[kk][coord] - positive_waves[kk]*n_ij[coord])*vol_intersection;
      }
    } // void evaluate(...) const

  private:
    void initialize_jacobians()
    {
      const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(RangeType(0)));
      std::vector< EigenMatrixType > jacobian_eigen;
      for (size_t ii = 0; ii < dimDomain; ++ii)
        jacobian_eigen.emplace_back(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian[ii])));
      calculate_jacobians(std::move(jacobian_eigen));
      jacobians_constructed_ = true;
    } // void initialize_jacobians()

    void reinitialize_jacobians(const RangeType& u_i,
                                const RangeType& u_j) const
    {
      // calculate jacobian as jacobian(0.5*(u_i+u_j)
      RangeType u_mean = u_i + u_j;
      u_mean *= RangeFieldType(0.5);
      const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(u_mean));
      std::vector< EigenMatrixType > jacobian_eigen;
      for (size_t ii = 0; ii < dimDomain; ++ii)
        jacobian_eigen.emplace_back(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian[ii])));
      // calculate jacobian_neg and jacobian_pos
      calculate_jacobians(std::move(jacobian_eigen));
    } // void reinitialize_jacobians(...)

    void calculate_jacobians(std::vector< EigenMatrixType >&& jacobian) const
    {
#if HAVE_EIGEN
      EigenMatrixType diag_jacobian_pos_tmp(dimRange, dimRange, RangeFieldType(0));
      EigenMatrixType diag_jacobian_neg_tmp(dimRange, dimRange, RangeFieldType(0));
      for (size_t ii = 0; ii < dimDomain; ++ii) {
        diag_jacobian_pos_tmp.scal(RangeFieldType(0));
        diag_jacobian_neg_tmp.scal(RangeFieldType(0));
        // create EigenSolver
        ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
            eigen_solver(jacobian[ii].backend());
        assert(eigen_solver.info() == ::Eigen::Success);
        const auto eigenvalues = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
        const auto eigenvectors = eigen_solver.eigenvectors(); // <- this should be an Eigen matrix of std::complex
        assert(boost::numeric_cast< size_t >(eigenvalues.size()) == dimRange);
        for (size_t jj = 0; jj < dimRange; ++jj) {
          // assert this is real
          assert(std::abs(eigenvalues[jj].imag()) < 1e-15);
          const RangeFieldType eigenvalue = eigenvalues[jj].real();
          if (eigenvalue < 0)
            diag_jacobian_neg_tmp.set_entry(jj, jj, eigenvalue);
          else
            diag_jacobian_pos_tmp.set_entry(jj, jj, eigenvalue);
        }
        const auto eigenvectors_inverse = eigenvectors.inverse();
        EigenMatrixType jacobian_neg_eigen(eigenvectors.real()*diag_jacobian_neg_tmp.backend()*eigenvectors_inverse.real());
        EigenMatrixType jacobian_pos_eigen(eigenvectors.real()*diag_jacobian_pos_tmp.backend()*eigenvectors_inverse.real());
        jacobian_neg_[ii] = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_neg_eigen));
        jacobian_pos_[ii] = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_pos_eigen));
      }
#else
      static_assert(AlwaysFalse< FluxJacobianRangeType >::value, "You are missing eigen!");
#endif
    } // void calculate_jacobians(...)

  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  static FluxJacobianRangeType jacobian_neg_;
  static FluxJacobianRangeType jacobian_pos_;
  static bool jacobians_constructed_;
  const bool is_linear_;
}; // class GodunovNumericalBoundaryFlux

template < class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp, size_t domainDim >
typename internal::GodunovNumericalBoundaryFluxTraits< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, domainDim >::FluxJacobianRangeType
GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, domainDim >::jacobian_neg_(0);

template < class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp, size_t domainDim >
typename internal::GodunovNumericalBoundaryFluxTraits< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, domainDim >::FluxJacobianRangeType
GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, domainDim >::jacobian_pos_(0);

template < class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp, size_t domainDim >
bool
GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, domainDim >::jacobians_constructed_(false);

template< class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp >
class GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >
  : public NumericalBoundaryFluxInterface< internal::GodunovNumericalBoundaryFluxTraits< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp > >
{
public:
  typedef internal::GodunovNumericalBoundaryFluxTraits< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >  Traits;
  typedef typename Traits::BoundaryValueFunctionType                BoundaryValueFunctionType;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType                    FluxJacobianRangeType;
  typedef typename Traits::RangeType                                RangeType;
  typedef typename Traits::EigenMatrixType                          EigenMatrixType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  typedef typename DS::Functions::Affine< typename AnalyticalFluxType::EntityType,
                                          RangeFieldType, dimRange,
                                          RangeFieldType, dimRange, 1 >         AffineFunctionType;
  explicit GodunovNumericalBoundaryFlux(const AnalyticalFluxType& analytical_flux,
                     const BoundaryValueFunctionType& boundary_values,
                     const bool is_linear = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , is_linear_(is_linear)
  {
    if (!jacobians_constructed_ && is_linear_)
      initialize_jacobians();
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(boundary_values_.local_function(entity));
  }

  template< class IntersectionType >
  ResultType evaluate(
      const LocalfunctionTupleType& local_functions_tuple,
      const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& local_source_entity,
      const IntersectionType& intersection,
      const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& x_intersection) const
  {
      const auto x_intersection_entity_coords = intersection.geometryInInside().global(x_intersection);
      const RangeType u_i = local_source_entity.evaluate(intersection.geometryInInside().global(x_intersection_entity_coords));
      const RangeType u_j = std::get< 0 >(local_functions_tuple)->evaluate(x_intersection_entity_coords);

      if (!jacobians_constructed_ && is_linear_)
        initialize_jacobians();

      if (!is_linear_) { // use simple linearized Riemann solver, LeVeque p.316
        reinitialize_jacobians(u_i, u_j);
      }
      // get unit outer normal
      const RangeFieldType n_ij = intersection.unitOuterNormal(localPoint)[0];
      // calculate return vector
      ReturnType ret;
#if PAPERFLUX
      //get flux values, FluxRangeType should be FieldVector< ..., dimRange >
      const FluxRangeType f_u_i_plus_f_u_j = is_linear_
                                             ? analytical_flux_.evaluate(u_i + u_j)
                                             : analytical_flux_.evaluate(u_i) + analytical_flux_.evaluate(u_j);
      RangeType waves(RangeFieldType(0));
      if (n_ij > 0) {
        // flux = 0.5*(f_u_i + f_u_j + |A|*(u_i-u_j))*n_ij
        jacobian_abs_function_.evaluate(u_i - u_j, waves);
        ret[0].axpy(RangeFieldType(0.5), f_u_i_plus_f_u_j + waves);
      } else {
        // flux = 0.5*(f_u_i + f_u_j - |A|*(u_i-u_j))*n_ij
        jacobian_abs_function_.evaluate(u_j - u_i, waves);
        ret[0].axpy(RangeFieldType(-0.5), f_u_i_plus_f_u_j + waves);
      }
#else
      const FluxRangeType f_u_i = analytical_flux_.evaluate(u_i);
      if (n_ij > 0) {
        RangeType negative_waves(RangeFieldType(0));
        jacobian_neg_function_.evaluate(u_j - u_i, negative_waves);
        ret[0] = Dune::DynamicVector< RangeFieldType >(f_u_i + negative_waves);
      } else {
        RangeType positive_waves(RangeFieldType(0));
        jacobian_pos_function_.evaluate(u_i - u_j, positive_waves);
        ret[0] = Dune::DynamicVector< RangeFieldType >(positive_waves - f_u_i);
      }
#endif
  } // void evaluate(...) const

  private:
    void initialize_jacobians() const
    {
      const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(RangeType(0)));
      EigenMatrixType jacobian_eigen(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian, 15)));
      calculate_jacobians(std::move(jacobian_eigen));
      jacobians_constructed_ = true;
    } // void initialize_jacobians()

    void reinitialize_jacobians(const RangeType& u_i,
                                const RangeType& u_j) const
    {
      // calculate jacobian as jacobian(0.5*(u_i+u_j)
      RangeType u_mean = u_i + u_j;
      u_mean *= RangeFieldType(0.5);
      const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(u_mean));
      EigenMatrixType jacobian_eigen = DSC::fromString< EigenMatrixType >(DSC::toString(jacobian, 15));

      // calculate jacobian_neg and jacobian_pos
      calculate_jacobians(std::move(jacobian_eigen));
    } // void reinitialize_jacobians(...)

    void calculate_jacobians(EigenMatrixType&& jacobian) const
    {
#if HAVE_EIGEN
      EigenMatrixType diag_jacobian_pos_tmp(dimRange, dimRange, RangeFieldType(0));
      EigenMatrixType diag_jacobian_neg_tmp(dimRange, dimRange, RangeFieldType(0));
      diag_jacobian_pos_tmp.scal(RangeFieldType(0));
      diag_jacobian_neg_tmp.scal(RangeFieldType(0));
      // create EigenSolver
      ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
          eigen_solver(jacobian.backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto eigenvalues = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
      const auto eigenvectors = eigen_solver.eigenvectors(); // <- this should be an Eigen matrix of std::complex
      assert(boost::numeric_cast< size_t >(eigenvalues.size()) == dimRange);
      for (size_t jj = 0; jj < dimRange; ++jj) {
        // assert this is real
        assert(std::abs(eigenvalues[jj].imag()) < 1e-15);
        const RangeFieldType eigenvalue = eigenvalues[jj].real();
        if (eigenvalue < 0)
          diag_jacobian_neg_tmp.set_entry(jj, jj, eigenvalue);
        else
          diag_jacobian_pos_tmp.set_entry(jj, jj, eigenvalue);
        const auto eigenvectors_inverse = eigenvectors.inverse();
        EigenMatrixType jacobian_neg_eigen(eigenvectors.real()*diag_jacobian_neg_tmp.backend()*eigenvectors_inverse.real());
        EigenMatrixType jacobian_pos_eigen(eigenvectors.real()*diag_jacobian_pos_tmp.backend()*eigenvectors_inverse.real());
        // set jacobian_neg_ and jacobian_pos_
        jacobian_neg_ = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_neg_eigen, 15));
        jacobian_pos_ = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_pos_eigen, 15));
        // jacobian_abs_ = jacobian_pos_ - jacobian_neg_;
        jacobian_abs_ = jacobian_neg_;
        jacobian_abs_ *= RangeFieldType(-1.0);
        jacobian_abs_ += jacobian_pos_;
# if PAPERFLUX
        jacobian_abs_function_ = AffineFunctionType(jacobian_abs_, RangeType(0), true);
# else
        jacobian_neg_function_ = AffineFunctionType(jacobian_neg_, RangeType(0), true);
        jacobian_pos_function_ = AffineFunctionType(jacobian_pos_, RangeType(0), true);
# endif
      }
#else
      static_assert(AlwaysFalse< FluxJacobianRangeType >::value, "You are missing eigen!");
#endif
    } // void calculate_jacobians(...)

    const AnalyticalFluxType& analytical_flux_;
    const BoundaryValueFunctionType& boundary_values_;
    thread_local static FluxJacobianRangeType jacobian_neg_;
    thread_local static FluxJacobianRangeType jacobian_pos_;
    thread_local static FluxJacobianRangeType jacobian_abs_;
    thread_local static AffineFunctionType jacobian_neg_function_;
    thread_local static AffineFunctionType jacobian_pos_function_;
    thread_local static AffineFunctionType jacobian_abs_function_;
    thread_local static bool jacobians_constructed_;
    const bool is_linear_;
  }; // class GodunovNumericalBoundaryFlux< ..., 1 >

template < class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp >
thread_local typename GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::FluxJacobianRangeType
GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::jacobian_neg_{typename GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::FluxJacobianRangeType()};

template < class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp >
thread_local typename GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::FluxJacobianRangeType
GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::jacobian_pos_{typename GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::FluxJacobianRangeType()};

template < class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp >
thread_local typename GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::FluxJacobianRangeType
GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::jacobian_abs_{typename GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::FluxJacobianRangeType()};

template < class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp >
thread_local typename GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::AffineFunctionType
GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::jacobian_neg_function_(GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::AffineFunctionType(GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::FluxJacobianRangeType(0)));

template < class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp >
thread_local typename GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::AffineFunctionType
GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::jacobian_pos_function_(GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::AffineFunctionType(GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::FluxJacobianRangeType(0)));

template < class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp >
thread_local typename GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::AffineFunctionType
GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::jacobian_abs_function_(GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::AffineFunctionType(GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::FluxJacobianRangeType(0)));

template < class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp >
thread_local bool
GodunovNumericalBoundaryFlux< AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1 >::jacobians_constructed_(false);


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALFLUXES_GODUNOV_HH