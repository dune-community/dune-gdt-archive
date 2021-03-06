// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_SPACES_CG_FEM_HH
#define DUNE_GDT_TEST_SPACES_CG_FEM_HH

#include <dune/gdt/spaces/cg/fem.hh>

#include "grids.hh"

#if HAVE_DUNE_FEM

#define SPACE_CG_FEM_YASPGRID(dd, rr, pp) \
  Spaces::CG::FemBased< Yasp ## dd ## dLeafGridPartType, pp, double, rr >

#define SPACES_CG_FEM(pp) \
    SPACE_CG_FEM_YASPGRID(1, pp, 1) \
  , SPACE_CG_FEM_YASPGRID(2, 1, pp) \
  , SPACE_CG_FEM_YASPGRID(3, 1, pp)


# if HAVE_ALUGRID


#define SPACE_CG_FEM_ALUCONFORMGRID(dd, rr, pp) \
  Spaces::CG::FemBased< AluConform ## dd ## dLeafGridPartType, pp, double, rr >

#define SPACE_CG_FEM_ALUCUBEGRID(dd, rr, pp) \
  Spaces::CG::FemBased< AluCube ## dd ## dLeafGridPartType, pp, double, rr >

#define SPACES_CG_FEM_ALUGRID(pp) \
    SPACE_CG_FEM_ALUCONFORMGRID(2, 1, pp) \
  , SPACE_CG_FEM_ALUCONFORMGRID(3, 1, pp) \
  , SPACE_CG_FEM_ALUCUBEGRID(2, 1, pp) \
  , SPACE_CG_FEM_ALUCUBEGRID(3, 1, pp)


# endif // HAVE_ALUGRID

#endif // HAVE_DUNE_FEM

#endif // DUNE_GDT_TEST_SPACES_CG_FEM_HH
