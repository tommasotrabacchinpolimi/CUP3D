//
//  Cubism3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Ivica Kicic (kicici@ethz.ch)
//

#include "ImportExportUniform.h"
#include <Cubism/ImportExport.hh>
#include <stdexcept>

CubismUP_3D_NAMESPACE_BEGIN

namespace {

void resolveField(SimulationData &sim, int component,
                  bool &isScalar, int &localComp,
                  ScalarGrid *&scalar, VectorGrid *&vector)
{
  scalar = nullptr;
  vector = nullptr;
  isScalar = true;
  localComp = 0;
  switch (component) {
    case FE_CHI:  scalar = sim.chi;  return;
    case FE_P:    scalar = sim.pres; return;
    case FE_U:    isScalar = false; localComp = 0; vector = sim.vel;  return;
    case FE_V:    isScalar = false; localComp = 1; vector = sim.vel;  return;
    case FE_W:    isScalar = false; localComp = 2; vector = sim.vel;  return;
    case FE_TMPU: isScalar = false; localComp = 0; vector = sim.tmpV; return;
    case FE_TMPV: isScalar = false; localComp = 1; vector = sim.tmpV; return;
    case FE_TMPW: isScalar = false; localComp = 2; vector = sim.tmpV; return;
    default:
      throw std::invalid_argument("unknown field component");
  }
}

}  // namespace

void exportGridFieldToUniformMatrix(SimulationData &sim, int component, Real *out)
{
  bool isScalar = true;
  int localComp = 0;
  ScalarGrid *scalar = nullptr;
  VectorGrid *vector = nullptr;
  resolveField(sim, component, isScalar, localComp, scalar, vector);

  if (isScalar) {
    cubism::exportGridToUniformMatrix<ScalarLab>(
        scalar,
        [](const ScalarElement &element) { return element.s; },
        {0},
        out);
  } else {
    cubism::exportGridToUniformMatrix<VectorLab>(
        vector,
        [localComp](const VectorElement &element) { return element.u[localComp]; },
        {localComp},
        out);
  }
}

void importGridFieldFromUniformMatrix(
    SimulationData &sim, int component, const Real * __restrict__ in)
{
  bool isScalar = true;
  int localComp = 0;
  ScalarGrid *scalar = nullptr;
  VectorGrid *vector = nullptr;
  resolveField(sim, component, isScalar, localComp, scalar, vector);

  if (isScalar) {
    cubism::importGridFromUniformMatrix(
        scalar,
        [](ScalarElement &element, Real value) { element.s = value; },
        in);
  } else {
    cubism::importGridFromUniformMatrix(
        vector,
        [localComp](VectorElement &element, Real value) { element.u[localComp] = value; },
        in);
  }
}

CubismUP_3D_NAMESPACE_END
