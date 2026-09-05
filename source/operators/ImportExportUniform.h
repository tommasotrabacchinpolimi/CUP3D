//
//  Cubism3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Ivica Kicic (kicici@ethz.ch)
//

#pragma once

#include "../Definitions.h"
#include "../SimulationData.h"

CubismUP_3D_NAMESPACE_BEGIN

enum FieldExportId
{
  FE_CHI = 0,
  FE_U = 1,
  FE_V = 2,
  FE_W = 3,
  FE_P = 4,
  FE_TMPU = 5,
  FE_TMPV = 6,
  FE_TMPW = 7
};

/// Export one field component to a contiguous fully refined matrix.
void exportGridFieldToUniformMatrix(
    SimulationData &sim, int component, Real * __restrict__ out);

/// Import one field component from a contiguous fully refined matrix.
void importGridFieldFromUniformMatrix(
    SimulationData &sim, int component, const Real * __restrict__ in);

CubismUP_3D_NAMESPACE_END
