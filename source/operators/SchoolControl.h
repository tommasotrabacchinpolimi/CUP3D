//
//  CubismUP_3D
//  Wang setpoints → StefanFish kinematics: Turn(b) + period(a) on half-period clock.
//

#pragma once

#include "Operator.h"
#include "../Utils/WangSchoolPID.h"

#include <fstream>
#include <memory>
#include <vector>

CubismUP_3D_NAMESPACE_BEGIN

class StefanFish;

class SchoolControl : public Operator
{
public:
  SchoolControl(SimulationData & s);

  void operator()(const Real dt) override;
  std::string getName() override { return "SchoolControl"; }

private:
  WangSchoolPID controller;
  bool lengthScaleFromFish = true;
  bool warnedCorrectPosition = false;
  std::ofstream log;
  int logEvery = 20;
  // Turn(b) and period(a) are wave events — call sparsely, hold between calls.
  std::vector<Real> tLastAct;
  std::vector<Real> bHeld;
  std::vector<Real> aHeld;
};

CubismUP_3D_NAMESPACE_END
