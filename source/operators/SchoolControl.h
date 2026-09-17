//
//  CubismUP_3D
//  Wang setpoints → hybrid plant: appliedForce surge + StefanFish::act({b}) yaw.
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
  // Turn(b) is a wave-shift event — call sparsely, hold b between calls.
  std::vector<Real> tLastAct;
  std::vector<Real> bHeld;
};

CubismUP_3D_NAMESPACE_END
