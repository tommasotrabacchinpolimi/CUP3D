//
//  CubismUP_3D
//  School control via StefanFish kinematics only (no appliedForce/Torque):
//    yaw    → Turn(b)
//    surge  → period action a  (a<0 → shorter T → higher speed)
//  Actions applied on the 0.5*Tperiod clock (not every CFD step).
//

#include "SchoolControl.h"

#include "../Obstacles/ObstacleVector.h"
#include "../Obstacles/StefanFish.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

CubismUP_3D_NAMESPACE_BEGIN

SchoolControl::SchoolControl(SimulationData & s) : Operator(s)
{
  WangSchoolPIDParams p;
  p.k = s.schoolK;
  p.gammaAtt = s.schoolGammaAtt;
  p.gammaAli = s.schoolGammaAli;
  p.gammaR = s.schoolGammaR;
  p.decideS = s.schoolDecide;
  p.v0Cm = s.schoolV0Cm;
  p.kpPhi = s.schoolKpPhi;
  p.kdPhi = s.schoolKdPhi;
  p.kpV = s.schoolKpV;
  p.bMax = s.schoolBMax;
  p.aMax = s.schoolAMax;
  p.seed = s.schoolSeed;
  if (s.schoolLengthScale > 0) {
    p.lengthScale = s.schoolLengthScale;
    lengthScaleFromFish = false;
  }
  controller.setParams(p);
  logEvery = s.schoolLogEvery;

  if (sim.rank == 0 && logEvery > 0) {
    log.open("school_control.log");
    if (log) {
      log << "# step time fish_id x y phi phi_des speed v_des omega b b_held a a_held\n";
      log << "# mode=kinematics (Turn(b)+period(a) every 0.5*Tperiod)\n";
      log.flush();
    }
  }
}

void SchoolControl::operator()(const Real dt)
{
  if (!sim.schoolControl) return;

  if (sim.step == 0 && dt <= 0) return;

  auto & obs = sim.obstacle_vector->getObstacleVector();
  std::vector<StefanFish *> fish;
  fish.reserve(obs.size());
  for (const auto & o : obs) {
    auto * sf = dynamic_cast<StefanFish *>(o.get());
    if (sf) fish.push_back(sf);
  }
  const int n = (int)fish.size();
  if (n == 0) return;

  if ((int)tLastAct.size() != n) {
    tLastAct.assign(n, -1e9);
    bHeld.assign(n, 0);
    aHeld.assign(n, 0);
  }

  if (!warnedCorrectPosition) {
    for (auto * sf : fish) {
      if (sf->bCorrectPosition || sf->bCorrectPositionZ) {
        if (sim.rank == 0) {
          printf("[SchoolControl] WARNING: CorrectPosition is on; it fights act() turn bias. "
                 "Set CorrectPosition=0 CorrectPositionZ=0.\n");
          fflush(stdout);
        }
        warnedCorrectPosition = true;
        break;
      }
    }
    if (!warnedCorrectPosition) warnedCorrectPosition = true;
  }

  if (lengthScaleFromFish) {
    Real Lmean = 0;
    for (auto * sf : fish) Lmean += sf->length;
    Lmean /= (Real)n;
    controller.params.lengthScale = Lmean / 3.1;
    lengthScaleFromFish = false;
    if (sim.rank == 0) {
      printf("[SchoolControl] kinematics mode  N=%d lengthScale=L/3.1=%.6g (Lmean=%.4g)\n",
             n, (double)controller.params.lengthScale, (double)Lmean);
      printf("[SchoolControl] yaw=Turn(b), surge=period(a); both every 0.5*Tperiod.\n");
      fflush(stdout);
    }
  }

  std::vector<Real> x(n), y(n), phi(n), speed(n), omega(n);
  for (int i = 0; i < n; ++i) {
    auto * sf = fish[i];
    x[i] = sf->position[0];
    y[i] = sf->position[1];
    const auto ypr = sf->getYawPitchRoll();
    phi[i] = ypr[0];
    const Real ux = sf->transVel[0] - sim.uinf[0];
    const Real uy = sf->transVel[1] - sim.uinf[1];
    speed[i] = std::hypot(ux, uy);
    omega[i] = sf->angVel[2];
  }

  controller.compute(sim.time, x, y, phi, speed, omega);

  for (int i = 0; i < n; ++i) {
    auto * sf = fish[i];
    sf->bForcedInSimFrame[0] = false;
    sf->bForcedInSimFrame[1] = false;
    sf->bBlockRotation[2] = false;
    sf->appliedForce  = {{0, 0, 0}};
    sf->appliedTorque = {{0, 0, 0}};

    // Turn(b) and period(a) are wave-shift events — must NOT be called every CFD step.
    // Do not use StefanFish::act(): if z-velocity is forced (planar), act() zeros
    // actions[1] as "no pitching", which would wipe the period action.
    const Real actPeriod = std::max((Real)0.5 * sf->getLearnTPeriod(), (Real)1e-3);
    if (sim.time - tLastAct[i] >= actPeriod - (Real)1e-12) {
      bHeld[i] = controller.b(i);
      aHeld[i] = controller.a(i);
      auto * cFish = dynamic_cast<CurvatureDefinedFishData *>(sf->myFish);
      if (cFish == nullptr) {
        printf("SchoolControl: expected CurvatureDefinedFishData\n");
        abort();
      }
      const Real tNext = sim.time + actPeriod;
      cFish->action_curvature(sim.time, tNext, bHeld[i]);
      cFish->action_period(sim.time, tNext, aHeld[i]);
      tLastAct[i] = sim.time;
    }
  }

  if (sim.rank == 0 && log && logEvery > 0 && (sim.step % logEvery == 0)) {
    for (int i = 0; i < n; ++i) {
      log << sim.step << ' ' << sim.time << ' ' << i << ' '
          << x[i] << ' ' << y[i] << ' ' << phi[i] << ' ' << controller.phiDes(i) << ' '
          << speed[i] << ' ' << controller.vDes(i) << ' ' << omega[i] << ' '
          << controller.b(i) << ' ' << bHeld[i] << ' '
          << controller.a(i) << ' ' << aHeld[i] << '\n';
    }
    log.flush();
  }

  (void)dt;
}

CubismUP_3D_NAMESPACE_END
