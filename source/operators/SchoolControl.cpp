//
//  CubismUP_3D
//  Hybrid school control:
//    surge  → appliedForce along heading (track v_des) every step
//    yaw    → StefanFish::act({b}) on half-period clock (Turn is a wave event)
//  No appliedTorque; period action a unused (force owns speed).
//

#include "SchoolControl.h"

#include "../Obstacles/ObstacleVector.h"
#include "../Obstacles/StefanFish.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

CubismUP_3D_NAMESPACE_BEGIN

namespace {
inline Real clip(Real x, Real lo, Real hi)
{
  return std::max(lo, std::min(hi, x));
}
} // namespace

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
      log << "# step time fish_id x y phi phi_des speed v_des omega b b_held Fx Fy\n";
      log << "# mode=hybrid (F every step; act({b}) every 0.5*Tperiod)\n";
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
      printf("[SchoolControl] hybrid mode  N=%d lengthScale=L/3.1=%.6g (Lmean=%.4g)\n",
             n, (double)controller.params.lengthScale, (double)Lmean);
      printf("[SchoolControl] surge=F(v_des) every step; yaw=act({b}) every 0.5*Tperiod (Turn held).\n");
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

  const Real kpV = controller.params.kpV;
  const Real aMax = controller.params.aMax;

  for (int i = 0; i < n; ++i) {
    auto * sf = fish[i];
    sf->bForcedInSimFrame[0] = false;
    sf->bForcedInSimFrame[1] = false;
    sf->bBlockRotation[2] = false;

    const Real eV = controller.vDes(i) - speed[i];
    const Real aSurge = clip(kpV * eV, -aMax, aMax);
    const Real c = std::cos(phi[i]);
    const Real s = std::sin(phi[i]);
    sf->appliedForce[0] = sf->mass * aSurge * c;
    sf->appliedForce[1] = sf->mass * aSurge * s;
    sf->appliedForce[2] = 0;
    sf->appliedTorque = {{0, 0, 0}};

    // StefanFish API: take actions once every 0.5*getLearnTPeriod().
    // Turn(b) shifts the bend wave — must NOT be called every CFD step.
    const Real actPeriod = std::max((Real)0.5 * sf->getLearnTPeriod(), (Real)1e-3);
    if (sim.time - tLastAct[i] >= actPeriod - (Real)1e-12) {
      bHeld[i] = controller.b(i);
      const Real tNext = sim.time + actPeriod;
      fish[i]->act(tNext, std::vector<Real>{bHeld[i]});
      tLastAct[i] = sim.time;
    }
  }

  if (sim.rank == 0 && log && logEvery > 0 && (sim.step % logEvery == 0)) {
    for (int i = 0; i < n; ++i) {
      auto * sf = fish[i];
      log << sim.step << ' ' << sim.time << ' ' << i << ' '
          << x[i] << ' ' << y[i] << ' ' << phi[i] << ' ' << controller.phiDes(i) << ' '
          << speed[i] << ' ' << controller.vDes(i) << ' ' << omega[i] << ' '
          << controller.b(i) << ' ' << bHeld[i] << ' '
          << sf->appliedForce[0] << ' ' << sf->appliedForce[1] << '\n';
    }
    log.flush();
  }

  (void)dt;
}

CubismUP_3D_NAMESPACE_END
