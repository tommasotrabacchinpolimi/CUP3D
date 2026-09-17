//
//  CubismUP_3D
//  Wang social setpoints + planar PD → act({b, a}) for StefanFish.
//  Port of CUP3D/prototypes/simple_point_fish.py :: StefanPIDSchool
//  (control only; CFD is the plant).
//
//  Length scaling: Wang paper uses cm. Convert with lengthScale =
//  L_fish / 3.1 (default) so distances and speeds map into sim units.
//  Override with -schoolLengthScale.
//

#pragma once

#include "../Definitions.h"

#include <random>
#include <vector>

CubismUP_3D_NAMESPACE_BEGIN

struct WangSchoolPIDParams
{
  int k = 1;
  Real gammaAtt = 0.04;
  Real gammaAli = 0.20;
  Real gammaR = 0.20;
  // Wang cm defaults; multiplied by lengthScale → sim units
  Real lAttCm = 28.0;
  Real lAliCm = 28.0;
  Real dAttCm = 3.0;
  Real dAliCm = 6.0;
  Real v0Cm = 14.0; // cruise in Wang-cm; kick band scales as (6/14–22/14)*v0
  Real lengthScale = 0.2 / 3.1; // L=0.2 fish default
  Real decideS = 0.5;
  Real kpPhi = 1.2;
  Real kdPhi = 0.25;
  Real kpV = 0.08;
  Real bMax = 1.5;
  Real aMax = 0.5;
  unsigned seed = 71;
};

/** Outer school controller: refresh (phi_des, v_des) on decide clock, PD → b,a. */
class WangSchoolPID
{
public:
  WangSchoolPID() = default;
  explicit WangSchoolPID(const WangSchoolPIDParams & p) : params(p), rng(p.seed) {}

  void setParams(const WangSchoolPIDParams & p);
  void resize(int n);

  /** Update actions from current planar poses. Does not integrate kinematics. */
  void compute(Real time,
               const std::vector<Real> & x,
               const std::vector<Real> & y,
               const std::vector<Real> & phi,
               const std::vector<Real> & speed,
               const std::vector<Real> & omega);

  int size() const { return n; }
  Real b(int i) const { return B[i]; }
  Real a(int i) const { return A[i]; }
  Real phiDes(int i) const { return phi_des[i]; }
  Real vDes(int i) const { return v_des[i]; }

  WangSchoolPIDParams params;

private:
  static Real wrap(Real a);
  Real social(int i,
              const std::vector<Real> & x,
              const std::vector<Real> & y,
              const std::vector<Real> & phi) const;
  void pdAct(int i,
             const std::vector<Real> & phi,
             const std::vector<Real> & speed,
             const std::vector<Real> & omega);

  int n = 0;
  mutable std::mt19937 rng;
  std::vector<Real> phi_des, v_des, t_decide, tau_decide, B, A;
  bool timersInit = false;
};

CubismUP_3D_NAMESPACE_END
