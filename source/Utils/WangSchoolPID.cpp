//
//  CubismUP_3D
//  Port of StefanPIDSchool social + PD (no plant) from simple_point_fish.py
//

#include "WangSchoolPID.h"

#include <algorithm>
#include <cmath>
#include <numeric>

CubismUP_3D_NAMESPACE_BEGIN

namespace {
inline Real clip(Real x, Real lo, Real hi)
{
  return std::max(lo, std::min(hi, x));
}
} // namespace

Real WangSchoolPID::wrap(Real a)
{
  const Real twopi = 2.0 * M_PI;
  a = std::fmod(a + M_PI, twopi);
  if (a < 0) a += twopi;
  return a - M_PI;
}

void WangSchoolPID::setParams(const WangSchoolPIDParams & p)
{
  params = p;
  rng.seed(p.seed);
  timersInit = false;
}

void WangSchoolPID::resize(int nNew)
{
  if (nNew == n && (int)B.size() == nNew) return;
  n = nNew;
  phi_des.assign(n, 0);
  v_des.assign(n, params.v0Cm * params.lengthScale);
  t_decide.assign(n, 0);
  tau_decide.assign(n, params.decideS);
  B.assign(n, 0);
  A.assign(n, 0);
  timersInit = false;
}

Real WangSchoolPID::social(int i,
                           const std::vector<Real> & x,
                           const std::vector<Real> & y,
                           const std::vector<Real> & phi) const
{
  const Real s = params.lengthScale;
  const Real lAtt = params.lAttCm * s;
  const Real lAli = params.lAliCm * s;
  const Real dAtt = params.dAttCm * s;
  const Real dAli = params.dAliCm * s;

  std::vector<Real> dphi(n, 0.0);
  for (int j = 0; j < n; ++j) {
    if (j == i) continue;
    const Real dx = x[j] - x[i];
    const Real dy = y[j] - y[i];
    const Real d = std::max(std::hypot(dx, dy), (Real)1e-6);
    const Real psi = wrap(std::atan2(dy, dx) - phi[i]);
    const Real dph = wrap(phi[j] - phi[i]);

    const Real f_att = params.gammaAtt * (d / dAtt - 1.0) / (1.0 + (d / lAtt) * (d / lAtt));
    const Real f_ali = params.gammaAli * (d / dAli + 1.0) * std::exp(-(d / lAli) * (d / lAli));
    const Real o_att = 1.395 * std::sin(psi) * (1.0 - 0.33 * std::cos(psi));
    const Real e_att = 0.9326 * (1.0 - 0.48 * std::cos(dph) - 0.31 * std::cos(2.0 * dph));
    const Real e_ali = 0.9012 * (1.0 + 0.6 * std::cos(psi) - 0.32 * std::cos(2.0 * psi));
    const Real o_ali = 1.6385 * std::sin(dph) * (1.0 + 0.3 * std::cos(2.0 * dph));
    dphi[j] = f_att * o_att * e_att + f_ali * e_ali * o_ali;
  }

  Real socialSum = 0;
  if (params.k >= n - 1) {
    socialSum = std::accumulate(dphi.begin(), dphi.end(), (Real)0);
  } else {
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    const int k = std::min(params.k, n - 1);
    std::partial_sort(idx.begin(), idx.begin() + k, idx.end(),
                      [&](int a, int b) { return std::fabs(dphi[a]) > std::fabs(dphi[b]); });
    for (int t = 0; t < k; ++t) socialSum += dphi[idx[t]];
  }

  std::normal_distribution<Real> gauss(0.0, 1.0);
  return socialSum + params.gammaR * gauss(rng);
}

void WangSchoolPID::pdAct(int i,
                          const std::vector<Real> & phi,
                          const std::vector<Real> & speed,
                          const std::vector<Real> & omega)
{
  const Real ePhi = wrap(phi_des[i] - phi[i]);
  const Real b = params.kpPhi * ePhi - params.kdPhi * omega[i];
  B[i] = clip(b, -params.bMax, params.bMax);
  const Real eV = v_des[i] - speed[i];
  // a<0 → shorter period → higher speed
  A[i] = clip(-params.kpV * eV, -params.aMax, params.aMax);
}

void WangSchoolPID::compute(Real time,
                            const std::vector<Real> & x,
                            const std::vector<Real> & y,
                            const std::vector<Real> & phi,
                            const std::vector<Real> & speed,
                            const std::vector<Real> & omega)
{
  const int nIn = (int)x.size();
  if (nIn <= 0) return;
  resize(nIn);

  const Real v0 = params.v0Cm * params.lengthScale;
  std::uniform_real_distribution<Real> uni(0.0, 1.0);
  std::normal_distribution<Real> gauss(0.0, 1.0);

  if (!timersInit) {
    for (int i = 0; i < n; ++i) {
      phi_des[i] = phi[i];
      v_des[i] = v0;
      tau_decide[i] = params.decideS;
      t_decide[i] = time - uni(rng) * params.decideS;
    }
    timersInit = true;
  }

  for (int i = 0; i < n; ++i) {
    if (time - t_decide[i] >= tau_decide[i]) {
      phi_des[i] = wrap(phi[i] + social(i, x, y, phi));
      // Wang kick band is 6–22 cm/s around v0=14 cm/s. Scale that band with
      // the chosen v0 so lowering -schoolV0Cm actually lowers v_des (otherwise
      // a hard floor at 6*lengthScale keeps demanding unreachable speed).
      const Real vRaw = v0 + 0.25 * v0 * gauss(rng);
      const Real vLo = (6.0 / 14.0) * v0;
      const Real vHi = (22.0 / 14.0) * v0;
      v_des[i] = clip(vRaw, vLo, vHi);
      tau_decide[i] = clip(params.decideS + 0.12 * gauss(rng), 0.05, 1.2);
      t_decide[i] = time;
    }
    pdAct(i, phi, speed, omega);
  }
}

CubismUP_3D_NAMESPACE_END
