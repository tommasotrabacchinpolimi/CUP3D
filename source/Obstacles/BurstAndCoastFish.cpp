//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "BurstandCoastFish.h"
#include "FishLibrary.h"
#include "FishShapes.h"
#include <vector>
#include <Cubism/ArgumentParser.h>
#include <cmath>


CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

class BurstandCoastFishMidlineData : public FishMidlineData
{

 public:
  // L=length, T=period, phi=phase shift, _h=grid size, A=amplitude modulation
  BurstandCoastFishMidlineData(Real L, Real T, Real phi, Real _h, Real A) :
  FishMidlineData(L,T,phi,_h,A),
  c2(0.1625), c1(-0.0825), c0(0.02),
  lmbda(0.2 * L), lmbda_bout(2.0 * L), T_bout(T), c((2.0 * L) / T)
  {
  }

  virtual void computeMidline(const Real t, const Real dt) override;

 private:
  Real c2, c1, c0;
  Real lmbda, lmbda_bout, T_bout;
  Real c;

  Real a(Real x)
  {
    const Real n = x / length;
    return amplitudeFactor * (c2 * n * n + c1 * n + c0) * length;
  }
  
  Real da_dx(Real x)
  {
    return amplitudeFactor * (2.0 * c2 * x / length + c1);
  }
  
  Real Y(Real x)
  {
    x = std::fmod(x, lmbda_bout);
    if (x < 0) x += lmbda_bout;
    if (x >= 0 && x <= lmbda_bout - lmbda)
      return 0.0;
    else if (x >= lmbda_bout - lmbda && x <= lmbda_bout - 0.75 * lmbda)
      return -0.5 * (1.0 - std::cos((4.0 * M_PI * (x - lmbda_bout)) / lmbda));
    else if (x >= lmbda_bout - 0.75 * lmbda && x <= lmbda_bout - 0.25 * lmbda)
      return -std::sin((2.0 * M_PI * (x - lmbda_bout)) / lmbda);
    else if (x >= lmbda_bout - 0.25 * lmbda && x <= lmbda_bout)
      return 0.5 * (1.0 - std::cos((4.0 * M_PI * (x - lmbda_bout)) / lmbda));
    else
      return 0.0;
  }
  
  Real dY_dx(Real x)
  {
    x = std::fmod(x, lmbda_bout);
    if (x < 0) x += lmbda_bout;
    const Real theta = (x - lmbda_bout) / lmbda;
    const Real k = 2.0 * M_PI / lmbda;
    if (x >= 0 && x <= lmbda_bout - lmbda)
      return 0.0;
    else if (x >= lmbda_bout - lmbda && x <= lmbda_bout - 0.75 * lmbda)
      return -k * std::sin(4.0 * M_PI * theta);
    else if (x >= lmbda_bout - 0.75 * lmbda && x <= lmbda_bout - 0.25 * lmbda)
      return -k * std::cos(2.0 * M_PI * theta);
    else if (x >= lmbda_bout - 0.25 * lmbda && x <= lmbda_bout)
      return k * std::sin(4.0 * M_PI * theta);
    else
      return 0.0;
  }

  Real y(Real x, Real t)
  {
    return a(x) * Y(x - c * t);
  }

  Real dy_dx(Real x, Real t)
  {
    return da_dx(x) * Y(x - c * t) + a(x) * dY_dx(x - c * t);
  }

  Real dy_dt(Real x, Real t)
  {
    return a(x) * dY_dx(x - c * t) * (-c);
  }


  // returns an array x, such that x[i] is the x coordinate corresponding to the endl/(dl*i) arch length
  void computeX(std::vector<Real>& x, Real dl, Real endl, Real t)
  {
    x[0] = 0.0;
    for (int i = 1; i < x.size(); i++) {
      x[i] = x[i-1] + dl / std::sqrt(1 + dy_dx(x[i-1], t)*dy_dx(x[i-1], t));
    }
  }

};

void BurstandCoastFishMidlineData::computeMidline(const Real t,const Real dt)
{
  const Real ds = length/Nm;
  const Real dl = ds/100;
  std::vector<Real> x((size_t)std::floor(length / dl) + 2);
  computeX(x, dl, length, t);

  rX[0] = 0.0;
  rY[0] = y(0.0, t);
  vX[0] = 0.0;
  vY[0] = dy_dt(0.0, t);
  rZ[0] = 0.0;
  vZ[0] = 0.0;
  for (int i = 1; i < Nm; i++) {
    int index = (int)(rS[i] / dl);
    if (index < 0) index = 0;
    if (index >= (int)x.size()) index = (int)x.size() - 1;
    rX[i] = x[index];
    rY[i] = y(rX[i], t);
    const Real ft = dy_dt(rX[i], t);
    const Real fx = dy_dx(rX[i], t);
    const Real delta_y = rY[i] - rY[i-1];
    const Real delta_x = rX[i] - rX[i-1];
    vX[i] = (1.0 / (1.0 + delta_y / delta_x * fx))
          * (vX[i-1] - delta_y / delta_x * (ft - vY[i-1]));
    vY[i] = ft + fx * vX[i];
    rZ[i] = 0.0;
    vZ[i] = 0.0;
  }

  #pragma omp parallel for schedule(static)
  for(int i=0; i<Nm-1; i++) {
    const Real ds = rS[i+1]-rS[i];
    const Real tX = rX[i+1]-rX[i];
    const Real tY = rY[i+1]-rY[i];
    const Real tVX = vX[i+1]-vX[i];
    const Real tVY = vY[i+1]-vY[i];
    norX[i] = -tY/ds;
    norY[i] =  tX/ds;
    norZ[i] =  0.0;
    vNorX[i] = -tVY/ds;
    vNorY[i] =  tVX/ds;
    vNorZ[i] = 0.0;
    binX[i] =  0.0;
    binY[i] =  0.0;
    binZ[i] =  1.0;
    vBinX[i] = 0.0;
    vBinY[i] = 0.0;
    vBinZ[i] = 0.0;
  }
  norX[Nm-1] = norX[Nm-2];
  norY[Nm-1] = norY[Nm-2];
  norZ[Nm-1] = norZ[Nm-2];
  vNorX[Nm-1] = vNorX[Nm-2];
  vNorY[Nm-1] = vNorY[Nm-2];
  vNorZ[Nm-1] = vNorZ[Nm-2];
  binX[Nm-1] = binX[Nm-2];
  binY[Nm-1] = binY[Nm-2];
  binZ[Nm-1] = binZ[Nm-2];
  vBinX[Nm-1] = vBinX[Nm-2];
  vBinY[Nm-1] = vBinY[Nm-2];
  vBinZ[Nm-1] = vBinZ[Nm-2];
}

BurstandCoastFish::BurstandCoastFish(SimulationData&s, ArgumentParser&p) : Fish(s, p)
{
  const Real ampFac = p("-amplitudeFactor").asDouble(1.0);
  const Real Tperiod = p("-T").asDouble(1.0);
  const Real phaseShift = p("-phi").asDouble(0.0);

  BurstandCoastFishMidlineData* localFish = new BurstandCoastFishMidlineData(
    length, Tperiod, phaseShift, sim.hmin, ampFac);

  assert(myFish == nullptr);
  myFish = (FishMidlineData*) localFish;

  std::string heightName = p("-heightProfile").asString("baseline");
  std::string  widthName = p( "-widthProfile").asString("baseline");
  MidlineShapes::computeWidthsHeights(heightName, widthName, length,
    myFish->rS, myFish->height, myFish->width, myFish->Nm, sim.rank);

  if(!sim.rank)
    printf("BurstandCoastFish: N:%d, L:%f, T:%f, phi:%f, amplitude:%f\n",
        myFish->Nm, length, Tperiod, phaseShift, ampFac);
}

CubismUP_3D_NAMESPACE_END
