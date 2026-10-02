#ifndef CubismUP_3D_BurstandCoastFish_h
#define CubismUP_3D_BurstandCoastFish_h

#include "Fish.h"

CubismUP_3D_NAMESPACE_BEGIN

class BurstandCoastFishMidlineData;

class BurstandCoastFish: public Fish
{
 public:
  BurstandCoastFish(SimulationData&s, cubism::ArgumentParser&p);
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_BurstandCoastFish_h
