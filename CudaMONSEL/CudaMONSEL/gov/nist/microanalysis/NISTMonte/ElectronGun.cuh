#ifndef _ELECTRON_GUN_CUH_
#define _ELECTRON_GUN_CUH_

#include "Declarations.cuh"

namespace ElectronGun
{
   class ElectronGun
   {
   public:
      virtual void setBeamEnergy(double beamEnergy) = 0;
      virtual double getBeamEnergy() const = 0;
      virtual void setCenter(const double center[]) = 0;
      virtual PositionVecT getCenter() const = 0;
      virtual ElectronT createElectron() const = 0;
   };
}

typedef ElectronGun::ElectronGun ElectronGunT;

#endif