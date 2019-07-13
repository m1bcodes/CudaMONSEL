#ifndef _ELECTRON_GUN_CUH_
#define _ELECTRON_GUN_CUH_

#include "Declarations.cuh"

namespace ElectronGun
{
   class ElectronGun
   {
   public:
      virtual void setBeamEnergy(double beamEnergy) = 0;
      __host__ __device__ virtual double getBeamEnergy() const = 0;
      __host__ __device__ virtual void setCenter(const double center[]) = 0;
      virtual const double* getCenter() const = 0;
      __host__ __device__ virtual ElectronT* createElectron() const = 0;
   };
}

typedef ElectronGun::ElectronGun ElectronGunT;

#endif