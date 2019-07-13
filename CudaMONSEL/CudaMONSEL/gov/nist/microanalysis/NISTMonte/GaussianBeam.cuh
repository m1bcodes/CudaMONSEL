#ifndef _GAUSSIAN_BEAM_CUH_
#define _GAUSSIAN_BEAM_CUH_

#include "Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\ElectronGun.cuh"

namespace GaussianBeam
{
   class GaussianBeam : public ElectronGunT
   {
   public:
      GaussianBeam(double width);
      __host__ __device__ GaussianBeam(double width, double energy, const double center[]);

      void setWidth(double width);
      double getWidth();

      void setBeamEnergy(double beamEnergy) override;
      __host__ __device__ double getBeamEnergy() const override;
      __host__ __device__ void setCenter(const double center[]) override;
      const double* getCenter() const override;
      __host__ __device__ ElectronT* createElectron() const override;

   private:
      // transient private Random mRandom = new Random();

      double mCenter[3];
      double mBeamEnergy;
      double mWidth = 1.0e-9;
   };
}

#endif