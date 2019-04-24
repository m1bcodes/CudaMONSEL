#ifndef _GAUSSIAN_BEAM_CUH_
#define _GAUSSIAN_BEAM_CUH_

#include "Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\ElectronGun.cuh"

namespace GaussianBeam
{
   class GaussianBeam : public ElectronGunT
   {
   public:
      GaussianBeam::GaussianBeam(double width);
      void setWidth(double width);
      double getWidth();

      void setBeamEnergy(double beamEnergy) override;
      double getBeamEnergy() const override;
      void setCenter(const double center[]) override;
      PositionVecT getCenter() const override;
      ElectronT createElectron() const override;

   private:
      // transient private Random mRandom = new Random();

      PositionVecT mCenter;
      double mBeamEnergy;
      double mWidth = 1.0e-9;
   };
}

#endif