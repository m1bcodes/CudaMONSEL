#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\EPQLibrary\PhysicalConstants.cuh"
#include "gov\nist\microanalysis\NISTMonte\GaussianBeam.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"

namespace GaussianBeam
{
   GaussianBeam::GaussianBeam(double width)
   {
      mCenter = Math2::multiply(0.99 * MonteCarloSS::ChamberRadius, PositionVecT(Math2::MINUS_Z_AXIS, Math2::MINUS_Z_AXIS + 3));
      mWidth = width;
   }

   void GaussianBeam::setWidth(double width)
   {
      mWidth = width;
   }

   double GaussianBeam::getWidth()
   {
      return mWidth;
   }

   void GaussianBeam::setBeamEnergy(double beamEnergy)
   {
      mBeamEnergy = beamEnergy;
   }

   double GaussianBeam::getBeamEnergy() const
   {
      return mBeamEnergy;
   }

   void GaussianBeam::setCenter(const double center[])
   {
      mCenter.assign(center, center + 3);
   }

   PositionVecT GaussianBeam::getCenter() const
   {
      return mCenter;
   }

   ElectronT* GaussianBeam::createElectron() const
   {
      double initialPos[] = {
         mCenter[0],
         mCenter[1],
         mCenter[2]
      };
      double r = ::sqrt(-2. * ::log(Math2::random())) * mWidth;
      double th = 2.0 * PhysicalConstants::PI * (Math2::random());
      initialPos[0] += r * ::cos(th);
      initialPos[1] += r * ::sin(th);

      return new ElectronT(initialPos, mBeamEnergy);
   }
}