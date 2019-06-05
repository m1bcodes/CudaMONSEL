#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\EPQLibrary\PhysicalConstants.cuh"
#include "gov\nist\microanalysis\NISTMonte\GaussianBeam.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"

namespace GaussianBeam
{
   GaussianBeam::GaussianBeam(double width) :
      mWidth(width)
   {
      double mult[3];
      Math2::multiply3d(0.99 * MonteCarloSS::ChamberRadius, Math2::MINUS_Z_AXIS, mult);
      setCenter(mult);
   }

   GaussianBeam::GaussianBeam(double width, double energy, const double center[]) :
      mWidth(width),
      mBeamEnergy(energy)
   {
      setCenter(center);
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
      memcpy(mCenter, center, sizeof(double) * 3);
   }

   const double* GaussianBeam::getCenter() const
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

      return new ElectronT(initialPos, mBeamEnergy); // TODO: handle this
   }
}