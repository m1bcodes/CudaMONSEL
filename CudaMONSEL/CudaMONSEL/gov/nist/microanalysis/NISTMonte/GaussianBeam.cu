#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\EPQLibrary\PhysicalConstants.cuh"
#include "gov\nist\microanalysis\NISTMonte\GaussianBeam.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"

#include "Amphibian\random.cuh"

namespace GaussianBeam
{
   GaussianBeam::GaussianBeam(double width) :
      mWidth(width)
   {
      double mult[3];
      Math2::multiply3d(0.99 * MonteCarloSS::ChamberRadius, Math2::MINUS_Z_AXIS, mult);
      setCenter(mult);
   }

   __host__ __device__ GaussianBeam::GaussianBeam(double width, double energy, const double center[]) :
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

   __host__ __device__ double GaussianBeam::getBeamEnergy() const
   {
      return mBeamEnergy;
   }

   __host__ __device__ void GaussianBeam::setCenter(const double center[])
   {
      memcpy(mCenter, center, sizeof(double) * 3);
   }

   const double* GaussianBeam::getCenter() const
   {
      return mCenter;
   }

   __host__ __device__ ElectronT* GaussianBeam::createElectron() const
   {
      double initialPos[] = {
         mCenter[0],
         mCenter[1],
         mCenter[2]
      };
      double r = ::sqrt(-2. * ::log(Random::random())) * mWidth;
      double th = 2.0 * PhysicalConstants::PI * (Random::random());
      initialPos[0] += r * ::cos(th);
      initialPos[1] += r * ::sin(th);

      return new ElectronT(initialPos, mBeamEnergy); // TODO: handle this
   }
}