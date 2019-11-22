#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\EPQLibrary\PhysicalConstants.cuh"
#include "gov\nist\microanalysis\NISTMonte\GaussianBeam.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"

#include "Amphibian\random.cuh"

namespace GaussianBeam
{
   GaussianBeam::GaussianBeam(double width) :
      mWidth(width),
      mPhi(0.),
      mTheta(0.),
      mFocalLength(20.)
   {
      double mult[3];
      Math2::multiply3d(0.99 * MonteCarloSS::ChamberRadius, Math2::MINUS_Z_AXIS, mult);
      setCenter(mult);
   }

   __host__ __device__ GaussianBeam::GaussianBeam(const double width, const double energy, const double center[]) :
      mWidth(width),
      mBeamEnergy(energy),
      mPhi(0.),
      mTheta(0.),
      mFocalLength(20.)
   {
      setCenter(center);
   }

   __host__ __device__ GaussianBeam::GaussianBeam(const double width, const double energy, const double theta, const double phi, const double center[], const double focalLength) :
      mWidth(width),
      mBeamEnergy(energy),
      mPhi(phi),
      mTheta(theta),
      mFocalLength(focalLength)
   {
      if (mFocalLength <= 0.) printf("GaussianBeam::GaussianBeam: bad mFocalLength (%.5e)\n", mFocalLength);
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
      memcpy(mCenter, center, sizeof(mCenter[0]) * 3);
   }

   const double* GaussianBeam::getCenter() const
   {
      return mCenter;
   }

   //__host__ __device__ ElectronT* GaussianBeam::createElectron() const
   //{
   //   double initialPos[] = {
   //      mCenter[0],
   //      mCenter[1],
   //      mCenter[2]
   //   };
   //   double rand = 0.;
   //   while (rand == 0.) {
   //      rand = Random::random();
   //   }
   //   const double r = ::sqrt(-2. * ::log(rand)) * mWidth;
   //   const double th = 2.0 * PhysicalConstants::PI * (Random::random());
   //   initialPos[0] += r * ::cos(th);
   //   initialPos[1] += r * ::sin(th);

   //   ElectronT* newElectron = new ElectronT(initialPos, mTheta, mPhi, mBeamEnergy);
   //   if (!newElectron) printf("ElectronT* GaussianBeam::createElectron: failed creating electron.\n");
   //   return newElectron;
   //}

   __host__ __device__ ElectronT* GaussianBeam::createElectron() const
   {
      //const double focalLength = 20. / ToSI::GIGA;

      double p[] = {
         mCenter[0],
         mCenter[1],
         mCenter[2]
      };

      double rand = 0.;
      while (rand == 0.) {
         rand = Random::random();
      }

      const double r = ::sqrt(-2. * ::log(rand)) * mWidth;
      if (r != r) {
         printf("GaussianBeam::createElectron: r\n");
      }
      const double th = 2.0 * PhysicalConstants::PI * Random::random();
      p[0] += r * ::cos(th);
      p[1] += r * ::sin(th);

      double pf[] = {
         mCenter[0] + mFocalLength * ::sin(mTheta) * ::cos(mPhi),
         mCenter[1] + mFocalLength * ::sin(mTheta) * ::sin(mPhi),
         mCenter[2] + mFocalLength * ::cos(mTheta)
      };

      rand = 0.;
      while (rand == 0.) {
         rand = Random::random();
      }
      //const double rf = ::sqrt(-2. * ::log(rand)) * mWidth/2.;
      //const double thf = 2.0 * PhysicalConstants::PI * Random::random();
      const double rf = r/2.;
      const double thf = th;
      pf[0] += rf * ::cos(thf);
      pf[1] += rf * ::sin(thf);
      //pf[0] += r * ::cos(th);
      //pf[1] += r * ::sin(th);

      double u[3];
      Math2::minus3d(pf, p, u);
      Math2::normalize3d(u, u);

      const double phi = ::atan2(u[1], u[0]);
      const double l = ::sqrt(u[1] * u[1] + u[0] * u[0]);
      const double theta = ::atan2(l, u[2]);
      
      if (phi != phi) {
         printf("GaussianBeam::createElectron: phi\n");
      }
      if (theta != theta) {
         printf("GaussianBeam::createElectron: theta\n");
      }

      ElectronT* newElectron = new ElectronT(p, theta, phi, mBeamEnergy);
      if (!newElectron) printf("ElectronT* GaussianBeam::createElectron: failed creating electron.\n");
      return newElectron;
   }
}