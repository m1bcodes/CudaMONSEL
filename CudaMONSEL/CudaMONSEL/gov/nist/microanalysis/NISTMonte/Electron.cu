#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

namespace Electron
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __device__ static long lastID = 0; // ID of last generated electron
#else
   static long lastID = 0; // ID of last generated electron
#endif

   long getlastIdent()
   {
      return lastID;
   }

   __host__ __device__ void Electron::init(const double initialPos[], double theta, double phi, double kE)
   {
      memcpy(mPosition, initialPos, sizeof(double) * 3);
      memcpy(mPrevPosition, initialPos, sizeof(double) * 3);
      mScatteringElement = (nullptr);
      mCurrentRegion = (nullptr);
      mPrevRegion = (nullptr);
      mEnergy = (kE);
      previousEnergy = (kE);
      mTheta = (theta);
      mPhi = (phi);
      mStepCount = (0);
      mTrajectoryComplete = (false);
      ident = (++lastID);
      parentID = 0;
   }

   __host__ __device__ Electron::Electron(const double initialPos[], double kE)
   {
      init(initialPos, 0., 0., kE);
   }

   Electron::Electron(const double initialPos[], double theta, double phi, double kE)
   {
      init(initialPos, theta, phi, kE);
   }

   __host__ __device__ Electron::Electron(const Electron& parent, double theta, double phi, double kE)
   {
      const double* pos = parent.getPosition();
      init(pos, theta, phi, kE);
      mCurrentRegion = parent.getCurrentRegion();
      mPrevRegion = mCurrentRegion;
      parentID = parent.getIdent();
   }

   __host__ __device__ void Electron::setDirection(double theta, double phi)
   {
      mTheta = theta;
      mPhi = phi;
   }

   __host__ __device__ const double * Electron::getPosition() const
   {
      return mPosition;
   }

   __host__ __device__ void Electron::setPosition(const double newpos[])
   {
      memcpy(mPosition, newpos, sizeof(double) * 3);
   }

   __host__ __device__ const double * Electron::getPrevPosition() const
   {
      return mPrevPosition;
   }

   __host__ __device__ const RegionBaseT* Electron::getCurrentRegion() const
   {
      return mCurrentRegion;
   }

   const RegionBaseT* Electron::getPreviousRegion() const
   {
      return mPrevRegion;
   }

   __host__ __device__ double Electron::getEnergy() const
   {
      return mEnergy;
   }

   __host__ __device__ double Electron::getPreviousEnergy() const
   {
      return previousEnergy;
   }

   __host__ __device__ int Electron::getStepCount() const
   {
      return mStepCount;
   }

   double Electron::stepLength() const
   {
      return Math2::distance3d(mPrevPosition, mPosition);
   }

   __host__ __device__ void Electron::candidatePoint(const double dS, double res[]) const
   {
      const double st = ::sin(mTheta);
      // Calculate the new point as dS distance from mPosition
      res[0] = mPosition[0] + dS * ::cos(mPhi) * st;
      res[1] = mPosition[1] + dS * ::sin(mPhi) * st;
      res[2] = mPosition[2] + dS * ::cos(mTheta);
   }

   __host__ __device__ void Electron::updateDirection(const double dTheta, const double dPhi)
   {
      // The candidate point is computed by rotating the current trajectory back
      // to the z-axis, deflecting the z-axis by dTheta down from the z-axis and
      // dPhi around the z-axis, then finally rotating back to the original
      // trajectory.

      double ct = ::cos(mTheta), st = ::sin(mTheta);
      double cp = ::cos(mPhi), sp = ::sin(mPhi);
      double ca = ::cos(dTheta), sa = ::sin(dTheta);
      double cb = ::cos(dPhi);

      double xx = cb * ct * sa + ca * st;
      double yy = sa * ::sin(dPhi);
      double dx = cp * xx - sp * yy;
      double dy = cp * yy + sp * xx;
      double dz = ca * ct - cb * sa * st;

      mTheta = ::atan2(::sqrt(dx * dx + dy * dy), dz);
      mPhi = ::atan2(dy, dx);
   }

   __host__ __device__ void Electron::move(const double newPoint[], double dE)
   {
      // Update mPrevPosition and then mPosition
      memcpy(mPrevPosition, mPosition, sizeof(double) * 3);
      memcpy(mPosition, newPoint, sizeof(double) * 3);

      // Update the energy
      previousEnergy = mEnergy;
      mEnergy += dE;
      ++mStepCount;
   }

   __host__ __device__ void Electron::setEnergy(double newEnergy)
   {
      mEnergy = newEnergy;
   }

   __host__ __device__ void Electron::setPreviousEnergy(double newPreviousEnergy)
   {
      previousEnergy = newPreviousEnergy;
   }

   __host__ __device__ void Electron::setCurrentRegion(const RegionBaseT* reg)
   {
      mPrevRegion = mCurrentRegion;
      mCurrentRegion = reg;
   }

   __host__ __device__ const ElementT* Electron::getScatteringElement() const
   {
      return mScatteringElement;
   }

   __host__ __device__ void Electron::setScatteringElement(const ElementT* scatteringElement)
   {
      mScatteringElement = scatteringElement;
   }

   __host__ __device__ double Electron::getPhi() const
   {
      return mPhi;
   }

   __host__ __device__ double Electron::getTheta() const
   {
      return mTheta;
   }

   __host__ __device__ bool Electron::isTrajectoryComplete() const
   {
      return mTrajectoryComplete;
   }

   __host__ __device__ void Electron::setTrajectoryComplete(bool trajectoryComplete)
   {
      mTrajectoryComplete = trajectoryComplete;
   }

   __host__ __device__ long Electron::getIdent() const
   {
      return ident;
   }

   long Electron::getParentID() const
   {
      return parentID;
   }
}