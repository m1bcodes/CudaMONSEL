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

   __host__ __device__ void Electron::init(const double initialPos[], float theta, float phi, float kE)
   {
      memcpy(mPosition, initialPos, sizeof(mPosition[0]) * 3);
      memcpy(mPrevPosition, initialPos, sizeof(mPrevPosition[0]) * 3);
      //mPosition[0] = initialPos[0]; mPosition[1] = initialPos[1]; mPosition[2] = initialPos[2];
      //mPrevPosition[0] = initialPos[0]; mPrevPosition[1] = initialPos[1]; mPrevPosition[2] = initialPos[2];

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

   __host__ __device__ Electron::Electron(const double initialPos[], float kE) : initialEnergy(kE)
   {
      init(initialPos, 0., 0., kE);
   }

   __host__ __device__ Electron::Electron(const double initialPos[], float theta, float phi, float kE) : initialEnergy(kE)
   {
      init(initialPos, theta, phi, kE);
   }

   __host__ __device__ Electron::Electron(const Electron& parent, float theta, float phi, float kE) : initialEnergy(kE)
   {
      init(parent.getPosition(), theta, phi, kE);
      mCurrentRegion = parent.getCurrentRegion();
      mPrevRegion = mCurrentRegion;
      parentID = parent.getIdent();
   }

   __host__ __device__ void Electron::setDirection(float theta, float phi)
   {
      mTheta = theta;
      mPhi = phi;

      if (mTheta != mTheta)
         printf("RegionBase::setDirection: mTheta");
      if (mPhi != mPhi)
         printf("RegionBase::setDirection: mPhi");
   }

   __host__ __device__ const double * Electron::getPosition() const
   {
      return mPosition;
   }

   __host__ __device__ void Electron::setPosition(const double newpos[])
   {
      mPosition[0] = newpos[0];
      mPosition[1] = newpos[1];
      mPosition[2] = newpos[2];
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

   __host__ __device__ float Electron::getEnergy() const
   {
      return mEnergy;
   }

   __host__ __device__ float Electron::getPreviousEnergy() const
   {
      return previousEnergy;
   }

   __host__ __device__ int Electron::getStepCount() const
   {
      return mStepCount;
   }

   float Electron::stepLength() const
   {
      return Math2::distance3d(mPrevPosition, mPosition);
   }

   __host__ __device__ void Electron::candidatePoint(const float dS, double res[]) const
   {
      const float st = ::sinf(mTheta);
      // Calculate the new point as dS distance from mPosition
      res[0] = mPosition[0] + dS * ::cosf(mPhi) * st;
      res[1] = mPosition[1] + dS * ::sinf(mPhi) * st;
      res[2] = mPosition[2] + dS * ::cosf(mTheta);

      if (res[0] != res[0] || res[1] != res[1] || res[2] != res[2])
         printf("RegionBase::candidatePoint: (%.5e, %.5e, %.5e)\n", res[0], res[1], res[2]);
   }

   __host__ __device__ void Electron::updateDirection(const float dTheta, const float dPhi)
   {
      // The candidate point is computed by rotating the current trajectory back
      // to the z-axis, deflecting the z-axis by dTheta down from the z-axis and
      // dPhi around the z-axis, then finally rotating back to the original
      // trajectory.

      const float ct = ::cosf(mTheta), st = ::sinf(mTheta);
      const float cp = ::cosf(mPhi), sp = ::sinf(mPhi);
      const float ca = ::cosf(dTheta), sa = ::sinf(dTheta);
      const float cb = ::cosf(dPhi);

      const float xx = cb * ct * sa + ca * st;
      const float yy = sa * ::sinf(dPhi);
      const float dx = cp * xx - sp * yy;
      const float dy = cp * yy + sp * xx;
      const float dz = ca * ct - cb * sa * st;

      mTheta = ::atan2f(::sqrtf(dx * dx + dy * dy), dz);
      mPhi = ::atan2f(dy, dx);

      if (mTheta != mTheta)
         printf("RegionBase::updateDirection: mTheta");
      if (mPhi != mPhi)
         printf("RegionBase::updateDirection: mPhi");
   }

   __host__ __device__ void Electron::move(const double newPoint[], float dE)
   {
      // Update mPrevPosition and then mPosition
      //mPrevPosition[0] = mPosition[0]; mPrevPosition[1] = mPosition[1]; mPrevPosition[2] = mPosition[2];
      //mPosition[0] = newPoint[0]; mPosition[1] = newPoint[1]; mPosition[2] = newPoint[2];
      memcpy(mPrevPosition, mPosition, sizeof(mPrevPosition[0]) * 3);
      memcpy(mPosition, newPoint, sizeof(mPosition[0]) * 3);

      // Update the energy
      previousEnergy = mEnergy;
      mEnergy += dE;
      ++mStepCount;
   }

   __host__ __device__ void Electron::setEnergy(float newEnergy)
   {
      mEnergy = newEnergy;
   }

   __host__ __device__ void Electron::setPreviousEnergy(float newPreviousEnergy)
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

   __host__ __device__ float Electron::getPhi() const
   {
      return mPhi;
   }

   __host__ __device__ float Electron::getTheta() const
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

   __host__ __device__ float Electron::getInitialEnergy() const
   {
      return initialEnergy;
   }

   long Electron::getParentID() const
   {
      return parentID;
   }
}