#ifndef _ELECTRON_CUH_
#define _ELECTRON_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace Electron
{
   class Electron
   {
   public:
      __host__ __device__ Electron(const double initialPos[], float kE);
      __host__ __device__ Electron(const double initialPos[], float theta, float phi, float kE);
      __host__ __device__ Electron(const Electron& parent, float theta, float phi, float kE);

      __host__ __device__ void init(const double initialPos[], float theta, float phi, float kE);

      __host__ __device__ void setDirection(float theta, float phi);
      __host__ __device__ const double * getPosition() const;
      __host__ __device__ void setPosition(const double newpos[]);
      __host__ __device__ const double * getPrevPosition() const;
      __host__ __device__ const RegionBaseT* getCurrentRegion() const;
      const RegionBaseT* getPreviousRegion() const;
      __host__ __device__ float getEnergy() const;
      __host__ __device__ float getPreviousEnergy() const;
      __host__ __device__ int getStepCount() const;
      float stepLength() const;
      __host__ __device__ void candidatePoint(const float dS, double res[]) const;
      __host__ __device__ void updateDirection(const float dTheta, const float dPhi);
      __host__ __device__ void move(const double newPoint[], float dE);
      __host__ __device__ void setEnergy(float newEnergy);
      __host__ __device__ void setPreviousEnergy(float newPreviousEnergy);
      __host__ __device__ void setCurrentRegion(const RegionBaseT* reg);
      __host__ __device__ const ElementT* getScatteringElement() const;
      __host__ __device__ void setScatteringElement(const ElementT* scatteringElement);
      __host__ __device__ float getPhi() const;
      __host__ __device__ float getTheta() const;
      __host__ __device__ bool isTrajectoryComplete() const;
      __host__ __device__ void setTrajectoryComplete(bool trajectoryComplete);
      __host__ __device__ long getIdent() const;
      __host__ __device__ float getInitialEnergy() const;
      long getParentID() const;

   private:
      // The x,y & z coordinates of the electron
      double mPosition[3];

      // The location of the electron before the last call to updatePosition
      double mPrevPosition[3];

      // The direction of the current trajectory segment
      float mPhi, mTheta; // transient

      // The kinetic energy of the electron
      float mEnergy; // transient

      // Kinetic energy of the electron upon conclusion of the previous step
      float previousEnergy; // transient

      const float initialEnergy; // check if energy has increased, which is invalid

      int mStepCount; // transient

      RegionBaseT const * mCurrentRegion; // transient

      RegionBaseT const * mPrevRegion; // transient

      ElementT const * mScatteringElement; // transient

      bool mTrajectoryComplete; // transient

      long ident; // A unique identifying number to assist tracking, final

      long parentID; // 0 if from e-gun. Otherwise ID of parent.
   };

   extern long getlastIdent();
}

#endif