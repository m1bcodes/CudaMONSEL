#ifndef _ELECTRON_CUH_
#define _ELECTRON_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace Electron
{
   class Electron
   {
   public:
      __host__ __device__ Electron(const double initialPos[], double kE);
      Electron(const double initialPos[], double theta, double phi, double kE);
      __host__ __device__ Electron(const Electron& parent, double theta, double phi, double kE);

      __host__ __device__ void init(const double initialPos[], double theta, double phi, double kE);

      __host__ __device__ void setDirection(double theta, double phi);
      __host__ __device__ const double * getPosition() const;
      __host__ __device__ void setPosition(const double newpos[]);
      __host__ __device__ const double * getPrevPosition() const;
      __host__ __device__ const RegionBaseT* getCurrentRegion() const;
      const RegionBaseT* getPreviousRegion() const;
      __host__ __device__ double getEnergy() const;
      __host__ __device__ double getPreviousEnergy() const;
      __host__ __device__ int getStepCount() const;
      double stepLength() const;
      __host__ __device__ void candidatePoint(const double dS, double res[]) const;
      __host__ __device__ void updateDirection(const double dTheta, const double dPhi);
      __host__ __device__ void move(const double newPoint[], double dE);
      __host__ __device__ void setEnergy(double newEnergy);
      __host__ __device__ void setPreviousEnergy(double newPreviousEnergy);
      __host__ __device__ void setCurrentRegion(const RegionBaseT* reg);
      __host__ __device__ const ElementT* getScatteringElement() const;
      __host__ __device__ void setScatteringElement(const ElementT* scatteringElement);
      __host__ __device__ double getPhi() const;
      __host__ __device__ double getTheta() const;
      __host__ __device__ bool isTrajectoryComplete() const;
      __host__ __device__ void setTrajectoryComplete(bool trajectoryComplete);
      __host__ __device__ long getIdent() const;
      long getParentID() const;

   private:
      // The x,y & z coordinates of the electron
      double mPosition[3];

      // The location of the electron before the last call to updatePosition
      double mPrevPosition[3];

      // The direction of the current trajectory segment
      double mPhi, mTheta; // transient

      // The kinetic energy of the electron
      double mEnergy; // transient

      // Kinetic energy of the electron upon conclusion of the previous step
      double previousEnergy; // transient

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