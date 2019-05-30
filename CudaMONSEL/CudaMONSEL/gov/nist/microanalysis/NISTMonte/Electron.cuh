#ifndef _ELECTRON_CUH_
#define _ELECTRON_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace Electron
{
   class Electron
   {
   public:
      Electron(const double initialPos[], double kE);
      Electron(const double initialPos[], double theta, double phi, double kE);
      Electron(const Electron& parent, double theta, double phi, double kE);

      void Init(const double initialPos[], double theta, double phi, double kE);

      void setDirection(double theta, double phi);
      __host__ __device__ const double * getPosition() const;
      void setPosition(const double newpos[]);
      __host__ __device__ const double * getPrevPosition() const;
      const RegionBaseT* getCurrentRegion() const;
      const RegionBaseT* getPreviousRegion() const;
      double getEnergy() const;
      double getPreviousEnergy() const;
      int getStepCount() const;
      double stepLength() const;
      __host__ __device__ void candidatePoint(double dS, double res[]) const;
      void updateDirection(double dTheta, double dPhi);
      void move(const double newPoint[], double dE);
      void setEnergy(double newEnergy);
      void setPreviousEnergy(double newPreviousEnergy);
      __host__ __device__ void setCurrentRegion(const RegionBaseT* reg);
      const ElementT* getScatteringElement() const;
      void setScatteringElement(const ElementT* scatteringElement);
      double getPhi() const;
      double getTheta() const;
      bool isTrajectoryComplete() const;
      __host__ __device__ void setTrajectoryComplete(bool trajectoryComplete);
      long getIdent() const;
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