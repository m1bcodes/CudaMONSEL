#ifndef _ELECTRON_CUH_
#define _ELECTRON_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace Electron
{
   class Electron
   {
   public:
      Electron(double initialPos[], double kE);
      Electron(double initialPos[], double theta, double phi, double kE);
      Electron(const Electron& parent, double theta, double phi, double kE);

      void Init(double initialPos[], double theta, double phi, double kE);

      void setDirection(double theta, double phi);
      VectorXd getPosition() const;
      void setPosition(double newpos[]);
      VectorXd getPrevPosition() const;
      const RegionBaseT* getCurrentRegion() const;
      const RegionBaseT* getPreviousRegion() const;
      double getEnergy() const;
      double getPreviousEnergy() const;
      int getStepCount() const;
      double stepLength() const;
      VectorXd candidatePoint(double dS) const;
      void updateDirection(double dTheta, double dPhi);
      void move(double newPoint[], double dE);
      void setEnergy(double newEnergy);
      void setPreviousEnergy(double newPreviousEnergy);
      void setCurrentRegion(const RegionBaseT* reg);
      const ElementT* getScatteringElement() const;
      void setScatteringElement(const ElementT* scatteringElement);
      double getPhi() const;
      double getTheta() const;
      bool isTrajectoryComplete() const;
      void setTrajectoryComplete(bool trajectoryComplete);
      long getIdent() const;
      long getParentID() const;

   private:
      // The x,y & z coordinates of the electron
      VectorXd mPosition; // final transient

      // The location of the electron before the last call to updatePosition
      VectorXd mPrevPosition; // final transient

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

      long parentID = 0; // 0 if from e-gun. Otherwise ID of parent.
   };

   extern long getlastIdent();

   extern double DefaultPos[];
   extern Electron Default;
}

#endif