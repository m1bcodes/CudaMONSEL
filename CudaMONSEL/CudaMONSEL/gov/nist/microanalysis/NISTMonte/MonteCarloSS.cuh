#ifndef _MONTE_CARLO_SS_CU_
#define _MONTE_CARLO_SS_CU_

#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\NISTMonte\NullMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

namespace MonteCarloSS
{
   extern const double ChamberRadius;
   extern const double SMALL_DISP;

   class MonteCarloSS final
   {
      typedef std::stack<ElectronT*> ElectronStack;

   public:
      MonteCarloSS(ElectronGunT const * gun, RegionT * mChamber, ElectronT * electron);
      RegionT* getChamber();
      void initializeTrajectory();

      void takeStep();

      void trackSecondaryElectron(ElectronT* newElectron);
      bool allElectronsComplete();
      void runTrajectory();
      void runMultipleTrajectories(int n);
      double getBeamEnergy();
      void setElectronGun(ElectronGunT& gun);
      //void setBeamEnergy(double beamEnergy);
      VectorXd computeDetectorPosition(double elevation, double theta);
      Element::UnorderedSetT getElementSet() const;
      void rotate(double pivot[], double phi, double theta, double psi);
      void translate(double distance[]);
      const RegionBaseT* findRegionContaining(double point[]) const;

   private:
      RegionT * mChamber;
      ElectronGunT const * mGun;
      ElectronT * mElectron;
      ElectronStack mElectronStack;

      //virtual int GetId() = 0;
   };

   double dist(const double pos0[], const double pos1[]);

   const NullMaterialScatterModelT NULL_MSM;
}

#endif
