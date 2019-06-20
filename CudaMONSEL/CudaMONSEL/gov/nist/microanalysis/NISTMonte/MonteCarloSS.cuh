#ifndef _MONTE_CARLO_SS_CU_
#define _MONTE_CARLO_SS_CU_

#include "gov\nist\microanalysis\NISTMonte\NullMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

#include "Amphibian\vector.cuh"
#include "Amphibian\stack.cuh"

namespace MonteCarloSS
{
   extern const int ScatterEvent;
   extern const int NonScatterEvent;
   extern const int BackscatterEvent;
   extern const int ExitMaterialEvent;
   extern const int TrajectoryStartEvent;
   extern const int TrajectoryEndEvent;
   extern const int LastTrajectoryEvent;
   extern const int FirstTrajectoryEvent;
   extern const int StartSecondaryEvent;
   extern const int EndSecondaryEvent;
   extern const int PostScatterEvent;

   extern const int BeamEnergyChanged;

   extern const double ChamberRadius;
   extern const double SMALL_DISP;

   class MonteCarloSS final
   {
      typedef amp::stack<ElectronT*> ElectronStack;
      typedef amp::vector<ActionListenerT*> ActionListeners;

   public:
      MonteCarloSS(ElectronGunT const * gun, RegionT * mChamber, ElectronT * electron);
      RegionT* getChamber();
      void initializeTrajectory();

      void takeStep();

      void trackSecondaryElectron(ElectronT* newElectron);
      bool allElectronsComplete();
      void runTrajectory();
      //int getElectronGeneration() const;
      void addActionListener(ActionListenerT& sel);
      void removeActionListener(ActionListenerT& sel);
      void runMultipleTrajectories(int n);
      double getBeamEnergy() const;
      void setElectronGun(ElectronGunT& gun);
      //void setBeamEnergy(double beamEnergy);
      void computeDetectorPosition(double elevation, double theta, double res[]);
      const ElectronT& getElectron() const;
      Element::UnorderedSetT getElementSet() const;
      void rotate(double pivot[], double phi, double theta, double psi);
      void translate(double distance[]);
      const RegionBaseT* findRegionContaining(double point[]) const;

   private:
      void fireEvent(const int);
      
      RegionT * mChamber;
      ElectronGunT const * mGun;
      ElectronT * mElectron;
      ElectronStack mElectronStack;

      ActionListeners mEventListeners;
   };
}

#endif
