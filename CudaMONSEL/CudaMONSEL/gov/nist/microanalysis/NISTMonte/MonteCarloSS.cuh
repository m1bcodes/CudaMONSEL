#ifndef _MONTE_CARLO_SS_CU_
#define _MONTE_CARLO_SS_CU_

#include "gov\nist\microanalysis\NISTMonte\NullMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

#include "Amphibian\vector.cuh"
#include "Amphibian\stack.cuh"

namespace MonteCarloSS
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   extern __constant__ const int ScatterEvent;
   extern __constant__ const int NonScatterEvent;
   extern __constant__ const int BackscatterEvent;
   extern __constant__ const int ExitMaterialEvent;
   extern __constant__ const int TrajectoryStartEvent;
   extern __constant__ const int TrajectoryEndEvent;
   extern __constant__ const int LastTrajectoryEvent;
   extern __constant__ const int FirstTrajectoryEvent;
   extern __constant__ const int StartSecondaryEvent;
   extern __constant__ const int EndSecondaryEvent;
   extern __constant__ const int PostScatterEvent;

   extern __constant__ const int  BeamEnergyChanged;

   extern __constant__ const float ChamberRadius;
   extern __constant__ const float SMALL_DISP;
#else
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

   extern const float ChamberRadius;
   extern const float SMALL_DISP;
#endif

   class MonteCarloSS final
   {
      typedef amp::stack<ElectronT*> ElectronStack;
      typedef amp::vector<ActionListenerT*> ActionListeners;

   public:
      __host__ __device__ MonteCarloSS(ElectronGunT const * gun, RegionT * mChamber, ElectronT * electron);
      RegionT* getChamber();
      __host__ __device__ void initializeTrajectory();

      __host__ __device__ void takeStep();

      __host__ __device__ void trackSecondaryElectron(ElectronT* newElectron);
      __host__ __device__ bool allElectronsComplete();
      __host__ __device__ void runTrajectory();
      //int getElectronGeneration() const;
      __host__ __device__ void addActionListener(ActionListenerT& sel);
      __host__ __device__ void removeActionListener(ActionListenerT& sel);
      __host__ __device__ void runMultipleTrajectories(int n);
      __host__ __device__ double getBeamEnergy() const;
      void setElectronGun(ElectronGunT& gun);
      //void setBeamEnergy(double beamEnergy);
      void computeDetectorPosition(double elevation, double theta, double res[]);
      __host__ __device__ const ElectronT& getElectron() const;
      Element::UnorderedSetT getElementSet() const;
      void rotate(double pivot[], double phi, double theta, double psi);
      void translate(double distance[]);
      const RegionBaseT* findRegionContaining(double point[]) const;

   private:
      __host__ __device__ void fireEvent(const int);
      
      RegionT * mChamber;
      ElectronGunT const * mGun;
      ElectronT * mElectron;
      ElectronStack mElectronStack;

      ActionListeners mEventListeners;
   };
}

#endif
