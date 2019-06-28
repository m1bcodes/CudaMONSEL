#ifndef _ALGORITHM_USER_CUH_
#define _ALGORITHM_USER_CUH_

#include "gov\nist\microanalysis\EPQLibrary\Strategy.cuh"
#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace AlgorithmUser
{
   class AlgorithmUser
   {
   protected:
      __host__ __device__ AlgorithmUser();

      void addDefaultAlgorithm(char cls[], const AlgorithmClassT* ac);
      const AlgorithmClassT* getAlgorithm(char clsName[]) const;

      //virtual void initializeDefaultStrategy() = 0;
      virtual void initializeDefaultStrategy() = 0;

   private:
      StrategyT mLocalOverride;
   };

   extern const BetheElectronEnergyLossT* getDefaultBetheEnergyLoss();

   extern StrategyT mGlobalOverride;

   extern const EdgeEnergyT& sDefaultEdgeEnergy;
   //extern TransitionEnergy sDefaultTransitionEnergy = null;
   //extern MassAbsorptionCoefficient sDefaultMAC = null;
   //extern FluorescenceYieldMean sDefaultFluorescenceYieldMean = null;
   //extern FluorescenceYield sDefaultFluorescenceYield = null;
   extern BetheElectronEnergyLossT* sDefaultBetheEnergyLoss;
   //extern Bremsstrahlung.AngularDistribution sDefaultAngularDistribution = null;
   //extern CorrectionAlgorithm sDefaultCorrectionAlgorithm = null;

   const EdgeEnergyT& getDefaultEdgeEnergy();
}

#endif