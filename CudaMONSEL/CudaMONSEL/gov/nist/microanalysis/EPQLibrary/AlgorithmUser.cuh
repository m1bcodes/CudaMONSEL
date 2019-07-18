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

      __host__ __device__ void addDefaultAlgorithm(char cls[], AlgorithmClassT const * ac);
      __host__ __device__ const AlgorithmClassT* getAlgorithm(char clsName[]) const;

      //virtual void initializeDefaultStrategy() = 0;
      __host__ __device__ virtual void initializeDefaultStrategy() = 0;

   private:
      StrategyT mLocalOverride;
   };

   __host__ __device__ extern const BetheElectronEnergyLossT* getDefaultBetheEnergyLoss();

   extern StrategyT mGlobalOverride;
   extern __device__ StrategyT* d_mGlobalOverride;

   extern const EdgeEnergyT& sDefaultEdgeEnergy;
   //extern TransitionEnergy sDefaultTransitionEnergy;
   //extern MassAbsorptionCoefficient sDefaultMAC;
   //extern FluorescenceYieldMean sDefaultFluorescenceYieldMean;
   //extern FluorescenceYield sDefaultFluorescenceYield;
   extern BetheElectronEnergyLossT* sDefaultBetheEnergyLoss;
   //extern Bremsstrahlung.AngularDistribution sDefaultAngularDistribution;
   //extern CorrectionAlgorithm sDefaultCorrectionAlgorithm;

   const EdgeEnergyT& getDefaultEdgeEnergy();

   extern __global__ void initCuda();
}

#endif