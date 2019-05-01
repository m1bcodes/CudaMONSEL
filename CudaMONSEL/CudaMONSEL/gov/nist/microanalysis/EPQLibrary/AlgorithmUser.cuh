#ifndef _ALGORITHM_USER_CUH_
#define _ALGORITHM_USER_CUH_

namespace AlgorithmUser
{
   class AlgorithmUser
   {
   protected:
      AlgorithmUser();

      virtual void initializeDefaultStrategy() = 0;
      void initializeDefaultStrategyWrapper(); // https://stackoverflow.com/questions/8630160/call-to-pure-virtual-function-from-base-class-constructor

   //private:
      //Strategy mLocalOverride = null;
      //static Strategy mGlobalOverride = null;

      //static EdgeEnergy sDefaultEdgeEnergy = NULL;
      //static TransitionEnergy sDefaultTransitionEnergy = null;
      //static MassAbsorptionCoefficient sDefaultMAC = null;
      //static FluorescenceYieldMean sDefaultFluorescenceYieldMean = null;
      //static FluorescenceYield sDefaultFluorescenceYield = null;
      //static BetheElectronEnergyLoss sDefaultBetheEnergyLoss = null;
      //static Bremsstrahlung.AngularDistribution sDefaultAngularDistribution = null;
      //static CorrectionAlgorithm sDefaultCorrectionAlgorithm = null;
   };

   //EdgeEnergy getDefaultEdgeEnergy();
}

#endif