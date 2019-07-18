#include "gov\nist\microanalysis\EPQLibrary\AlgorithmUser.cuh"
#include "gov\nist\microanalysis\EPQLibrary\EdgeEnergy.cuh"
#include "gov\nist\microanalysis\EPQLibrary\BetheElectronEnergyLoss.cuh"

namespace AlgorithmUser
{
   StrategyT mGlobalOverride;
   __device__ StrategyT* d_mGlobalOverride = nullptr;

   __global__ void initCuda()
   {
      d_mGlobalOverride = new StrategyT();
   }

   const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::SuperSet;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::DTSA;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::Chantler2005;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::NISTxrtdb;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::Wernish84;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::DHSIonizationEnergy;

   //TransitionEnergy sDefaultTransitionEnergy = nullptr;
   //MassAbsorptionCoefficient sDefaultMAC = nullptr;
   //FluorescenceYieldMean sDefaultFluorescenceYieldMean = nullptr;
   //FluorescenceYield sDefaultFluorescenceYield = nullptr;
   BetheElectronEnergyLossT* sDefaultBetheEnergyLoss = nullptr;
   //Bremsstrahlung.AngularDistribution sDefaultAngularDistribution = nullptr;
   //CorrectionAlgorithm sDefaultCorrectionAlgorithm = nullptr;

   __device__ BetheElectronEnergyLossT* d_sDefaultBetheEnergyLoss = nullptr;

   // TODO: figure out what to do with this, only need to call when it is overriden
   // https://stackoverflow.com/questions/8630160/call-to-pure-virtual-function-from-base-class-constructor
   __host__ __device__ AlgorithmUser::AlgorithmUser()
   {
      //initializeDefaultStrategy();
   }

   const EdgeEnergyT& getDefaultEdgeEnergy()
   {
      return sDefaultEdgeEnergy;
   }

   //static CorrectionAlgorithm getDefaultCorrectionAlgorithm()
   //{
   //   return sDefaultCorrectionAlgorithm == nullptr ? CorrectionAlgorithm.XPPExtended : sDefaultCorrectionAlgorithm;
   //}

   // static TransitionEnergy getDefaultTransitionEnergy() {
   //   return sDefaultTransitionEnergy == nullptr ? TransitionEnergy.SuperSet : sDefaultTransitionEnergy;
   //}

   // static MassAbsorptionCoefficient getDefaultMAC() {
   //   return sDefaultMAC == nullptr ? MassAbsorptionCoefficient.Chantler2005 : sDefaultMAC;
   //}

   // static FluorescenceYieldMean getDefaultFluorescenceYieldMean() {
   //   return sDefaultFluorescenceYieldMean == nullptr ? FluorescenceYieldMean.DefaultMean : sDefaultFluorescenceYieldMean;
   //}

   // static FluorescenceYield getDefaultFluorescenceYield() {
   //   return sDefaultFluorescenceYield == nullptr ? FluorescenceYield.DefaultShell : sDefaultFluorescenceYield;
   //}

   __host__ __device__ const BetheElectronEnergyLossT* getDefaultBetheEnergyLoss()
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      return d_sDefaultBetheEnergyLoss == nullptr ? BetheElectronEnergyLoss::d_JoyLuo1989 : d_sDefaultBetheEnergyLoss;
#else
      return sDefaultBetheEnergyLoss == nullptr ? &BetheElectronEnergyLoss::JoyLuo1989 : sDefaultBetheEnergyLoss;
#endif
   }

   // static Bremsstrahlung.AngularDistribution getDefaultAngularDistribution() {
   //   return sDefaultAngularDistribution == nullptr ? Bremsstrahlung.AngularDistribution.Isotropic : sDefaultAngularDistribution;
   //}

   //Strategy getActiveStrategy() {
   //   final Strategy res = mLocalOverride != nullptr ? (Strategy)mLocalOverride.clone() : new Strategy();
   //   if (mGlobalOverride != nullptr)
   //      res.apply(mGlobalOverride);
   //   return res;
   //}

   // AlgorithmClass getAlgorithm(Class< ? > cls) {
   //   AlgorithmClass res = nullptr;
   //   if (mGlobalOverride != nullptr)
   //      res = mGlobalOverride.getAlgorithm(cls);
   //   if ((res == nullptr) && (mLocalOverride != nullptr))
   //      res = mLocalOverride.getAlgorithm(cls);
   //   return res;
   //}

   __host__ __device__ const AlgorithmClassT* AlgorithmUser::getAlgorithm(char clsName[]) const
   {
      const AlgorithmClassT* res = nullptr;
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      if (d_mGlobalOverride->empty())
         res = d_mGlobalOverride->getAlgorithm(clsName);
#else
      if (mGlobalOverride.empty())
         res = mGlobalOverride.getAlgorithm(clsName);
#endif
      if ((res == nullptr) && !(mLocalOverride.empty()))
         res = mLocalOverride.getAlgorithm(clsName);
      return res;
   }


   void clearGlobalOverride()
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      d_mGlobalOverride->clear();
      //sDefaultTransitionEnergy = nullptr;
      //sDefaultEdgeEnergy = nullptr;
      //sDefaultMAC = nullptr;
      //sDefaultFluorescenceYield = nullptr;
      //sDefaultFluorescenceYieldMean = nullptr;
      d_sDefaultBetheEnergyLoss = nullptr;
      //sDefaultAngularDistribution = nullptr;
#else
      mGlobalOverride.clear();
      //sDefaultTransitionEnergy = nullptr;
      //sDefaultEdgeEnergy = nullptr;
      //sDefaultMAC = nullptr;
      //sDefaultFluorescenceYield = nullptr;
      //sDefaultFluorescenceYieldMean = nullptr;
      sDefaultBetheEnergyLoss = nullptr;
      //sDefaultAngularDistribution = nullptr;
#endif
   }

   void applyGlobalOverride(StrategyT* strat)
   {
      if (strat != nullptr) {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
         d_mGlobalOverride->addAll(*strat);
         //sDefaultTransitionEnergy = (TransitionEnergy)strat.getAlgorithm(TransitionEnergy.class);
         //sDefaultEdgeEnergy = (EdgeEnergy)strat.getAlgorithm(EdgeEnergy.class);
         //sDefaultMAC = (MassAbsorptionCoefficient)strat.getAlgorithm(MassAbsorptionCoefficient.class);
         //sDefaultFluorescenceYield = (FluorescenceYield)strat.getAlgorithm(FluorescenceYield.class);
         //sDefaultFluorescenceYieldMean = (FluorescenceYieldMean)strat.getAlgorithm(FluorescenceYieldMean.class);
         d_sDefaultBetheEnergyLoss = (BetheElectronEnergyLossT*)strat->getAlgorithm("BetheElectronEnergyLoss");
         //sDefaultAngularDistribution = (Bremsstrahlung.AngularDistribution) strat.getAlgorithm(Bremsstrahlung.AngularDistribution.class);
         //sDefaultCorrectionAlgorithm = (CorrectionAlgorithm)strat.getAlgorithm(CorrectionAlgorithm.class);
#else
         mGlobalOverride.addAll(*strat);
         //sDefaultTransitionEnergy = (TransitionEnergy)strat.getAlgorithm(TransitionEnergy.class);
         //sDefaultEdgeEnergy = (EdgeEnergy)strat.getAlgorithm(EdgeEnergy.class);
         //sDefaultMAC = (MassAbsorptionCoefficient)strat.getAlgorithm(MassAbsorptionCoefficient.class);
         //sDefaultFluorescenceYield = (FluorescenceYield)strat.getAlgorithm(FluorescenceYield.class);
         //sDefaultFluorescenceYieldMean = (FluorescenceYieldMean)strat.getAlgorithm(FluorescenceYieldMean.class);
         sDefaultBetheEnergyLoss = (BetheElectronEnergyLossT*)strat->getAlgorithm("BetheElectronEnergyLoss");
         //sDefaultAngularDistribution = (Bremsstrahlung.AngularDistribution) strat.getAlgorithm(Bremsstrahlung.AngularDistribution.class);
         //sDefaultCorrectionAlgorithm = (CorrectionAlgorithm)strat.getAlgorithm(CorrectionAlgorithm.class);
#endif
      }
      else
         clearGlobalOverride();
   }

   //static Strategy getGlobalStrategy() {
   //   Strategy res = new Strategy();
   //   res.addAll(mGlobalOverride);
   //   res.addAlgorithm(TransitionEnergy.class, getDefaultTransitionEnergy());
   //   res.addAlgorithm(EdgeEnergy.class, getDefaultEdgeEnergy());
   //   res.addAlgorithm(MassAbsorptionCoefficient.class, getDefaultMAC());
   //   res.addAlgorithm(FluorescenceYield.class, getDefaultFluorescenceYield());
   //   res.addAlgorithm(FluorescenceYieldMean.class, getDefaultFluorescenceYieldMean());
   //   res.addAlgorithm(BetheElectronEnergyLoss.class, getDefaultBetheEnergyLoss());
   //   res.addAlgorithm(Bremsstrahlung.AngularDistribution.class, getDefaultAngularDistribution());
   //   res.addAlgorithm(CorrectionAlgorithm.class, getDefaultCorrectionAlgorithm());
   //   return res;
   //}

   __host__ __device__ void AlgorithmUser::addDefaultAlgorithm(char cls[], AlgorithmClassT const * ac)
   {
      if (ac == nullptr) printf("AlgorithmUser::addDefaultAlgorithm: nullptr class");
      mLocalOverride.addAlgorithm(cls, ac);
   }

   // void documentStrategy(Writer wr)
   //   throws IOException{
   //   for (final String cls : mGlobalOverride.getStrategyMap().keySet()) {
   //      wr.write(getAlgorithm(cls).toString());
   //      wr.write('\n');
   //   }
   //}

   //   Strategy getEffectiveStrategyHelper() {
   //   final Strategy strat = new Strategy();
   //   if (mLocalOverride != nullptr) {
   //      strat.addAll(mLocalOverride);
   //      final Collection<AlgorithmClass> algs = mLocalOverride.getAlgorithms();
   //      for (final AlgorithmClass ac : algs)
   //         strat.addAll(((AlgorithmUser)ac).getEffectiveStrategyHelper());
   //   }
   //   return strat;
   //}

   // Strategy getEffectiveStrategy() {
   //   final Strategy strat = getEffectiveStrategyHelper();
   //   strat.apply(mGlobalOverride);
   //   return strat;
   //}
}
