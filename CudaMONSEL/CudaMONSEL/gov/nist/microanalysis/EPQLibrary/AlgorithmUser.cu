#include "gov\nist\microanalysis\EPQLibrary\AlgorithmUser.cuh"
#include "gov\nist\microanalysis\EPQLibrary\EdgeEnergy.cuh"
#include "gov\nist\microanalysis\EPQLibrary\BetheElectronEnergyLoss.cuh"

namespace AlgorithmUser
{
   StrategyT mGlobalOverride;

   const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::SuperSet;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::DTSA;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::Chantler2005;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::NISTxrtdb;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::Wernish84;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::DHSIonizationEnergy;

   //TransitionEnergy sDefaultTransitionEnergy = NULL;
   //MassAbsorptionCoefficient sDefaultMAC = NULL;
   //FluorescenceYieldMean sDefaultFluorescenceYieldMean = NULL;
   //FluorescenceYield sDefaultFluorescenceYield = NULL;
   BetheElectronEnergyLossT* sDefaultBetheEnergyLoss = NULL;
   //Bremsstrahlung.AngularDistribution sDefaultAngularDistribution = NULL;
   //CorrectionAlgorithm sDefaultCorrectionAlgorithm = NULL;

   // TODO: figure out what to do with this
   // https://stackoverflow.com/questions/8630160/call-to-pure-virtual-function-from-base-class-constructor
   AlgorithmUser::AlgorithmUser()
   {
      //initializeDefaultStrategy();
   }

   const EdgeEnergyT& getDefaultEdgeEnergy()
   {
      return sDefaultEdgeEnergy;
   }

   //static CorrectionAlgorithm getDefaultCorrectionAlgorithm()
   //{
   //   return sDefaultCorrectionAlgorithm == NULL ? CorrectionAlgorithm.XPPExtended : sDefaultCorrectionAlgorithm;
   //}

   // static TransitionEnergy getDefaultTransitionEnergy() {
   //   return sDefaultTransitionEnergy == NULL ? TransitionEnergy.SuperSet : sDefaultTransitionEnergy;
   //}

   // static MassAbsorptionCoefficient getDefaultMAC() {
   //   return sDefaultMAC == NULL ? MassAbsorptionCoefficient.Chantler2005 : sDefaultMAC;
   //}

   // static FluorescenceYieldMean getDefaultFluorescenceYieldMean() {
   //   return sDefaultFluorescenceYieldMean == NULL ? FluorescenceYieldMean.DefaultMean : sDefaultFluorescenceYieldMean;
   //}

   // static FluorescenceYield getDefaultFluorescenceYield() {
   //   return sDefaultFluorescenceYield == NULL ? FluorescenceYield.DefaultShell : sDefaultFluorescenceYield;
   //}

   const BetheElectronEnergyLossT* getDefaultBetheEnergyLoss()
   {
      return sDefaultBetheEnergyLoss == nullptr ? &BetheElectronEnergyLoss::JoyLuo1989 : sDefaultBetheEnergyLoss;
   }

   // static Bremsstrahlung.AngularDistribution getDefaultAngularDistribution() {
   //   return sDefaultAngularDistribution == NULL ? Bremsstrahlung.AngularDistribution.Isotropic : sDefaultAngularDistribution;
   //}

   //Strategy getActiveStrategy() {
   //   final Strategy res = mLocalOverride != NULL ? (Strategy)mLocalOverride.clone() : new Strategy();
   //   if (mGlobalOverride != NULL)
   //      res.apply(mGlobalOverride);
   //   return res;
   //}

   // AlgorithmClass getAlgorithm(Class< ? > cls) {
   //   AlgorithmClass res = NULL;
   //   if (mGlobalOverride != NULL)
   //      res = mGlobalOverride.getAlgorithm(cls);
   //   if ((res == NULL) && (mLocalOverride != NULL))
   //      res = mLocalOverride.getAlgorithm(cls);
   //   return res;
   //}

   const AlgorithmClassT* AlgorithmUser::getAlgorithm(char clsName[]) const
   {
      const AlgorithmClassT* res = nullptr;
      if (mGlobalOverride.empty())
         res = mGlobalOverride.getAlgorithm(clsName);
      if ((res == nullptr) && !(mLocalOverride.empty()))
         res = mLocalOverride.getAlgorithm(clsName);
      return res;
   }


   void clearGlobalOverride()
   {
      mGlobalOverride.clear();
      //sDefaultTransitionEnergy = NULL;
      //sDefaultEdgeEnergy = NULL;
      //sDefaultMAC = NULL;
      //sDefaultFluorescenceYield = NULL;
      //sDefaultFluorescenceYieldMean = NULL;
      sDefaultBetheEnergyLoss = NULL;
      //sDefaultAngularDistribution = NULL;
   }

   void applyGlobalOverride(StrategyT* strat)
   {
      if (strat != NULL) {
         mGlobalOverride.addAll(*strat);
         //sDefaultTransitionEnergy = (TransitionEnergy)strat.getAlgorithm(TransitionEnergy.class);
         //sDefaultEdgeEnergy = (EdgeEnergy)strat.getAlgorithm(EdgeEnergy.class);
         //sDefaultMAC = (MassAbsorptionCoefficient)strat.getAlgorithm(MassAbsorptionCoefficient.class);
         //sDefaultFluorescenceYield = (FluorescenceYield)strat.getAlgorithm(FluorescenceYield.class);
         //sDefaultFluorescenceYieldMean = (FluorescenceYieldMean)strat.getAlgorithm(FluorescenceYieldMean.class);
         sDefaultBetheEnergyLoss = (BetheElectronEnergyLossT*)strat->getAlgorithm("BetheElectronEnergyLoss");
         //sDefaultAngularDistribution = (Bremsstrahlung.AngularDistribution) strat.getAlgorithm(Bremsstrahlung.AngularDistribution.class);
         //sDefaultCorrectionAlgorithm = (CorrectionAlgorithm)strat.getAlgorithm(CorrectionAlgorithm.class);
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

   void AlgorithmUser::addDefaultAlgorithm(char cls[], AlgorithmClassT const * const ac)
   {
      if (ac == nullptr) printf("AlgorithmUser::addDefaultAlgorithm: NULL class");
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
   //   if (mLocalOverride != NULL) {
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
