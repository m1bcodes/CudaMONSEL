#include "gov\nist\microanalysis\EPQLibrary\AlgorithmUser.cuh"
#include "gov\nist\microanalysis\EPQLibrary\EdgeEnergy.cuh"

namespace AlgorithmUser
{
   // TODO: figure out what to do with this
   // https://stackoverflow.com/questions/8630160/call-to-pure-virtual-function-from-base-class-constructor
   AlgorithmUser::AlgorithmUser()
   {
      //initializeDefaultStrategy();
   }

   const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::SuperSet;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::DTSA;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::Chantler2005;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::NISTxrtdb;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::Wernish84;
   //const EdgeEnergyT& sDefaultEdgeEnergy = EdgeEnergy::DHSIonizationEnergy;

   const EdgeEnergyT& getDefaultEdgeEnergy()
   {
      return sDefaultEdgeEnergy;
   }

   //static CorrectionAlgorithm getDefaultCorrectionAlgorithm()
   //{
   //   return sDefaultCorrectionAlgorithm == null ? CorrectionAlgorithm.XPPExtended : sDefaultCorrectionAlgorithm;
   //}

   // static TransitionEnergy getDefaultTransitionEnergy() {
   //   return sDefaultTransitionEnergy == null ? TransitionEnergy.SuperSet : sDefaultTransitionEnergy;
   //}

   // static MassAbsorptionCoefficient getDefaultMAC() {
   //   return sDefaultMAC == null ? MassAbsorptionCoefficient.Chantler2005 : sDefaultMAC;
   //}

   // static FluorescenceYieldMean getDefaultFluorescenceYieldMean() {
   //   return sDefaultFluorescenceYieldMean == null ? FluorescenceYieldMean.DefaultMean : sDefaultFluorescenceYieldMean;
   //}

   // static FluorescenceYield getDefaultFluorescenceYield() {
   //   return sDefaultFluorescenceYield == null ? FluorescenceYield.DefaultShell : sDefaultFluorescenceYield;
   //}

   // static BetheElectronEnergyLoss getDefaultBetheEnergyLoss() {
   //   return sDefaultBetheEnergyLoss == null ? BetheElectronEnergyLoss.JoyLuo1989 : sDefaultBetheEnergyLoss;
   //}

   // static Bremsstrahlung.AngularDistribution getDefaultAngularDistribution() {
   //   return sDefaultAngularDistribution == null ? Bremsstrahlung.AngularDistribution.Isotropic : sDefaultAngularDistribution;
   //}

   //Strategy getActiveStrategy() {
   //   final Strategy res = mLocalOverride != null ? (Strategy)mLocalOverride.clone() : new Strategy();
   //   if (mGlobalOverride != null)
   //      res.apply(mGlobalOverride);
   //   return res;
   //}

   // AlgorithmClass getAlgorithm(Class< ? > cls) {
   //   AlgorithmClass res = null;
   //   if (mGlobalOverride != null)
   //      res = mGlobalOverride.getAlgorithm(cls);
   //   if ((res == null) && (mLocalOverride != null))
   //      res = mLocalOverride.getAlgorithm(cls);
   //   return res;
   //}

   // AlgorithmClass getAlgorithm(String clsName) {
   //   AlgorithmClass res = null;
   //   if (mGlobalOverride != null)
   //      res = mGlobalOverride.getAlgorithm(clsName);
   //   if ((res == null) && (mLocalOverride != null))
   //      res = mLocalOverride.getAlgorithm(clsName);
   //   return res;
   //}


   //static void applyGlobalOverride(Strategy strat)
   //{
   //   if (strat != null) {
   //      mGlobalOverride = strat;
   //      sDefaultTransitionEnergy = (TransitionEnergy)strat.getAlgorithm(TransitionEnergy.class);
   //      sDefaultEdgeEnergy = (EdgeEnergy)strat.getAlgorithm(EdgeEnergy.class);
   //      sDefaultMAC = (MassAbsorptionCoefficient)strat.getAlgorithm(MassAbsorptionCoefficient.class);
   //      sDefaultFluorescenceYield = (FluorescenceYield)strat.getAlgorithm(FluorescenceYield.class);
   //      sDefaultFluorescenceYieldMean = (FluorescenceYieldMean)strat.getAlgorithm(FluorescenceYieldMean.class);
   //      sDefaultBetheEnergyLoss = (BetheElectronEnergyLoss)strat.getAlgorithm(BetheElectronEnergyLoss.class);
   //      sDefaultAngularDistribution = (Bremsstrahlung.AngularDistribution) strat.getAlgorithm(Bremsstrahlung.AngularDistribution.class);
   //      sDefaultCorrectionAlgorithm = (CorrectionAlgorithm)strat.getAlgorithm(CorrectionAlgorithm.class);
   //   }
   //   else
   //      clearGlobalOverride();
   //}

   //static  Strategy getGlobalStrategy() {
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

   //static  void clearGlobalOverride() {
   //   mGlobalOverride = null;
   //   sDefaultTransitionEnergy = null;
   //   sDefaultEdgeEnergy = null;
   //   sDefaultMAC = null;
   //   sDefaultFluorescenceYield = null;
   //   sDefaultFluorescenceYieldMean = null;
   //   sDefaultBetheEnergyLoss = null;
   //   sDefaultAngularDistribution = null;
   //}

   //protected void addDefaultAlgorithm(Class< ? > cls, AlgorithmClass ac) {
   //   assert ac != null;
   //   if (mLocalOverride == null)
   //      mLocalOverride = new Strategy();
   //   mLocalOverride.addAlgorithm(cls, ac);
   //}

   // void documentStrategy(Writer wr)
   //   throws IOException{
   //   for (final String cls : mGlobalOverride.getStrategyMap().keySet()) {
   //      wr.write(getAlgorithm(cls).toString());
   //      wr.write('\n');
   //   }
   //}

   //   Strategy getEffectiveStrategyHelper() {
   //   final Strategy strat = new Strategy();
   //   if (mLocalOverride != null) {
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
