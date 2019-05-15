// file: gov\nist\microanalysis\EPQLibrary\BetheElectronEnergyLoss.cuh

#ifndef _BETHE_ELECTRON_ENERGY_LOSS_CUH_
#define _BETHE_ELECTRON_ENERGY_LOSS_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AlgorithmClass.cuh"

namespace BetheElectronEnergyLoss
{
   class BetheElectronEnergyLoss : public AlgorithmClassT
   {
      void initializeDefaultStrategy() override;

   public:
      BetheElectronEnergyLoss(StringT name, const ReferenceT& ref);

      AlgorithmClassT const * const * getAllImplementations() const override;

      virtual double compute(const ElementT& elm, double eB) const = 0;
   };

   class JoyLuoBetheElectronEnergyLoss : public BetheElectronEnergyLoss
   {
   public:
      JoyLuoBetheElectronEnergyLoss();
      double compute(const ElementT& el, double eB) const override;

   private:
      VectorXd mK;
      const double K;
   };

   class Bethe30ModElectronEnergyLoss : public BetheElectronEnergyLoss
   {
   public:
      Bethe30ModElectronEnergyLoss();
      double compute(const ElementT& el, double eB) const override;

   private:
      const double K;
   };

   class Bethe30ElectronEnergyLoss : public BetheElectronEnergyLoss
   {
   public:
      Bethe30ElectronEnergyLoss();
      double compute(const ElementT& el, double eB) const override;

   private:
      const double K;
   };

   class StragglingModified : public BetheElectronEnergyLoss
   {
   public:
      StragglingModified(const BetheElectronEnergyLoss& base, double percent);
      double compute(const ElementT& el, double eB) const override;

   private:
      const BetheElectronEnergyLoss& mBethe;
      const double mPercent;
   };

   extern const BetheElectronEnergyLoss& JoyLuo1989;
   extern const BetheElectronEnergyLoss& Bethe1930;
   extern const BetheElectronEnergyLoss& Bethe1930Strict;
}

#endif