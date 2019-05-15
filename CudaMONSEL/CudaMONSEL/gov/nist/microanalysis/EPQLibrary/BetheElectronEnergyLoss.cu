#include "gov\nist\microanalysis\EPQLibrary\BetheElectronEnergyLoss.cuh"

#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Reference.cuh"
#include "gov\nist\microanalysis\EPQLibrary\FromSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\MeanIonizationPotential.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

namespace BetheElectronEnergyLoss
{
   BetheElectronEnergyLoss::BetheElectronEnergyLoss(StringT name, const ReferenceT& ref) : AlgorithmClass("Stopping power", name, ref)
   {
   }

   JoyLuoBetheElectronEnergyLoss::JoyLuoBetheElectronEnergyLoss() : BetheElectronEnergyLoss("Joy-Luo", Reference::GoldsteinBook), K(-(785 * ToSI::EV * ::pow(ToSI::CM, 3.0)) / (ToSI::ANGSTROM * ToSI::GRAM))
   {
      mK.resize(Element::elmEndOfElements);
      for (int z = Element::elmH; z < Element::elmEndOfElements; ++z)
         mK[z] = 0.731 + 0.0688 * ::log10(z);
   }

   double JoyLuoBetheElectronEnergyLoss::compute(const ElementT& el, double eB) const
   {
      const int z = el.getAtomicNumber();
      const MeanIonizationPotentialT* mip = (MeanIonizationPotentialT*)getAlgorithm("MeanIonizationPotential");
      const double j = mip->compute(el);
      const double j_star = j / (1.0 + (mK[z] * j / eB));
      return ((K * z) / (el.getAtomicWeight() * FromSI::eV(eB))) * ::log((1.166 * eB) / j_star);
   }

   const JoyLuoBetheElectronEnergyLoss JoyLuo1989Ref;
   const BetheElectronEnergyLoss& JoyLuo1989 = JoyLuo1989Ref;

   /**
   * Bethe1930 - The original expression of Bethe for the stopping power
   * adjusted so that even when the electron energy falls below about J/1.166,
   * the electron continues to decelerate (albeit slowly).
   */

   Reference::CrudeReference Leipzig1930CR("Bethe H. Ann. Phys. (Leipzig) 1930; 5: 325");
   Bethe30ModElectronEnergyLoss::Bethe30ModElectronEnergyLoss() : BetheElectronEnergyLoss("Bethe(Modified)", Leipzig1930CR), K(-(785 * ToSI::EV * ::pow(ToSI::CM, 3.0)) / (ToSI::ANGSTROM * ToSI::GRAM))
   {
   }

   double Bethe30ModElectronEnergyLoss::compute(const ElementT& el, double eB) const
   {
      const MeanIonizationPotentialT* mip = (MeanIonizationPotentialT*)getAlgorithm("MeanIonizationPotential");
      const double e_eV = FromSI::eV(eB);
      const double j = FromSI::eV(mip->compute(el));
      double f = 1.166 * e_eV / j;
      // The low energy modification...
      if (f < 1.1)
         f = 1.1;
      return ((K * el.getAtomicNumber()) / (el.getAtomicWeight() * e_eV)) * ::log(f);
   }

   const Bethe30ModElectronEnergyLoss Bethe1930Ref;
   const BetheElectronEnergyLoss& Bethe1930 = Bethe1930Ref;

   /**
   * Bethe1930Strict - The original expression of Bethe for the stopping power.
   * Below eB = J/1.166 the energy loss goes positive (the electron shows
   * unphysical acceleration.)
   */

   Bethe30ElectronEnergyLoss::Bethe30ElectronEnergyLoss() : BetheElectronEnergyLoss("Bethe(Modified)", Leipzig1930CR), K(-(785 * ToSI::EV * ::pow(ToSI::CM, 3.0)) / (ToSI::ANGSTROM * ToSI::GRAM))
   {
   }

   double Bethe30ElectronEnergyLoss::compute(const ElementT& el, double eB) const
   {
      const MeanIonizationPotentialT* mip = (MeanIonizationPotentialT*)getAlgorithm("MeanIonizationPotential");
      const double e_eV = FromSI::eV(eB);
      const double j = FromSI::eV(mip->compute(el));
      const double f = 1.166 * e_eV / j;
      return ((K * el.getAtomicNumber()) / (el.getAtomicWeight() * e_eV)) * ::log(f);
   }

   const Bethe30ElectronEnergyLoss Bethe1930StrictRef;
   const BetheElectronEnergyLoss& Bethe1930Strict = Bethe1930StrictRef;

   /**
   * <p>
   * Modifies an existing BetheElectronEnergyLoss model to add variation in the
   * amount of energy lost per step. The class takes the nominal energy loss
   * and modifies it by a certain randomized fractional amount to emulate the
   * way sometimes an electron may loose slightly more than average or slightly
   * less than average.
   * </p>
   * <p>
   * Copyright: Pursuant to title 17 Section 105 of the United States Code this
   * software is not subject to copyright protection and is in the public
   * domain
   * </p>
   * <p>
   * Institution: National Institute of Standards and Technology
   * </p>
   *
   * @author nritchie
   * @version 1.0
   */
   StragglingModified::StragglingModified(const BetheElectronEnergyLoss& base, double percent) : 
      BetheElectronEnergyLoss("Straggling[" + base.getName() + "]", base.getReferenceObj()),
      mBethe(base),
      mPercent(percent)
   {
   }

   double StragglingModified::compute(const ElementT& elm, double eB) const
   {
      const double bee = mBethe.compute(elm, eB);
      return ::fmin(0.0, bee * (1.0 + Math2::generateGaussianNoise(0, 1) * mPercent));
   }

   static const AlgorithmClassT* mAllImplementations[] = {
      &Bethe1930,
      &Bethe1930Strict,
      &JoyLuo1989
   };

   AlgorithmClassT const * const * BetheElectronEnergyLoss::getAllImplementations() const
   {
      return mAllImplementations;
   }

   void BetheElectronEnergyLoss::initializeDefaultStrategy()
   {
      addDefaultAlgorithm("MeanIonizationPotential", &MeanIonizationPotential::Berger83);
   }
}