#include "gov\nist\nanoscalemetrology\JMONSEL\JoyLuoNieminenCSD.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\EPQLibrary\MeanIonizationPotential.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\SEmaterial.cuh"

namespace JoyLuoNieminenCSD
{
   JoyLuoNieminenCSD::JoyLuoNieminenCSD(SEmaterialT& mat, double breakE) : mat(mat)
   {
      if (breakE > (mat.getWorkfunction() + ToSI::eV(1.)))
         this->breakE = breakE;
      else
         printf("JoyLuoNieminenCSD::JoyLuoNieminenCSD: Supplied breakpoint energy is too small.");
      setMaterial(&mat);
   }

   JoyLuoNieminenCSD::JoyLuoNieminenCSD(MaterialT& mat, double bh, double breakE) : mat(mat)
   {
      if (breakE > (bh + ToSI::eV(1.)))
         this->breakE = breakE;
      else
         printf("JoyLuoNieminenCSD::JoyLuoNieminenCSD: Supplied breakpoint energy is too small.");
      setMaterial(mat, bh);
   }

   void JoyLuoNieminenCSD::setMaterial(const MaterialT& mat, double bh)
   {
      this->mat = mat;
      this->bh = bh;
      init();
   }

   void JoyLuoNieminenCSD::setMaterial(const SEmaterialT* mat)
   {
      this->mat = *mat;
      bh = -mat->getEnergyCBbottom();
      init();
   }

   void JoyLuoNieminenCSD::init()
   {
      nce = mat.getElementCount();
      if (nce == 0)
         return;
      recipJ.resize(nce);
      coef.resize(nce);
      beta.resize(nce);
      bhplus1eV = bh + ToSI::eV(1.);
      if (breakE < bhplus1eV)
         breakE = bhplus1eV;
      //final Object[] el = mat.getElementSet().toArray();
      
      /* Why can't I cast the above array to Element[]? */
      int i = 0;
      for (auto el : mat.getElementSet()) {
         recipJ[i] = 1.166 / (el->getAtomicNumber() < 13 ? MeanIonizationPotential::Wilson41.compute(*el) : MeanIonizationPotential::Sternheimer64.compute(*el));
         beta[i] = 1. - (recipJ[i] * bhplus1eV);
         /*
         * The constant in the following expression is in units of J^2 m^2/kg.
         * It is appropriate for density in kg/m^3, energy in Joules, and
         * resulting continuous energy loss in J/m.
         */
         coef[i] = 2.01507E-28 * mat.getDensity() * mat.weightFraction(*el, true) * el->getAtomicNumber() / el->getAtomicWeight();
      }
      /*
      * In the original MONSEL, the CSD routine knows the minEforTracking and
      * can use it to optimize--quitting early if the electron energy falls
      * below this during a step (since the electron will be dropped anyway).
      * In this version, minEforTracking is just a synonym for the barrier
      * height.
      */
      minEforTracking = bh;

      // Determine the proportionality constant
      gamma = 0.;
      for (int i = 0; i < nce; i++)
         gamma += coef[i] * ::log((recipJ[i] * breakE) + beta[i]);
      gamma /= ::pow(breakE, 3.5);
   }

   static double maxlossfraction = 0.1;

   double JoyLuoNieminenCSD::compute(double len, const ElectronT* pe) const
   {
      return compute(len, pe->getEnergy());
   }

   double JoyLuoNieminenCSD::compute(double len, double kE) const
   {
      /*
      * No energy loss in vacuum, and don't waste time on any with energy
      * already low enough to be dropped.
      */
      if ((nce == 0) || (kE < minEforTracking) || (kE <= 0.))
         return 0.;
      /*
      * For energies below the break use the Nieminen formula, which can be
      * integrated exactly to give the energy loss for distance len.
      */
      if (kE <= breakE)
         return (kE / ::pow(1. + (1.5 * gamma * len * kE * ::sqrt(kE)), 2. / 3.)) - kE;

      // Otherwise, Joy/Luo formula
      double loss = 0.;
      for (int i = 0; i < nce; i++)
         loss += coef[i] * ::log((recipJ[i] * kE) + beta[i]);
      loss *= len / kE;

      if (loss <= (maxlossfraction * kE))
         return -loss;

      /*
      * We're here if the loss was too big to take all in one step. Recursively
      * divide the step.
      */
      loss = compute(len / 2., kE); // loss from 1st half of step
      return loss + compute(len / 2., kE + loss); // add loss from second half
   }

   double JoyLuoNieminenCSD::getBreakE() const
   {
      return breakE;
   }

   /**
   * Sets the
   *
   * @param breakE The breakE to set.
   */
   void JoyLuoNieminenCSD::setBreakE(double breakE)
   {
      this->breakE = breakE;
      setMaterial(mat, bh);
   }

   /**
   * @return - a string in the form "JoyLuoNieminenCSD(material,barrier
   *         height,breakE)", where material, barrier height, and breakE are
   *         the parameters either supplied to the constructor or, in the case
   *         of barrier height, possibly ascertained from the material
   *         property.
   */
   StringT JoyLuoNieminenCSD::toString()
   {
      return "JoyLuoNieminenCSD(" + mat.toString() + "," + std::to_string(bh) + "," + std::to_string(breakE) + ")";
   }
}
