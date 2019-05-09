#ifndef _JOY_LUO_NIEMINEN_CSD_CUH_
#define _JOY_LUO_NIEMINEN_CSD_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\SlowingDownAlg.cuh"

namespace JoyLuoNieminenCSD
{
   class JoyLuoNieminenCSD : public SlowingDownAlgT
   {
   public:
      JoyLuoNieminenCSD(SEmaterialT& mat, double breakE);
      JoyLuoNieminenCSD(MaterialT& mat, double bh, double breakE);

      void setMaterial(const MaterialT& mat, double bh);

      void init();

      void setMaterial(const SEmaterialT* mat) override;
      double compute(double d, const ElectronT* pe) const override;

      double compute(double len, double kE) const;

      double getBreakE() const;
      void setBreakE(double breakE);

      StringT toString();

   private:
      MaterialT& mat;

      double bh; // barrier height

      int nce; // # constituent elements

      /* Combinations of variables that we need */
      VectorXd recipJ; // 1.166/J where J=ionization energy of const.
      // elem.

      VectorXd coef; // 78500 rho c[i]Z[i]/A[i] = leading coefficient in

      VectorXd beta; // 1=1.166(wf+1eV)/J[i]

      double bhplus1eV; // BH + 1eV

      double minEforTracking;

      double breakE; // Dividing line between Joy/Luo and Nieminen
      double gamma; // The proportionality constant in Nieminen's formula
   };
}

#endif