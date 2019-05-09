#ifndef _FITTED_INEL_SM_CUH_
#define _FITTED_INEL_SM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ScatterMechanism.cuh"

namespace GanachaudMokraniPolaronTrapSM
{
   class GanachaudMokraniPolaronTrapSM : public ScatterMechanismT
   {
   public:
      GanachaudMokraniPolaronTrapSM(double prefactor, double extinction);

      double scatterRate(const ElectronT* pe) override;
      ElectronT* scatter(ElectronT* pe) override;
      void setMaterial(const MaterialT* mat) override;

      StringT toString() const;

   private:
      double prefactor;
      double extinction;
      double CUTOFF;
   };
}

#endif