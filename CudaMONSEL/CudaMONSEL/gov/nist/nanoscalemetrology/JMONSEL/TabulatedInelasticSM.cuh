#ifndef _TABULATED_INELASTIC_SM_CUH_
#define _TABULATED_INELASTIC_SM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ScatterMechanism.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\SlowingDownAlg.cuh"

namespace TabulatedInelasticSM
{
   class TabulatedInelasticSM : public ScatterMechanismT
   {
   public:
      TabulatedInelasticSM();

      double scatterRate(const ElectronT* pe) override;
      ElectronT* scatter(ElectronT* pe) override;
      void setMaterial(const MaterialT* mat) override;

      StringT toString() const;

   private:
      const int methodSE;
      double energyOffset = 0.;
      double eFermi;
   };
}

#endif