#ifndef _FITTED_INEL_SM_CUH_
#define _FITTED_INEL_SM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ScatterMechanism.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\SlowingDownAlg.cuh"

namespace FittedInelSM
{
   class FittedInelSM : public ScatterMechanismT
   {
   public:
      FittedInelSM(const SEmaterialT& mat, double energySEgen, const SlowingDownAlgT& sdAlg);

      double scatterRate(const ElectronT* pe) override;
      ElectronT* scatter(ElectronT* pe) override;
      void setMaterial(const MaterialT* mat) override;

      StringT toString() const;

   private:
      const SlowingDownAlgT& sdAlg;
      const double energySEgen; // Average energy for SE generation
      double eFermi;
   };
}

#endif