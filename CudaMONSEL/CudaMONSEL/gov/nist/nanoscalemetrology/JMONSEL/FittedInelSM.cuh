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
      __host__ __device__ FittedInelSM(const SEmaterialT& mat, double energySEgen, const SlowingDownAlgT& sdAlg);

      __host__ __device__ double scatterRate(const ElectronT* pe) override;
      __host__ __device__ ElectronT* scatter(ElectronT* pe) override;
      __host__ __device__ void setMaterial(const MaterialT* mat) override;

      StringT toString() const;

   private:
      const SlowingDownAlgT& sdAlg;
      const double energySEgen; // Average energy for SE generation
      double eFermi;
   };
}

#endif