#ifndef _EXP_QM_BARRIER_SM_CUH_
#define _EXP_QM_BARRIER_SM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\BarrierScatterMechanism.cuh"

namespace ExpQMBarrierSM
{
   class ExpQMBarrierSM : public BarrierScatterMechanismT
   {
   public:
      __host__ __device__ ExpQMBarrierSM(const MaterialT* mat);
      ExpQMBarrierSM(const MaterialT* mat, double lambda);

      __host__ __device__ ElectronT* barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const override;
      StringT toString() const;

   private:
      __host__ __device__ double generalBarrierT(double rootPerpE, double rootDiff) const;

      const double u0; // Barrier height of this material to vacuum.
      const double lambda; // Barrier width
      const bool classical;
      const double lambdaFactor;
      const MaterialT* mat;
   };
}

#endif