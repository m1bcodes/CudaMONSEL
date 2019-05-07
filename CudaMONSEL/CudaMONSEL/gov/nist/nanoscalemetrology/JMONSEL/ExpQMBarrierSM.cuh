#ifndef _EXP_QM_BARRIER_SM_CUH_
#define _EXP_QM_BARRIER_SM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\BarrierScatterMechanism.cuh"

namespace ExpQMBarrierSM
{
   class ExpQMBarrierSM : public BarrierScatterMechanismT
   {
   public:
      ExpQMBarrierSM(const MaterialT* mat);
      ExpQMBarrierSM(const MaterialT* mat, double lambda);

      ElectronT* barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const override;
      StringT toString() const;

   private:
      double generalBarrierT(double rootPerpE, double rootDiff) const;

      double u0; // Barrier height of this material to vacuum.
      double lambda; // Barrier width
      bool classical;
      double lambdaFactor;
      const MaterialT* mat;
   };
}

#endif