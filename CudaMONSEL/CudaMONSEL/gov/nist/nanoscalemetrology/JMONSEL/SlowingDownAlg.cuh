#ifndef _SLOW_DOWN_ALG_CUH_
#define _SLOW_DOWN_ALG_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace SlowingDownAlg
{
   class SlowingDownAlg
   {
   public:
      __host__ __device__ virtual void setMaterial(const SEmaterialT* mat) = 0;
      __host__ __device__ virtual double compute(double d, const ElectronT* pe) const = 0;

      __host__ __device__ virtual StringT toString() const { return "class SlowingDownAlg"; }
   };
}

#endif