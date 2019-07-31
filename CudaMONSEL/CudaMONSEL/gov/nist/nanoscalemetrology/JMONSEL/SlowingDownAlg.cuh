#ifndef _SLOW_DOWN_ALG_CUH_
#define _SLOW_DOWN_ALG_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace SlowingDownAlg
{
   typedef float data_type;

   class SlowingDownAlg
   {
   public:
      __host__ __device__ virtual void setMaterial(const SEmaterialT* mat) = 0;
      __host__ __device__ virtual data_type compute(data_type d, const ElectronT* pe) const = 0;

      __host__ __device__ virtual StringT toString() const { return "class SlowingDownAlg"; }
   };
}

#endif