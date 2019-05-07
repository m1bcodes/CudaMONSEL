#ifndef _SLOW_DOWN_ALG_CUH_
#define _SLOW_DOWN_ALG_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace SlowingDownAlg
{
   class SlowingDownAlg
   {
   public:
      virtual void setMaterial(const SEmaterialT* mat) = 0;
      virtual double compute(double d, const ElectronT* pe) const = 0;
   };
}

#endif