#ifndef _SCATTER_MECHANISM_CUH_
#define _SCATTER_MECHANISM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace ScatterMechanism
{
   class ScatterMechanism
   {
   public:
      virtual double scatterRate(const ElectronT* pe) = 0;
      virtual ElectronT* scatter(ElectronT* pe) = 0; // WARNING: needs destructor
      virtual void setMaterial(const MaterialT* mat) = 0;
   };
}

#endif