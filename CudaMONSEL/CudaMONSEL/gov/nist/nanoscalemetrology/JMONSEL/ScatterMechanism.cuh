#ifndef _SCATTER_MECHANISM_CUH_
#define _SCATTER_MECHANISM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace ScatterMechanism
{
   class ScatterMechanism
   {
   public:
      __host__ __device__ virtual double scatterRate(const ElectronT* pe) = 0;
      __host__ __device__ virtual ElectronT* scatter(ElectronT* pe) = 0; // note: needs destructor
      __host__ __device__ virtual void setMaterial(const MaterialT* mat) = 0;
   };
}

#endif