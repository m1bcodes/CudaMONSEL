#ifndef _SCATTER_MECHANISM_CUH_
#define _SCATTER_MECHANISM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace ScatterMechanism
{
   typedef float data_type;

   class ScatterMechanism
   {
   public:
      __host__ __device__ virtual data_type scatterRate(const ElectronT* pe) = 0;
      __host__ __device__ virtual ElectronT* scatter(ElectronT* pe) = 0; // note: needs destructor
      __host__ __device__ virtual void setMaterial(const MaterialT* mat) = 0;
   };
}

#endif