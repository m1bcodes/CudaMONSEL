#ifndef _I_TRANSFORM_CUH_
#define _I_TRANSFORM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace ITransform
{
   class ITransform
   {
   public:
      virtual __host__ __device__ void rotate(const double pivot[], double phi, double theta, double psi) = 0;
      virtual __host__ __device__ void translate(const double distance[]) = 0;
   };
}

#endif