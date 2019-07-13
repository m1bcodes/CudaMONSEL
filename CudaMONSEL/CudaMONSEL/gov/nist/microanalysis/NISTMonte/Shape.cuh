#ifndef _SHAPE_CUH_
#define _SHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace Shape
{
   class Shape
   {
   public:
      __host__ __device__ virtual bool contains(const double pos[]) const = 0;
      __host__ __device__ virtual double getFirstIntersection(const double pos0[], const double pos1[]) = 0;
      __host__ __device__ virtual StringT toString() const = 0;

      __host__ __device__ virtual bool isNormalShape() const
      {
         return false;
      }
   };
}

#endif