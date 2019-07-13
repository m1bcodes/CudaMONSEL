#ifndef _NORMAL_SHAPE_CUH_
#define _NORMAL_SHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"

namespace NormalShape
{
   class NormalShape : public ShapeT
   {
   public:
      __host__ __device__ virtual bool contains(const double pos0[], const double pos1[]) const = 0;
      __host__ __device__ virtual const double* getFirstNormal(const double pos0[], const double pos1[]) = 0;
      __host__ __device__ virtual const double* getPreviousNormal() const = 0;

      __host__ __device__ bool isNormalShape() const override {
         return true;
      };
   };
}

#endif