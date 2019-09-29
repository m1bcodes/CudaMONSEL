#ifndef _NORMAL_COMPLEMENT_SHAPE_CUH_
#define _NORMAL_COMPLEMENT_SHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ITransform.cuh"
#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalShape.cuh"

namespace NormalComplementShape
{
   class NormalComplementShape : public ITransformT, public NormalShapeT
   {
   public:
      __host__ __device__ NormalComplementShape(NormalShapeT& shapeA);

      __host__ __device__ bool contains(const double pos[]) const override;
      __host__ __device__ double getFirstIntersection(const double pos0[], const double pos1[]) override;

      __host__ __device__ void rotate(const double pivot[], double phi, double theta, double psi) override;
      __host__ __device__ void translate(const double distance[]) override;

      __host__ __device__ const double* getFirstNormal(const double pos0[], const double pos1[]) override;
      __host__ __device__ bool contains(const double pos0[], const double pos1[]) const override;
      __host__ __device__ const double* getPreviousNormal() const override;

      __host__ __device__ StringT toString() const override;

   private:
      NormalShapeT& shapeA;
      double nv[4]; // Most recent normal vector
   };
}
#endif