#ifndef _NORMAL_INTERSECTION_SHAPE_CUH_
#define _NORMAL_INTERSECTION_SHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ITransform.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalShape.cuh"

namespace NormalIntersectionShape
{
   class NormalIntersectionShape : public NormalShapeT, public ITransformT
   {
   public:
      __host__ __device__ NormalIntersectionShape(NormalShapeT& shapeA, NormalShapeT& shapeB);

      __host__ __device__ bool contains(const double pos[]) const override;
      __host__ __device__ double getFirstIntersection(const double pos0[], const double pos1[]) override;
      __host__ __device__ StringT toString() const override;

      __host__ __device__ void rotate(const double pivot[], double phi, double theta, double psi) override;
      __host__ __device__ void translate(const double distance[]) override;

      __host__ __device__ bool contains(const double pos0[], const double pos1[]) const  override;
      __host__ __device__ const double* getFirstNormal(const double pos0[], const double pos1[]) override;
      __host__ __device__ const double* getPreviousNormal() const override;

   private:
      NormalShapeT& shapeA;
      NormalShapeT& shapeB;
      double result[4];
   };

   //extern __host__ __device__ void rotateNormalShape(const double pivot[], const double phi, const double theta, const double psi, NormalShapeT& shape);
   //extern __host__ __device__ void translateNormalShape(const double distance[], NormalShapeT& shape);
}

#endif