#include "gov\nist\nanoscalemetrology\JMONSEL\NormalComplementShape.cuh"

#include "gov\nist\nanoscalemetrology\JMONSEL\NormalShapeTransformer.cuh"

namespace NormalComplementShape
{
   __host__ __device__ NormalComplementShape::NormalComplementShape(NormalShapeT& shapeA) : shapeA(shapeA)
   {
   }

   __host__ __device__ bool NormalComplementShape::contains(const double pos[]) const
   {
      return !((ShapeT&)shapeA).contains(pos);
   }

   __host__ __device__ double NormalComplementShape::getFirstIntersection(const double pos0[], const double pos1[])
   {
      //memcpy(nv, getFirstNormal(pos0, pos1), sizeof(nv[0]) * 4);
      return getFirstNormal(pos0, pos1)[3];
   }

   __host__ __device__ void NormalComplementShape::rotate(const double pivot[], double phi, double theta, double psi)
   {
      //if (!(shapeA instanceof ITransform)) throw new EPQFatalException(shapeA.toString() + " does not support transformation.");
      //((ITransformT&)shapeA).rotate(pivot, phi, theta, psi);

      NormalShapeTransformer::rotate(pivot, phi, theta, psi, shapeA);
   }

   __host__ __device__ void NormalComplementShape::translate(const double distance[])
   {
      //if (!(shapeA instanceof ITransform)) throw new EPQFatalException(shapeA.toString() + " does not support transformation.");
      //((ITransformT&)shapeA).translate(distance);

      NormalShapeTransformer::translate(distance, shapeA);
   }

   __host__ __device__ const double* NormalComplementShape::getFirstNormal(const double pos0[], const double pos1[])
   {
      memcpy(nv, shapeA.getFirstNormal(pos0, pos1), sizeof(nv[0]) * 4); // Get normal vector of
      // shapeA
      nv[0] *= -1.; // Reverse direction of normal vector
      nv[1] *= -1.;
      nv[2] *= -1.;
      return nv;
   }

   __host__ __device__ bool NormalComplementShape::contains(const double pos0[], const double pos1[]) const
   {
      return !shapeA.contains(pos0, pos1);
   }

   __host__ __device__ const double* NormalComplementShape::getPreviousNormal() const
   {
      return nv;
   }

   __host__ __device__ StringT NormalComplementShape::toString() const
   {
      return "NormalComplementShape";
   }
}