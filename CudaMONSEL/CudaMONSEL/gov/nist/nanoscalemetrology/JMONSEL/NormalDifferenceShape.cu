#include "gov\nist\nanoscalemetrology\JMONSEL\NormalDifferenceShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalComplementShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalIntersectionShape.cuh"

namespace NormalDifferenceShape
{
   __host__ __device__ NormalDifferenceShape::NormalDifferenceShape(NormalShapeT& shapeA, NormalShapeT& shapeB) : shapeA(shapeA), shapeB(shapeB)
   {
   }

   __host__ __device__ bool NormalDifferenceShape::contains(const double pos0[], const double pos1[]) const
   {
      return shapeA.contains(pos0, pos1) && !shapeB.contains(pos0, pos1);
   }

   __host__ __device__ const double* NormalDifferenceShape::getFirstNormal(const double pos0[], const double pos1[])
   {
      NormalComplementShapeT ncs(shapeB);
      NormalIntersectionShapeT nis(shapeA, ncs);
      memcpy(nv, nis.getFirstNormal(pos0, pos1), sizeof(nv[0]));
      return nv;
   }

   __host__ __device__ bool NormalDifferenceShape::contains(const double pos[]) const
   {
      return ((ShapeT&)shapeA).contains(pos) && !((ShapeT&)shapeB).contains(pos);
   }

   __host__ __device__ double NormalDifferenceShape::getFirstIntersection(const double pos0[], const double pos1[])
   {
      //memcpy(nv, getFirstNormal(pos0, pos1), sizeof(nv[0]));
      return getFirstNormal(pos0, pos1)[3];
   }

   __host__ __device__ void NormalDifferenceShape::rotate(const double pivot[], double phi, double theta, double psi)
   {
      //if (!(shapeA instanceof ITransform)) throw new EPQFatalException(shapeA.toString() + " does not support transformation.");
      //if (!(shapeB instanceof ITransform)) throw new EPQFatalException(shapeB.toString() + " does not support transformation.");
      ((ITransformT&)shapeA).rotate(pivot, phi, theta, psi);
      ((ITransformT&)shapeB).rotate(pivot, phi, theta, psi);
   }

   __host__ __device__ void NormalDifferenceShape::translate(const double distance[])
   {
      //if (!(shapeA instanceof ITransform)) throw new EPQFatalException(shapeA.toString() + " does not support transformation.");
      //if (!(shapeB instanceof ITransform)) throw new EPQFatalException(shapeB.toString() + " does not support transformation.");
      ((ITransformT&)shapeA).translate(distance);
      ((ITransformT&)shapeB).translate(distance);
   }

   __host__ __device__ const double* NormalDifferenceShape::getPreviousNormal() const
   {
      return nv;
   }

   __host__ __device__ StringT NormalDifferenceShape::toString() const
   {
      return "NormalDifferenceShape";
   }
}