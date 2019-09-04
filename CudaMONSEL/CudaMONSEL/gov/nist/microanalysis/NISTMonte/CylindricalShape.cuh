#ifndef _CYLINDICAL_SHAPE_CUH_
#define _CYLINDICAL_SHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"
#include "gov\nist\microanalysis\EPQlibrary\ITransform.cuh"

namespace CylindricalShape
{
   class CylindricalShape : public ShapeT, public ITransformT
   {
   public:
      __host__ __device__ CylindricalShape(const double end0[], const double end1[], double radius);
      __host__ __device__ CylindricalShape(const CylindricalShape&);

      __host__ __device__ bool contains(const double pos[]) const override;
      __host__ __device__ double getFirstIntersection(const double pos0[], const double pos1[]) override;
      __host__ __device__ StringT toString() const override;

      __host__ __device__ void rotate(const double pivot[], double phi, double theta, double psi) override;
      __host__ __device__ void translate(const double distance[]) override;

      __host__ __device__ const double* getEnd0() const;
      __host__ __device__ const double* getEnd1();

      __host__ __device__ double getRadius() const;
      double getLength() const;

   private:
      __host__ __device__ double closestPointOnAxis(const double p[]) const;
      __host__ __device__ double distanceSqr(const double p[], double u) const;

      double mEnd0[3]; // The position of the center of one end cap
      double mEnd1[3]; // redundant... = mEnd0 + mDelta
      double mDelta[3]; // The length and direction of the axis
      double mRadius2; // The sqr(radius) of the cylinder
      double mLen2; // Cache the length squared...
      double mDelta2;
   };
}

#endif