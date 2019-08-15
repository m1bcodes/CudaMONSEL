#ifndef _SPHERE_CUH_
#define _SPHERE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ITransform.cuh"

namespace Sphere
{
   class Sphere : public ShapeT, public ITransformT//, TrajectoryVRML.IRender
   {
   public:
      __host__ __device__ Sphere(const double center[], double radius);

      __host__ __device__ bool contains(const double pos[]) const override;
      double getRadius() const;
      __host__ __device__ double getFirstIntersection(const double pos0[], const double pos1[]) override;

      void getInitialPoint(int res[]) const;
      void getPointAt(double phi, double theta, double frac, double res[]) const;

      __host__ __device__ void rotate(const double pivot[], double phi, double theta, double psi) override;
      __host__ __device__ void translate(const double distance[]) override;

      const double* getCenter() const;
      __host__ __device__ StringT toString() const override;

   private:
      const double mRadius; // meters
      double mCenter[3]; // x,y & z in meters
   };
}

#endif