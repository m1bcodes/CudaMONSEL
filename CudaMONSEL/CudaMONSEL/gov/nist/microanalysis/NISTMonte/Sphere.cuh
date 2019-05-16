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
      Sphere(const double center[], double radius);

      bool contains(const double pos[]) const override;
      double getRadius() const;
      double getFirstIntersection(const double pos0[], const double pos1[]) override;

      VectorXd getInitialPoint() const;
      VectorXd getPointAt(double phi, double theta, double frac) const;

      void rotate(const double pivot[], double phi, double theta, double psi) override;
      void translate(const double distance[]) override;

      VectorXd getCenter() const;
      StringT toString() const override;

   private:
      const double mRadius; // meters
      VectorXd mCenter; // x,y & z in meters
   };
}

#endif