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
      CylindricalShape(const double end0[], const double end1[], double radius);
      CylindricalShape(const CylindricalShape&);

      bool contains(const double pos[]) const override;
      double getFirstIntersection(const double pos0[], const double pos1[]) override;
      StringT toString() const override;

      void rotate(const double pivot[], double phi, double theta, double psi) override;
      void translate(const double distance[]) override;

      VectorXd getEnd0() const;
      VectorXd getEnd1() const;

      double getRadius() const;
      double getLength() const;

   private:
      double closestPointOnAxis(const double p[]) const;
      double distanceSqr(const double p[], double u) const;

      VectorXd mEnd0; // The position of the center of one end cap
      VectorXd mDelta; // The length and direction of the axis
      double mRadius2; // The sqr(radius) of the cylinder
      double mLen2; // Cache the length squared...
      double mDelta2;
   };
}

#endif