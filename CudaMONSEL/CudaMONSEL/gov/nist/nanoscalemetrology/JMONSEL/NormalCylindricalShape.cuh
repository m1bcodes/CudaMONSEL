#ifndef _NORMAL_CYLINDRICAL_SHAPE_CUH_
#define _NORMAL_CYLINDRICAL_SHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\CylindricalShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalShape.cuh"

namespace NormalCylindricalShape
{
   class NormalCylindricalShape : public CylindricalShapeT, public NormalShapeT
   {
   public:
      NormalCylindricalShape(const double end0[], const double end1[], double radius);

      bool contains(const double pos[]) const override;
      double getFirstIntersection(const double pos0[], const double pos1[]) override;
      StringT toString() const override;

      bool isNormalShape() const override;

      bool contains(const double pos0[], const double pos1[]) const override;
      VectorXd getFirstNormal(const double pos0[], const double pos1[]) override;
      VectorXd getPreviousNormal() const override;

      void rotate(const double pivot[], double phi, double theta, double psi) override;
      void translate(const double distance[]) override;

   private:
      double end0[3]; // Vector position of the 1st end cap center
      double axis[3]; // Vector from c to center of 2nd end cap
      double radius2; // Square of radius
      double normalizedaxis[3]; // cache axis normalized to unit length
      double mLen2; // Square of axis length
      double mLen; // axis length
      VectorXd nv; // Most recent normal vector
   };
}

#endif