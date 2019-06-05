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
      NormalIntersectionShape(NormalShapeT& shapeA, NormalShapeT& shapeB);

      bool contains(const double pos[]) const override;
      double getFirstIntersection(const double pos0[], const double pos1[]) override;
      StringT toString() const override;

      void rotate(const double pivot[], double phi, double theta, double psi) override;
      void translate(const double distance[]) override;

      bool contains(const double pos0[], const double pos1[]) const  override;
      const double* getFirstNormal(const double pos0[], const double pos1[]) override;
      const double* getPreviousNormal() const override;

   private:
      NormalShapeT& shapeA;
      NormalShapeT& shapeB;
      double result[4];
   };
}

#endif