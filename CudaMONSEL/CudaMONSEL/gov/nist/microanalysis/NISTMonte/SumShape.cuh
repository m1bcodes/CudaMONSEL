// file: gov\nist\microanalysis\NISTMonte\SumShape.cuh

#ifndef _SUM_SHAPE_CUH_
#define _SUM_SHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"
#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ITransform.cuh"

#include "Amphibian\vector.cuh"

namespace SumShape
{
   class SumShape : public ShapeT, public ITransformT//, TrajectoryVRML.IRender
   {
   public:
      SumShape(ShapeT* const shapes[], int len);
      SumShape(ShapeT* const a, ShapeT* const b);

      bool contains(const double pos[]) const override;
      double getFirstIntersection(const double pos0[], const double pos1[]) override;
      StringT toString() const override;
      bool isNormalShape() const override;

      void rotate(const double pivot[], double phi, double theta, double psi) override;
      void translate(const double distance[]) override;

      std::vector<ShapeT*> getShapes() const;

   private:
      std::vector<ShapeT*> mShapes;
   };
}

#endif