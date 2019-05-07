#ifndef _NORMAL_SHAPE_CUH_
#define _NORMAL_SHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"

namespace NormalShape
{
   class NormalShape : public ShapeT
   {
   public:
      virtual bool contains(const double pos0[], const double pos1[]) const = 0;
      virtual PositionVecT getFirstNormal(double pos0[], double pos1[]) const = 0;
      virtual PositionVecT getPreviousNormal() const = 0;

      bool isNormalShape() const override {
         return true;
      };
   };
}

#endif