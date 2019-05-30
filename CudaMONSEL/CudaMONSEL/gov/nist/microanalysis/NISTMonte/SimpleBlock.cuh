#ifndef _SIMPLE_BLOCK_CUH_
#define _SIMPLE_BLOCK_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"

namespace SimpleBlock
{
   class SimpleBlock : public ShapeT //,public TrajectoryVRML.IRender
   {
   public:
      SimpleBlock(const double corner0[], const double corner1[]);

      bool contains(const double pos[]) const;
      double getFirstIntersection(const double pos0[], const double pos1[]);
      StringT toString() const;

      const double* getCorner0() const;
      const double* getCorner1() const;

   private:
      double mCorner0[3];
      double mCorner1[3];
   };
}

#endif