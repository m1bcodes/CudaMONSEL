#ifndef _SHAPE_CUH_
#define _SHAPE_CUH_

namespace Shape
{
   class Shape
   {
   public:
      virtual bool contains(const double pos[]) const = 0;
      virtual double getFirstIntersection(const double pos0[], const double pos1[]) = 0;
      virtual char const * toString() const = 0;
      virtual bool isNormalShape() const { return false; }
   };
}

#endif