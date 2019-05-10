#ifndef _I_TRANSFORM_CUH_
#define _I_TRANSFORM_CUH_

namespace ITransform
{
   class ITransform
   {
   public:
      virtual void rotate(const double pivot[], double phi, double theta, double psi) = 0;
      virtual void translate(const double distance[]) = 0;
   };
}

#endif