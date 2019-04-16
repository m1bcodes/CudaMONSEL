#ifndef _I_TRANSFORM_CUH_
#define _I_TRANSFORM_CUH_

class ITransform
{
   virtual void rotate(double pivot[], double phi, double theta, double psi) = 0;
   virtual void translate(double distance[]) = 0;
};

#endif