#ifndef _MATH_2_CUH_
#define _MATH_2_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include <cstdlib> 
#include <ctime>

namespace Math2
{
   extern PositionVecT minus(PositionVecT a, PositionVecT b);
   extern PositionVecT divide(PositionVecT a, double b);
   extern PositionVecT plus(PositionVecT a, PositionVecT b);
   extern PositionVecT normalize(PositionVecT p);
   extern PositionVecT multiply(double a, PositionVecT b);
   extern double distance(PositionVecT p1, PositionVecT p2);
   extern double distance3d(const double p1[], const double p2[]);
   extern double dot(PositionVecT a, PositionVecT b);
   extern double magnitude3d(PositionVecT p);
   extern double sqr(double x);
   extern double magnitude(PositionVecT p);
   extern double random();

   extern double toRadians(double deg);

   extern double MINUS_Z_AXIS[];

   extern const double PI;
}

#endif