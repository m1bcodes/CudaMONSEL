#ifndef _MATH_2_CUH_
#define _MATH_2_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include <cstdlib> 
#include <ctime>

namespace Math2
{
   PositionVecT minus(PositionVecT a, PositionVecT b);
   PositionVecT divide(PositionVecT a, double b);
   PositionVecT plus(PositionVecT a, PositionVecT b);
   PositionVecT normalize(PositionVecT p);
   PositionVecT multiply(double a, PositionVecT b);
   double distance(PositionVecT p1, PositionVecT p2);
   double distance3d(const double p1[], const double p2[]);
   double dot(PositionVecT a, PositionVecT b);
   double magnitude3d(PositionVecT p);
   double sqr(double x);
   double magnitude(PositionVecT p);
   double random();

   double toRadians(double deg);

   extern double MINUS_Z_AXIS[];

   const double PI = 3.14159265358979323846;
}

#endif