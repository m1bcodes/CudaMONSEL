#ifndef _MATH_2_CUH_
#define _MATH_2_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include <cstdlib> 
#include <ctime>

namespace Math2
{
   extern PositionVecT minus(const const PositionVecT a, const PositionVecT b);
   extern PositionVecT divide(const PositionVecT a, double b);
   extern PositionVecT plus(const PositionVecT a, const PositionVecT b);
   extern PositionVecT normalize(const PositionVecT p);
   extern PositionVecT multiply(double a, const PositionVecT b);
   extern double distance(const PositionVecT p1, const PositionVecT p2);
   extern double distance3d(const double p1[], const double p2[]);
   extern double dot(const PositionVecT a, const PositionVecT b);
   extern VectorXd cross(const VectorXd a, const VectorXd b);
   extern double sqr(double x);
   extern double magnitude(const PositionVecT p);
   extern double random();

   extern double toRadians(double deg);

   extern const double PI;
   extern const double ORIGIN_3D[];
   extern const double ONE[];
   extern const double X_AXIS[];
   extern const double Y_AXIS[];
   extern const double Z_AXIS[];
   extern const double MINUS_X_AXIS[];
   extern const double MINUS_Y_AXIS[];
   extern const double MINUS_Z_AXIS[];
}

#endif