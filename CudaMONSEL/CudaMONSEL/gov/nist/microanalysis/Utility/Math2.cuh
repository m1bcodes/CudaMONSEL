#ifndef _MATH_2_CUH_
#define _MATH_2_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include <cstdlib> 
#include <ctime>

namespace Math2
{
   extern VectorXd minus(const VectorXd& a, const VectorXd& b);
   extern VectorXd divide(const VectorXd& a, double b);
   extern VectorXd plus(const VectorXd& a, const VectorXd& b);
   extern VectorXd normalize(const VectorXd& p);
   extern VectorXd multiply(double a, const VectorXd& b);
   extern double distance(const VectorXd& p1, const VectorXd& p2);
   extern double distance3d(const VectorXd& p1, const VectorXd& p2);
   extern double dot(const VectorXd& a, const VectorXd& b);
   extern VectorXd cross(const VectorXd& a, const VectorXd& b);
   extern double sqr(double x);
   extern double magnitude(const VectorXd& p);
   extern VectorXd pointBetween(const VectorXd& a, const VectorXd& b, double f);
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