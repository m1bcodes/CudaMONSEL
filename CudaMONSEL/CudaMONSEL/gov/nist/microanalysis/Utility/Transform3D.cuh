#ifndef _TRANSFORM_3D_CUH_
#define _TRANSFORM_3D_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace Transform3D
{
   extern MatrixXd rot(double ph, double th, double ps);
   extern VectorXd rotate(const double vector[], double phi, double th, double psi);
   extern VectorXd translate(const double point[], const double distance[], bool negate);
   extern VectorXd rotate(const double point[], const double pivot[], double phi, double theta, double psi);

   extern void rotate3d(const double vector[], double phi, double th, double psi, double[]);
   extern void translate3d(const double point[], const double distance[], bool negate, double[]);
   extern void rotate3d(const double point[], const double pivot[], double phi, double theta, double psi, double[]);
}

#endif