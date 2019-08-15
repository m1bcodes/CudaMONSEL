#ifndef _TRANSFORM_3D_CUH_
#define _TRANSFORM_3D_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace Transform3D
{
   extern __host__ __device__ MatrixXd rot(double ph, double th, double ps);
   extern __host__ __device__ VectorXd rotate(const double vector[], double phi, double th, double psi);
   extern __host__ __device__ VectorXd translate(const double point[], const double distance[], bool negate);
   extern __host__ __device__ VectorXd rotate(const double point[], const double pivot[], double phi, double theta, double psi);

   extern __host__ __device__ void rotate3d(const double vector[], double phi, double th, double psi, double[]);
   extern __host__ __device__ void translate3d(const double point[], const double distance[], bool negate, double[]);
   extern __host__ __device__ void rotate3d(const double point[], const double pivot[], double phi, double theta, double psi, double[]);
}

#endif