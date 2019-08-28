#ifndef _TRANSFORM_3D_CUH_
#define _TRANSFORM_3D_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace Transform3D
{
   extern __host__ __device__ MatrixXd rot(const double ph, const double th, const double ps);
   extern __host__ __device__ VectorXd rotate(const double vector[], const double phi, const double th, const double psi);
   extern __host__ __device__ VectorXd translate(const double point[], const double distance[], bool negate);
   extern __host__ __device__ VectorXd rotate(const double point[], const double pivot[], const double phi, const double theta, const double psi);

   extern __host__ __device__ void rotate3d(const double vector[], const double phi, const double th, const double psi, double[]);
   extern __host__ __device__ void translate3d(const double point[], const double distance[], bool negate, double[]);
   extern __host__ __device__ void rotate3d(const double point[], const double pivot[], const double phi, const double theta, const double psi, double[]);
}

#endif