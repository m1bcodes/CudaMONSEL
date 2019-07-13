// file : gov\nist\nanoscalemetrology\JMONSELutils\NULagrangeInterpolation.cuh

#ifndef _NU_LAGRNAGE_INTERPOLATION_CUH_
#define _NU_LAGRNAGE_INTERPOLATION_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace NULagrangeInterpolation
{
   extern __host__ __device__ VectorXd d1(const double fsamp[], int fsamplen, const double xsamp[], int xsamplen, int order, double x);
   extern __host__ __device__ VectorXd d2(const MatrixXd& f, const MatrixXd& xsamp, int order, const double x[], int xlen);
   extern __host__ __device__ VectorXd d3(const Matrix3DXd& f, const MatrixXd& xsamp, int order, const double x[], int xlen);
   extern __host__ __device__ VectorXd d4(const Matrix4DXd& f, const MatrixXd& xsamp, int order, const double x[], int xlen);
}

#endif