#ifndef _U_LAGRANGE_INTERPOLATEION_CUH_
#define _U_LAGRANGE_INTERPOLATEION_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace ULagrangeInterpolation
{
   extern VectorXd d1(const double * const f, int len, double x0, double xinc, int order, double x);
   extern VectorXd d2(double const * const * const f, int flen0, int flen1, const double x0[], int x0len, const double xinc[], int xinclen, int order, const double x[], int xlen);
   extern VectorXd d3(double const * const * const * const f, int flen0, int flen1, int flen2, const double x0[], int x0len, const double xinc[], int xinclen, int order, const double x[], int xlen);
   extern VectorXd d4(double const * const * const * const * const f, int flen0, int flen1, int flen2, int flen3, double x0[], int x0len, double xinc[], int xinclen, int order, double x[], int xlen);

   extern VectorXd d1(const VectorXd& f, double x0, double xinc, int order, double x);
   extern VectorXd d2(const MatrixXd& f, const double x0[], int x0len, const double xinc[], int xinclen, int order, const double x[], int xlen);

   extern __host__ __device__ VectorXf d1(const VectorXf& f, const float x0, const float xinc, const int order, const float x);
   extern __host__ __device__ VectorXf d2(const MatrixXf& f, const float x0[], int x0len, const float xinc[], int xinclen, int order, const float x[], int xlen);
}

#endif