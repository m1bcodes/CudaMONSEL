// file : gov\nist\nanoscalemetrology\JMONSELutils\NULagrangeInterpolation.cuh

#ifndef _NU_LAGRNAGE_INTERPOLATION_CUH_
#define _NU_LAGRNAGE_INTERPOLATION_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace NULagrangeInterpolation
{
   extern VectorXd d1(double fsamp[], int fsamplen, double xsamp[], int xsamplen, int order, double x);
   extern VectorXd d2(MatrixXd f, MatrixXd xsamp, int order, double x[], int xlen);
   extern VectorXd d3(Matrix3DXd f, MatrixXd xsamp, int order, double x[], int xlen);
   extern VectorXd d4(Matrix4DXd f, MatrixXd xsamp, int order, double x[], int xlen);
}

#endif