#ifndef _NORMAL_SHAPE_TRANSFORMER_CUH_
#define _NORMAL_SHAPE_TRANSFORMER_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace NormalShapeTransformer
{
   extern __host__ __device__ void translate(const double distance[], NormalShapeT& shape);
   extern __host__ __device__ void rotate(const double pivot[], double phi, double theta, double psi, NormalShapeT& shape);
};

#endif