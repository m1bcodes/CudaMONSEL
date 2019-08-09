#ifndef _LINES_ON_LAYERS_0_CUH_
#define _LINES_ON_LAYERS_0_CUH_

#include <cuda_runtime.h>

namespace LinesOnLayers
{
   extern void loadNUTable();
   extern void transferNUTableToCuda();

   //extern void initCuda();
   //extern void runCuda();

   extern __global__ void initCuda();
   extern __global__ void runCuda(float* result);

   extern __host__ __device__ void runSinglePixel(const unsigned int r, const unsigned int c, float* result);
   extern __host__ __device__ void initRange();
   extern __host__ __device__ void destroyRange();

   extern void runSinglePixelThread(int id, const unsigned int r, const unsigned int c, float* result);

   extern void testLineProjection();
}

#endif