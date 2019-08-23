#ifndef _LINES_ON_LAYERS_0_CUH_
#define _LINES_ON_LAYERS_0_CUH_

#include <cuda_runtime.h>

namespace LinesOnLayers
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   extern __device__ unsigned int ysize;
   extern __device__ unsigned int xsize;
#else
   extern unsigned int ysize;
   extern unsigned int xsize;
#endif
   extern void loadNUTable();
   extern void transferNUTableToCuda();

   extern __global__ void initCuda();
   extern __global__ void runCuda(float* result);

   extern __host__ __device__ void runSinglePixel(const unsigned int r, const unsigned int c, float* result);
   extern __host__ __device__ void initRange();

   extern void runSinglePixelThread(int id, const unsigned int r, const unsigned int c, float* result);

   extern void testLineProjection();
}

#endif