#ifndef _LINES_ON_LAYERS_0_CUH_
#define _LINES_ON_LAYERS_0_CUH_

#include <cuda_runtime.h>

namespace LinesOnLayers
{
   extern void copyDataToCuda();

   //extern __global__ void initCuda();
   //extern __global__ void runCuda();

   extern __global__ void initCuda();
   extern __global__ void runCuda();
}

#endif