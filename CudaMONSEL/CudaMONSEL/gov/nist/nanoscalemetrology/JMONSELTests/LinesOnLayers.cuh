#ifndef _LINES_ON_LAYERS_CUH_
#define _LINES_ON_LAYERS_CUH_

#include <cuda_runtime.h>

namespace LinesOnLayers
{
   extern void run();
   __global__ extern void runCuda();
}

#endif