#ifndef _LINES_ON_LAYERS_0_CUH_
#define _LINES_ON_LAYERS_0_CUH_

#include <cuda_runtime.h>

namespace LinesOnLayers
{
   extern void run();
   extern void run1();

   extern void initCuda();
   extern __global__ void runCuda();

   class LinesOnLayers
   {
   public:
   private:
   };
}

#endif