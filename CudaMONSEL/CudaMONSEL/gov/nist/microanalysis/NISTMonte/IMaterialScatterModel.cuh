#ifndef _I_MATERIAL_SCATTER_MODEL_CUH_
#define _I_MATERIAL_SCATTER_MODEL_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

#include <cuda_runtime.h>

namespace IMaterialScatterModel
{
   class IMaterialScatterModel
   {
   public:
      __host__ __device__ virtual const MaterialT& getMaterial() const = 0;
      __host__ __device__ virtual double getMinEforTracking() const = 0;
      __host__ __device__ virtual void setMinEforTracking(double minEforTracking) = 0;
      __host__ __device__ virtual double randomMeanPathLength(ElectronT& pe) = 0;
      __host__ __device__ virtual ElectronT* scatter(ElectronT& pe) = 0;
      __host__ __device__ virtual ElectronT* barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const = 0;
      __host__ __device__ virtual double calculateEnergyLoss(double len, const ElectronT& pe) const = 0;
   };
}

#endif