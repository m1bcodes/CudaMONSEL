#ifndef _NULL_MATERIAL_SCATTER_MODEL_CUH_
#define _NULL_MATERIAL_SCATTER_MODEL_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\IMaterialScatterModel.cuh"

namespace NullMaterialScatterModel
{
   class NullMaterialScatterModel : public IMaterialScatterModelT
   {
   public:
      __host__ __device__ NullMaterialScatterModel();

      __host__ __device__ ElectronT* barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const override;

      __host__ __device__ double calculateEnergyLoss(double len, const ElectronT& pe) const override;

      __host__ __device__ const MaterialT& getMaterial() const override;

      __host__ __device__ double randomMeanPathLength(ElectronT& pe) override;

      __host__ __device__ ElectronT* scatter(ElectronT& pe) override;

      __host__ __device__ double getMinEforTracking() const override;

      __host__ __device__ void setMinEforTracking(double minEforTracking) override;

   private:
      double minEforTracking;
   };
}
#endif