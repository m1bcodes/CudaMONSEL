#ifndef _BASIC_MATERIAL_MODEL_CUH_
#define _BASIC_MATERIAL_MODEL_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\IMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AlgorithmUser.cuh"

namespace BasicMaterialModel
{
   class BasicMaterialModel : public AlgorithmUserT, public IMaterialScatterModelT
   {
      __host__ __device__ void initializeDefaultStrategy() override;

   public:
      BasicMaterialModel(const MaterialT& mat);

      __host__ __device__ const MaterialT& getMaterial() const override;
      __host__ __device__ double randomMeanPathLength(ElectronT& pe) override;

      __host__ __device__ ElectronT* scatter(ElectronT& pe) override;

      __host__ __device__ ElectronT* barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const override;

      __host__ __device__ double calculateEnergyLoss(double len, const ElectronT& pe) const override;

      __host__ __device__ double getMinEforTracking() const override;
      __host__ __device__ void setMinEforTracking(double minEforTracking) override;

   private:
      const MaterialT& mMaterial;
      double minEforTracking;
   };
}

#endif