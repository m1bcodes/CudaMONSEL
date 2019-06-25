#ifndef _BASIC_MATERIAL_MODEL_CUH_
#define _BASIC_MATERIAL_MODEL_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\IMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AlgorithmUser.cuh"

namespace BasicMaterialModel
{
   class BasicMaterialModel : public AlgorithmUserT, public IMaterialScatterModelT
   {
      void initializeDefaultStrategy() override;

   public:
      BasicMaterialModel(const MaterialT& mat);

      const MaterialT& getMaterial() const override;
      double randomMeanPathLength(ElectronT& pe) override;

      ElectronT* scatter(ElectronT& pe) override;

      ElectronT* barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const override;

      double calculateEnergyLoss(double len, const ElectronT& pe) const override;

      double getMinEforTracking() const override;
      __host__ __device__ void setMinEforTracking(double minEforTracking) override;

   private:
      const MaterialT& mMaterial;
      double minEforTracking;
   };
}

#endif