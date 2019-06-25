#ifndef _NULL_MATERIAL_SCATTER_MODEL_CUH_
#define _NULL_MATERIAL_SCATTER_MODEL_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\IMaterialScatterModel.cuh"

namespace NullMaterialScatterModel
{
   class NullMaterialScatterModel : public IMaterialScatterModelT
   {
   public:
      NullMaterialScatterModel();

      ElectronT* barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const override;

      double calculateEnergyLoss(double len, const ElectronT& pe) const override;

      const MaterialT& getMaterial() const override;

      double randomMeanPathLength(ElectronT& pe) override;

      ElectronT* scatter(ElectronT& pe) override;

      double getMinEforTracking() const override;

      __host__ __device__ void setMinEforTracking(double minEforTracking) override;

   private:
      double minEforTracking;
   };
}
#endif