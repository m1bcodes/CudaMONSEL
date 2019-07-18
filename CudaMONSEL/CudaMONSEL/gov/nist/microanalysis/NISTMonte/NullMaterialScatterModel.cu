#include "gov\nist\microanalysis\NISTMonte\NullMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"

namespace NullMaterialScatterModel
{
   __host__ __device__ NullMaterialScatterModel::NullMaterialScatterModel() : minEforTracking(ToSI::eV(0.1))
   {
   }

   __host__ __device__ ElectronT* NullMaterialScatterModel::barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const
   {
      //minEforTracking = ToSI::eV(0.1);
      pe->setCurrentRegion(nextRegion);
      pe->setScatteringElement(nullptr);
      return nullptr;
   }

   __host__ __device__ double NullMaterialScatterModel::calculateEnergyLoss(double len, const ElectronT& pe) const
   {
      return 0.0;
   }

   __host__ __device__ const MaterialT& NullMaterialScatterModel::getMaterial() const
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      return *Material::d_Default;
#else
      return Material::Default;
#endif
   }

   __host__ __device__ double NullMaterialScatterModel::randomMeanPathLength(ElectronT& pe)
   {
      return 1.0;
   }

   __host__ __device__ ElectronT* NullMaterialScatterModel::scatter(ElectronT& pe)
   {
      return nullptr;
   }

   __host__ __device__ double NullMaterialScatterModel::getMinEforTracking() const
   {
      return minEforTracking;
   }

   __host__ __device__ void NullMaterialScatterModel::setMinEforTracking(double minEforTracking)
   {
      this->minEforTracking = minEforTracking;
   }
}