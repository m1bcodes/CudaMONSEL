#include "gov\nist\microanalysis\NISTMonte\NullMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"

namespace NullMaterialScatterModel
{
   ElectronT& NullMaterialScatterModel::barrierScatter(ElectronT& pe, const RegionBaseT& nextRegion)
   {
      minEforTracking = ToSI::eV(0.1);
      pe.setCurrentRegion(nextRegion);
      pe.setScatteringElement(NULL);
      return Electron::Default;
   }

   double NullMaterialScatterModel::calculateEnergyLoss(double len, ElectronT& pe) const
   {
      return 0.0;
   }

   const MaterialT& NullMaterialScatterModel::getMaterial() const
   {
      return Material::Default;
   }

   double NullMaterialScatterModel::randomMeanPathLength(ElectronT& pe) const
   {
      return 1.0;
   }

   ElectronT& NullMaterialScatterModel::scatter(ElectronT& pe)
   {
      return Electron::Default;
   }

   double NullMaterialScatterModel::getMinEforTracking() const
   {
      return minEforTracking;
   }

   void NullMaterialScatterModel::setMinEforTracking(double minEforTracking)
   {
      this->minEforTracking = minEforTracking;
   }
}