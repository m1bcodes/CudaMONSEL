#include "gov\nist\microanalysis\NISTMonte\NullMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"

namespace NullMaterialScatterModel
{

   NullMaterialScatterModel::NullMaterialScatterModel() : minEforTracking(ToSI::eV(0.1))
   {
   }

   ElectronT* NullMaterialScatterModel::barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const
   {
      //minEforTracking = ToSI::eV(0.1);
      pe->setCurrentRegion(nextRegion);
      pe->setScatteringElement(NULL);
      return NULL;
   }

   double NullMaterialScatterModel::calculateEnergyLoss(double len, const ElectronT& pe) const
   {
      return 0.0;
   }

   const MaterialT& NullMaterialScatterModel::getMaterial() const
   {
      return Material::Default;
   }

   double NullMaterialScatterModel::randomMeanPathLength(ElectronT& pe)
   {
      return 1.0;
   }

   ElectronT* NullMaterialScatterModel::scatter(ElectronT& pe)
   {
      return NULL;
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