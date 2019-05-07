#ifndef _I_MATERIAL_SCATTER_MODEL_CUH_
#define _I_MATERIAL_SCATTER_MODEL_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace IMaterialScatterModel
{
   class IMaterialScatterModel
   {
   public:
      virtual const MaterialT& getMaterial() const = 0;
      virtual double getMinEforTracking() const = 0;
      virtual void setMinEforTracking(double minEforTracking) = 0;
      virtual double randomMeanPathLength(ElectronT& pe) = 0;
      virtual ElectronT* scatter(ElectronT& pe) = 0;
      virtual ElectronT* barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const = 0;
      virtual double calculateEnergyLoss(double len, const ElectronT& pe) const = 0;
   };
}

#endif