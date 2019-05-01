#ifndef _BASIC_MATERIAL_MODEL_CUH_
#define _BASIC_MATERIAL_MODEL_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\IMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AlgorithmUser.cuh"

namespace BasicMaterialModel
{
   class BasicMaterialModel : public AlgorithmUserT, public IMaterialScatterModelT
   {
   public:
      BasicMaterialModel(const MaterialT& mat);

      const MaterialT& getMaterial() const;
      //double randomMeanPathLength(const ElectronT& pe) const;

   private:
      const MaterialT& mMaterial;
      double minEforTracking;
   };
}

#endif