#ifndef _SELECTABLE_ELASTIC_SM_CUH_
#define _SELECTABLE_ELASTIC_SM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ScatterMechanism.cuh"

#include "Amphibian\vector.cuh"

namespace SelectableElasticSM
{
   class SelectableElasticSM : public ScatterMechanismT
   {
      typedef amp::vector<const RandomizedScatterT*> RandomizedScatterList;

   public:
      __host__ __device__ SelectableElasticSM(const MaterialT& mat, const RandomizedScatterFactoryT& rsf);
      __host__ __device__ SelectableElasticSM(const MaterialT& mat);
      //SelectableElasticSM(const SelectableElasticSM&); // clonable

      __host__ __device__ double scatterRate(const ElectronT* pe) override;
      __host__ __device__ ElectronT* scatter(ElectronT* pe) override;
      __host__ __device__ void setMaterial(const MaterialT* mat) override;

   private:
      __host__ __device__ void setCache(double kE);

      RandomizedScatterList rse;
      // Set scatter class default to NISTMottScatteringAngle
      const RandomizedScatterFactoryT& rsf;

      VectorXd scalefactor; // weight fraction/atomic weight
      /* We use cross sections divided by atomic weight */
      VectorXd cumulativeScaledCrossSection;

      double totalScaledCrossSection;
      int nce; // # constituent elements

      double densityNa; // Avagadro's # * density for this material

      double cached_kE; // Initialize to impossible value
   };
}

#endif