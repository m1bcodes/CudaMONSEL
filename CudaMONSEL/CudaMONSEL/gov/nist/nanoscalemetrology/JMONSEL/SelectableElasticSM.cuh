#ifndef _SELECTABLE_ELASTIC_SM_CUH_
#define _SELECTABLE_ELASTIC_SM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ScatterMechanism.cuh"

#include "Amphibian\vector.cuh"

namespace SelectableElasticSM
{
   typedef ::ScatterMechanism::data_type data_type;

   class SelectableElasticSM : public ScatterMechanismT
   {
      typedef amp::vector<const RandomizedScatterT*> RandomizedScatterList;

   public:
      __host__ __device__ SelectableElasticSM(const MaterialT& mat, const RandomizedScatterFactoryT& rsf);
      __host__ __device__ SelectableElasticSM(const MaterialT& mat);
      //SelectableElasticSM(const SelectableElasticSM&); // clonable

      __host__ __device__ data_type scatterRate(const ElectronT* pe) override;
      __host__ __device__ ElectronT* scatter(ElectronT* pe) override;
      __host__ __device__ void setMaterial(const MaterialT* mat) override;

   private:
      __host__ __device__ void setCache(data_type kE);

      RandomizedScatterList rse;
      // Set scatter class default to NISTMottScatteringAngle
      const RandomizedScatterFactoryT& rsf;

      VectorXd scalefactor; // weight fraction/atomic weight
      /* We use cross sections divided by atomic weight */
      VectorXd cumulativeScaledCrossSection;

      data_type totalScaledCrossSection;
      int nce; // # constituent elements

      data_type densityNa; // Avagadro's # * density for this material

      data_type cached_kE; // Initialize to impossible value
   };
}

#endif