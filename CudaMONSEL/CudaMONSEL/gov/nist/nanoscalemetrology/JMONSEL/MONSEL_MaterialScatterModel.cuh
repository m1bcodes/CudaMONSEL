#ifndef _MONSEL_MATERIAL_SCATTER_MODEL_CUH_
#define _MONSEL_MATERIAL_SCATTER_MODEL_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\IMaterialScatterModel.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ZeroCSD.cuh"

#include "Amphibian\vector.cuh"

namespace MONSEL_MaterialScatterModel
{
   class MONSEL_MaterialScatterModel : public IMaterialScatterModelT
   {
      typedef amp::vector<ScatterMechanismT*> ScatterMechanismList;

   public:
      __host__ __device__ MONSEL_MaterialScatterModel(const SEmaterialT* mat, const BarrierScatterMechanismT* bsm, SlowingDownAlgT* sda);

      const BarrierScatterMechanismT* MONSEL_MaterialScatterModel::getBarrierSM() const;
      __host__ __device__ void setCSD(SlowingDownAlgT* csd);
      SlowingDownAlgT* getCSD();
      __host__ __device__ bool addScatterMechanism(ScatterMechanismT* mech);
      //bool removeScatterMechanism(const ScatterMechanismT* mech);
      //bool removeScatterMechanism(ScatterMechanismT* mech);

      __host__ __device__ const MaterialT& getMaterial() const override;
      __host__ __device__ double getMinEforTracking() const override;
      __host__ __device__ void setMinEforTracking(double minEforTracking) override;
      __host__ __device__ double randomMeanPathLength(ElectronT& pe) override;
      __host__ __device__ ElectronT* scatter(ElectronT& pe) override;
      __host__ __device__ ElectronT* barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const override;
      __host__ __device__ double calculateEnergyLoss(double len, const ElectronT& pe) const override;

   private:
      __host__ __device__ void setCache(const ElectronT* pe);

      const SEmaterialT* mat; // The material for which this is the scatter model

      double minEforTracking;

      // Following variables are associated with scattering
      double cached_eK = -1.; // Initialize to an impossible eK.

      VectorXd cached_cumulativeScatterRate;

      double totalScatterRate;

      ScatterMechanismList scatterArray;

      int nscattermech;

      /*
      * The continuous slowing down algorithm to use for this material
      */
      SlowingDownAlgT* csd;

      /*
      * The barrier scattering mechanism to use for this material
      */
      const BarrierScatterMechanismT* barrierSM;

      /*
      * A set of scatter mechanisms at work in this material. Use of a set type
      * instead of a list type is meant to discourage duplicate entries. (Each
      * scatter mechanism should appear only once. Use of a set collection
      * prevents the same object from being added multiple times. Note that
      * duplicates are still possible inasmuch as we must rely upon the user not
      * to add distinct objects that are nevertheless instances of the same
      * underlying mechanism.) LinkedHashSet is used rather than HashSet or
      * TreeSet because it guarantees iteration order to be the same as order of
      * insertion and because of better performance in iterations.
      */

      /*
      * TODO: Maintaining a separate LinkedHashSet and scatterArray is cumbersome.
      * I should probably convert this to an ArrayList (so I can use the
      * .get(index) method) and get rid of the scatterArray.
      */
   };
}

#endif