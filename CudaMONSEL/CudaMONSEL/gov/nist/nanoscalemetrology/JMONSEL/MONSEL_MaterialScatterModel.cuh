#ifndef _MONSEL_MATERIAL_SCATTER_MODEL_CUH_
#define _MONSEL_MATERIAL_SCATTER_MODEL_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\IMaterialScatterModel.cuh"

namespace MONSEL_MaterialScatterModel
{
   class MONSEL_MaterialScatterModel : public IMaterialScatterModelT
   {
      typedef std::vector<ScatterMechanismT*> ScatterMechanismList;

   public:
      MONSEL_MaterialScatterModel(const SEmaterialT* mat, const BarrierScatterMechanismT* bsm);

      const BarrierScatterMechanismT* MONSEL_MaterialScatterModel::getBarrierSM() const;
      void setCSD(SlowingDownAlgT* csd);
      SlowingDownAlgT* getCSD();
      bool addScatterMechanism(ScatterMechanismT* mech);
      bool removeScatterMechanism(const ScatterMechanismT* mech);

      const MaterialT& getMaterial() const override;
      double getMinEforTracking() const override;
      void setMinEforTracking(double minEforTracking) override;
      double randomMeanPathLength(ElectronT& pe) override;
      ElectronT* scatter(ElectronT& pe) override;
      ElectronT* barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const override;
      double calculateEnergyLoss(double len, const ElectronT& pe) const override;

   private:
      void setCache(const ElectronT* pe);

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