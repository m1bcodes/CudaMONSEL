#ifndef _SELECTABLE_ELASTIC_SM_CUH_
#define _SELECTABLE_ELASTIC_SM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ScatterMechanism.cuh"

namespace SelectableElasticSM
{
   class SelectableElasticSM : public ScatterMechanismT
   {
      typedef std::vector<const RandomizedScatterT*> RandomizedScatterList;
   public:
      SelectableElasticSM(const MaterialT& mat, const RandomizedScatterFactoryT& rsf);
      SelectableElasticSM(const MaterialT& mat);
      //SelectableElasticSM(const SelectableElasticSM&); // clonable

      double scatterRate(const ElectronT* pe) override;
      ElectronT* scatter(ElectronT* pe) override;
      void setMaterial(const MaterialT* mat) override;

   private:
      void setCache(double kE);

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