#include "gov\nist\nanoscalemetrology\JMONSEL\SelectableElasticSM.cuh"
#include "gov\nist\microanalysis\EPQLibrary\NISTMottScatteringAngle.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
#include "gov\nist\microanalysis\EPQLibrary\PhysicalConstants.cuh"

#include "Amphibian\random.cuh"

namespace SelectableElasticSM
{
   __host__ __device__ SelectableElasticSM::SelectableElasticSM(const MaterialT& mat, const RandomizedScatterFactoryT& rsf) : rsf(rsf), cached_kE(-1.)
   {
      setMaterial(&mat);
   }

   __host__ __device__ SelectableElasticSM::SelectableElasticSM(const MaterialT& mat) :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      rsf(*NISTMottScatteringAngle::d_Factory), cached_kE(-1.)
#else
      rsf(NISTMottScatteringAngle::Factory), cached_kE(-1.f)
#endif
   {
      setMaterial(&mat);
   }

   //SelectableElasticSM::SelectableElasticSM(const SelectableElasticSM& sm) : rsf(sm.rsf), 
   //{
   //   setMaterial(sm.mat);
   //}

   __host__ __device__ void SelectableElasticSM::setCache(data_type kE)
   {
      /*
      * Algorithm: 1. Get scaled cross section (cross section times weight
      * fraction divided by atomic weight) for each element in this material 2.
      * From this, determine the total scaled cross section. 3. Cache these for
      * later use.
      */
      totalScaledCrossSection = 0.;
      for (int i = 0; i < nce; i++) {
         totalScaledCrossSection += rse[i]->totalCrossSection(kE) * scalefactor[i];
         cumulativeScaledCrossSection[i] = totalScaledCrossSection;
      }
      // Remember kinetic energy for which the cache was created
      cached_kE = kE;
   }

   __host__ __device__ data_type SelectableElasticSM::scatterRate(const ElectronT* pe)
   {
      setCache(pe->getEnergy()); // computes totalScaledCrossSection for this eK
      return totalScaledCrossSection * densityNa;
   }

   __host__ __device__ ElectronT* SelectableElasticSM::scatter(ElectronT* pe)
   {
      const data_type kE = pe->getPreviousEnergy();
      if (kE != cached_kE)
         setCache(kE);
      // Decide which element we scatter from
      const data_type r = Random::random() * totalScaledCrossSection;
      int index = 0; // Index is first index

      // Increment index and mechanism until cumulative scatter rate exceeds r
      while (cumulativeScaledCrossSection[index] < r)
         index++;

      const data_type alpha = rse[index]->randomScatteringAngle(kE);
      if (alpha != alpha) printf("%s\n", rse[index]->toString().c_str());
      const data_type beta = 2.f * Math2::PI * Random::random();
      pe->updateDirection(alpha, beta);
      pe->setScatteringElement(&(rse[index]->getElement()));
      return nullptr; // This mechanism is elastic. No SE.
   }

   __host__ __device__ void SelectableElasticSM::setMaterial(const MaterialT* mat)
   {
      nce = mat->getElementCount();
      densityNa = mat->getDensity() * PhysicalConstants::AvagadroNumber;
      if (nce > 0) {
         // Element[] elements = (Element[]) mat.getElementSet().toArray();
         const Element::UnorderedSetT& elements = mat->getElementSet();
         rse.resize(nce);
         scalefactor.resize(nce);
         cumulativeScaledCrossSection.resize(nce);

         int i = 0;
         for (auto &elm : elements) {
            rse[i] = &rsf.get(*elm);
            // The factor of 1000 in the next line is to convert atomic
            // weight in g/mole to kg/mole.
            scalefactor[i] = (1000.f * mat->weightFraction(*elm, true)) / elm->getAtomicWeight();
            i++;
         }
      }
   }
}