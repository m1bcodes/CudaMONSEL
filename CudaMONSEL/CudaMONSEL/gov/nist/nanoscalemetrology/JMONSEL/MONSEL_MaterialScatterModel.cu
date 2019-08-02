#include "gov\nist\nanoscalemetrology\JMONSEL\MONSEL_MaterialScatterModel.cuh"

#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

#include "gov\nist\nanoscalemetrology\JMONSEL\SEmaterial.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ScatterMechanism.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\BarrierScatterMechanism.cuh"

#include "Amphibian\random.cuh"

namespace MONSEL_MaterialScatterModel
{
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//   __constant__ static ZeroCSDT sZeroCSD;
//#else
//   static ZeroCSDT sZeroCSD;
//#endif

   __host__ __device__ MONSEL_MaterialScatterModel::MONSEL_MaterialScatterModel(const SEmaterialT* mat, const BarrierScatterMechanismT* bsm, SlowingDownAlgT* sda) :
      mat(mat), barrierSM(bsm), csd(sda), minEforTracking(mat->getElementCount() == 0 ? -INFINITY : (-mat->getEnergyCBbottom() < 0.f ? 0.f : -mat->getEnergyCBbottom()))
   {
      // TODO: take note of difference wrt original JMONSEL code
      //this->mat = mat.clone();
      //barrierSM = new ExpQMBarrierSM(this.mat);

      /*
      * By default the minimum energy for tracking is negative infinity in
      * vacuum and the bottom of the conduction band in non-vacuum. This means
      * we never drop an electron from the simulation when it is in vacuum, but
      * rather wait for it to either enter a material or be detected.
      */
      //if (mat->getElementCount() == 0)
      //   minEforTracking = -INFINITY;
      //else {
         /*
         * Default to barrier height or 0. whichever is larger. (Barrier can be
         * negative in some negative electronic affinity materials, but in this
         * case we assume a negative kinetic energy, which means energy less
         * than the conduction band bottom, refers to a trapped state.
         */
      //   minEforTracking = -mat->getEnergyCBbottom();
      //   if (minEforTracking < 0.f) minEforTracking = 0.f;
      //}
   }

   __host__ __device__ const MaterialT& MONSEL_MaterialScatterModel::getMaterial() const
   {
      return *mat;
   }

   __host__ __device__ void MONSEL_MaterialScatterModel::setCache(const ElectronT* pe)
   {
      // Remember kinetic energy for which the cache was created
      cached_eK = pe->getEnergy();
      /*
      * Algorithm: 1. Get scatter rate (1/(mean free path) for each mechanism
      * active in this material 2. From this, determine the total scatter rate.
      * 3. Cache the scatter rates for later use.
      */
      totalScatterRate = 0.;
      if (scatterArray.empty())
         return;
      int index = 0;
      for (auto &mech : scatterArray) {
         totalScatterRate += mech->scatterRate(pe);
         // Cache cumulative scatterRate
         cached_cumulativeScatterRate[index++] = totalScatterRate;
      }
   }

   __host__ __device__ double MONSEL_MaterialScatterModel::randomMeanPathLength(ElectronT& pe)
   {
      setCache(&pe);
      /*
      * Return at most the chamber diameter for the free path. This is long
      * enough to guarantee the electron hits the chamber wall from any
      * starting position, provided it doesn't hit something else first. A long
      * path will most likely be shortened to a boundary by a routine that
      * computes the boundary distance as a multiple of this free path. If the
      * free path is too long (or infinite) the multiple will be too small,
      * leading to numerical imprecision.
      */
      const double maxFreePath = 2. * MonteCarloSS::ChamberRadius;
      if (totalScatterRate != 0.) {
         const double freepath = -::log(Random::random()) / totalScatterRate;
         if (freepath != freepath) printf("MONSEL_MaterialScatterModel::randomMeanPathLength: freepath is NaN\n");
         return freepath > maxFreePath ? maxFreePath : freepath;
      }
      /*
      * I still on very rare occasions get situations where an electron gets
      * stuck in an infinite loop. E.g., after several days and millions of
      * electrons I had one that struck within 5.e-8 nm of a corner. To avoid
      * having these lock up an otherwise good simulation, just drop the
      * electron.
      */
      if (pe.getStepCount() > 1000000)
         pe.setTrajectoryComplete(true);
      return maxFreePath;
   }

   __host__ __device__ ElectronT* MONSEL_MaterialScatterModel::scatter(ElectronT& pe)
   {
      /*
      * 1. Determine which of the scatter mechanisms active in this material
      * caused this scatter event. (Randomly choose one, with probability
      * weighted by scatter rate.) 2. Obtain a scattering angle appropriate to
      * this mechanism
      */

      const double eK = pe.getPreviousEnergy();
      if (eK != cached_eK) {
         /*
         * scattering is based on energy at conclusion of previous step, but
         * scatterRate methods used to set the cache use the current energy.
         * This is inelegant, but for now I hold my nose and ...
         */
         const double eKsaved = pe.getEnergy();
         pe.setEnergy(eK);
         setCache(&pe);
         pe.setEnergy(eKsaved);
      }

      /*
      * There is the possibility (e.g., vacuum) that no scatter mechanisms are
      * assigned to this material. In this case the totalScatterRate will be 0.
      * In this event we do nothing.
      */
      if (totalScatterRate == 0.)
         return NULL;

      // Find the scatter mechanism that produced this scattering event
      // Generate a random # between 0 and total cumulative scatter rate

      const double r = Random::random() * totalScatterRate;
      int index = 0; // Index is first index

      // Increment index and mechanism until cumulative scatter rate exceeds r
      while (cached_cumulativeScatterRate[index] < r)
         index++;

      return scatterArray[index]->scatter(&pe);
   }

   __host__ __device__ ElectronT* MONSEL_MaterialScatterModel::barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const
   {
      return barrierSM->barrierScatter(pe, nextRegion);
   }

   const BarrierScatterMechanismT* MONSEL_MaterialScatterModel::getBarrierSM() const
   {
      return barrierSM;
   }

   __host__ __device__ void MONSEL_MaterialScatterModel::setCSD(SlowingDownAlgT* csd)
   {
      this->csd = csd;
      // TODO I should clone it first. User is likely to use this in more than
      // one place.
      this->csd->setMaterial(mat);
   }

   SlowingDownAlgT* MONSEL_MaterialScatterModel::getCSD()
   {
      return csd;
   }

   __host__ __device__ double MONSEL_MaterialScatterModel::calculateEnergyLoss(double len, const ElectronT& pe) const
   {
      return csd->compute(len, &pe);
   }

   __host__ __device__ bool MONSEL_MaterialScatterModel::addScatterMechanism(ScatterMechanismT* mech)
   {
      // final ScatterMechanism mechCopy = mech.clone(); // Make a copy of the
      // scatter
      // mechanism
      // mechCopy.setMaterial(mat); // Initialize the copy for this material
      // if(scatterSet.add(mechCopy)) {
      if (amp::find(scatterArray.begin(), scatterArray.end(), mech) == scatterArray.end()) {
         scatterArray.push_back(mech);
         nscattermech = scatterArray.size();
         cached_cumulativeScatterRate.resize(nscattermech);
         cached_eK = -1.;// Force recompute of cache on next call
         return true;
      }
      return false;
   }

   //LinkedHashSet<ScatterMechanism> getScatterSet() {
   //   return (LinkedHashSet<ScatterMechanism>) Collections.unmodifiableSet(scatterSet);
   //   /*
   //   * TODO This may not be good enough. A caller can get the unmodifiable set
   //   * (a set of references to scattermechanisms) and call the init routine on
   //   * one or more of its elements with a different material than the one with
   //   * which it was originally initialized. This will change its internal
   //   * state. If I fix this I may not also need to override clone();
   //   */
   //}

   //bool MONSEL_MaterialScatterModel::removeScatterMechanism(ScatterMechanismT* mech)
   //{
   //   ScatterMechanismList::const_iterator itr = amp::find(scatterArray.begin(), scatterArray.end(), mech); // TODO: fix slow find
   //   if (itr == scatterArray.end()) return false;

   //   scatterArray.erase(itr);
   //   nscattermech = scatterArray.size();
   //   cached_cumulativeScatterRate.resize(nscattermech);
   //   cached_eK = -1.; // Force recompute of cache on next call
   //   return true;
   //}

   __host__ __device__ double MONSEL_MaterialScatterModel::getMinEforTracking() const
   {
      return minEforTracking;
   }

   __host__ __device__ void MONSEL_MaterialScatterModel::setMinEforTracking(double minEforTracking)
   {
      this->minEforTracking = minEforTracking;
   }
}