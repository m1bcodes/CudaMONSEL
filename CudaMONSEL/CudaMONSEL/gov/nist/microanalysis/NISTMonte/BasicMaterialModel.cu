//// File: gov\nist\microanalysis\NISTMonte\BasicMaterialModel.cu
//
#include "gov\nist\microanalysis\NISTMonte\BasicMaterialModel.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AlgorithmUser.cuh"
//#include "gov\nist\microanalysis\EPQLibrary\EPQException.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Gas.cuh"
#include "gov\nist\microanalysis\EPQLibrary\GasScatteringCrossSection.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
#include "gov\nist\microanalysis\EPQLibrary\NISTMottScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatterFactory.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\BetheElectronEnergyLoss.cuh"

#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"

#include "gov\nist\microanalysis\Utility\Math2.cuh"

#include "Amphibian\random.cuh"

namespace BasicMaterialModel
{
   BasicMaterialModel::BasicMaterialModel(const MaterialT& mat) : mMaterial(mat), minEforTracking(ToSI::eV(50.0))
   {
      initializeDefaultStrategy();
      //if (mMaterial instanceof Gas) {
      //   addDefaultAlgorithm(RandomizedScatterFactory.class, GasScatteringCrossSection.Factory);
      //}
   }
   
   __host__ __device__ const MaterialT& BasicMaterialModel::getMaterial() const
   {
      return mMaterial;
   }

   __host__ __device__ double BasicMaterialModel::randomMeanPathLength(ElectronT& pe)
   {
      // Ref: Heinrich 1981 p 458
      double kE = pe.getEnergy();
      double minMfp = 1.0;
      const ElementT* bestEl = nullptr;
      double den = mMaterial.getDensity();
      const RandomizedScatterFactoryT* rsf = (RandomizedScatterFactoryT*)getAlgorithm("RandomizedScatterFactory");
      if (!(rsf != NULL)) printf("BasicMaterialModel::randomMeanPathLength: rsf == NULL\n");
      for (auto &el : mMaterial.getElementSet()) {
         const double mfp = (el->getMass() * Random::expRand()) / (den * mMaterial.weightFraction(*el, true) * rsf->get(*el).totalCrossSection(kE));
         if (mfp < minMfp) {
            minMfp = mfp;
            bestEl = el;
         }
      }
      pe.setScatteringElement(bestEl);
      return minMfp;
   }

   __host__ __device__ ElectronT* BasicMaterialModel::scatter(ElectronT& pe)
   {
      const ElementT* se = pe.getScatteringElement();
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      if ((se != nullptr) && !(*se == *Element::dNone))
#else
      if ((se != nullptr) && !(*se == Element::None))
#endif
      {
         const RandomizedScatterFactoryT* rsf = (RandomizedScatterFactoryT*)getAlgorithm("RandomizedScatterFactory");
         if (rsf == nullptr) printf("BasicMaterialModel::scatter: rsf == NULL\n");
         const double alpha = rsf->get(*se).randomScatteringAngle(pe.getEnergy());
         const double beta = 2.0 * Math2::PI * Random::random();
         // Update the primary electron's direction angles
         // We could pe.setEnergy() here, but this is an elastic model
         // so it is not necessary.
         pe.updateDirection(alpha, beta);
      }
      return nullptr; // --no SE generation in this model
   }

   __host__ __device__ ElectronT* BasicMaterialModel::barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const
   {
      pe->setCurrentRegion(nextRegion);
      pe->setScatteringElement(nullptr);
      return nullptr;
   }

   __host__ __device__ double BasicMaterialModel::calculateEnergyLoss(double len, const ElectronT& pe) const
   {
      // See Heinrich 1981 pp 226-227
      const double kE = pe.getEnergy();
      double res = 0.0;
      for (const auto &el : mMaterial.getElementSet())
         res += ::AlgorithmUser::getDefaultBetheEnergyLoss()->compute(*el, kE) * mMaterial.weightFraction(*el, true);
      return res * mMaterial.getDensity() * len;
   }

   __host__ __device__ double BasicMaterialModel::getMinEforTracking() const
   {
      return minEforTracking;
   }

   __host__ __device__ void BasicMaterialModel::setMinEforTracking(double minEforTracking)
   {
      this->minEforTracking = minEforTracking;
   }

   __host__ __device__ void BasicMaterialModel::initializeDefaultStrategy()
   {
      // addDefaultAlgorithm(RandomizedScatterFactory.class,
      // ScreenedRutherfordScatteringAngle.Factory);
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      addDefaultAlgorithm("RandomizedScatterFactory", NISTMottScatteringAngle::d_Factory);
#else
      addDefaultAlgorithm("RandomizedScatterFactory", &NISTMottScatteringAngle::Factory);
#endif
   }
}