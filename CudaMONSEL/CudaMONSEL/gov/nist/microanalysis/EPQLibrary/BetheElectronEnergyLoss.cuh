// file: gov\nist\microanalysis\EPQLibrary\BetheElectronEnergyLoss.cuh

#ifndef _BETHE_ELECTRON_ENERGY_LOSS_CUH_
#define _BETHE_ELECTRON_ENERGY_LOSS_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AlgorithmClass.cuh"

namespace BetheElectronEnergyLoss
{
   class BetheElectronEnergyLoss : public AlgorithmClassT
   {
      __host__ __device__ void initializeDefaultStrategy() override;

   public:
      __host__ __device__ BetheElectronEnergyLoss(StringT name, const ReferenceT& ref);

      AlgorithmClassT const * const * getAllImplementations() const override;

      __host__ __device__ virtual float compute(const ElementT& elm, float eB) const = 0;
   };

   class JoyLuoBetheElectronEnergyLoss : public BetheElectronEnergyLoss
   {
   public:
      __host__ __device__ JoyLuoBetheElectronEnergyLoss();
      __host__ __device__ float compute(const ElementT& el, float eB) const override;

   private:
      VectorXf mK;
      const float K;
   };

   class Bethe30ModElectronEnergyLoss : public BetheElectronEnergyLoss
   {
   public:
      __host__ __device__ Bethe30ModElectronEnergyLoss();
      __host__ __device__ float compute(const ElementT& el, float eB) const override;

   private:
      const float K;
   };

   class Bethe30ElectronEnergyLoss : public BetheElectronEnergyLoss
   {
   public:
      __host__ __device__ Bethe30ElectronEnergyLoss();
      __host__ __device__ float compute(const ElementT& el, float eB) const override;

   private:
      const float K;
   };

   class StragglingModified : public BetheElectronEnergyLoss
   {
   public:
      __host__ __device__ StragglingModified(const BetheElectronEnergyLoss& base, float percent);
      __host__ __device__ float compute(const ElementT& el, float eB) const override;

   private:
      const BetheElectronEnergyLoss& mBethe;
      const float mPercent;
   };

   extern const BetheElectronEnergyLoss& JoyLuo1989;
   extern const BetheElectronEnergyLoss& Bethe1930;
   extern const BetheElectronEnergyLoss& Bethe1930Strict;

   extern __device__ const BetheElectronEnergyLoss* d_JoyLuo1989;

   extern __global__ void initCuda();
}

#endif