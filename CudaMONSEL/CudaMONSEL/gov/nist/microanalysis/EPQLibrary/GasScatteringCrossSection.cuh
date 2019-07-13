#ifndef _GAS_SCATTERING_CROSS_SECTION_CUH_
#define _GAS_SCATTERING_CROSS_SECTION_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatter.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatterFactory.cuh"

namespace GasScatteringCrossSection
{
   class GasScatteringCrossSection : public RandomizedScatterT
   {
   public:
      __host__ __device__ explicit GasScatteringCrossSection(const ElementT& elm);
      //GasScatteringCrossSection(const GasScatteringCrossSection& gscs);

      const RandomizedScatterT& getElasticModel();
      __host__ __device__ double ratioInelasticOverElastic() const;

      __host__ __device__ const ElementT& getElement() const override;
      __host__ __device__ double totalCrossSection(const double energy) const override;
      __host__ __device__ double randomScatteringAngle(const double energy) const override;

   private:
      const ElementT& mElement;
      const ScreenedRutherfordScatteringAngleT& mElastic;
   };
   
   class GasScatteringRandomizedScatterFactory : public RandomizedScatterFactoryT
   {
   public:
      __host__ __device__ GasScatteringRandomizedScatterFactory();

      __host__ __device__ const RandomizedScatterT& get(const ElementT& elm) const override;
   };

   extern const RandomizedScatterFactoryT& Factory;

   extern void init();
   extern __global__ void initCuda();
}

#endif