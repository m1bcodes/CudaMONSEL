// file: gov\nist\microanalysis\EPQLibrary\ScreenedRutherfordScatteringAngle.cuh
#ifndef _SCREENED_RUTHERFORD_SCATTERING_ANGLE_CUH_
#define _SCREENED_RUTHERFORD_SCATTERING_ANGLE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatter.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatterFactory.cuh"

namespace ScreenedRutherfordScatteringAngle
{
   class ScreenedRutherfordScatteringAngle : public RandomizedScatterT
   {
   public:
      __host__ __device__ explicit ScreenedRutherfordScatteringAngle(const ElementT& elm);

      StringT toString() const;

      __host__ __device__ const ElementT& getElement() const override;
      __host__ __device__ double totalCrossSection(const double energy) const override;
      __host__ __device__ double randomScatteringAngle(const double energy) const override;

   private:
      const ElementT& mElement;
   };

   class ScreenedRutherfordRandomizedScatterFactory : public RandomizedScatterFactoryT
   {
   public:
      __host__ __device__ ScreenedRutherfordRandomizedScatterFactory();

      __host__ __device__ const RandomizedScatterT& get(const ElementT& elm) const override;
   };

   extern const RandomizedScatterFactoryT& Factory;

   __host__ __device__ extern const ScreenedRutherfordScatteringAngle& getSRSA(int an);

   extern void init();
   extern __global__ void initCuda();
}

#endif