// package gov.nist.microanalysis.EPQLibrary.BrowningEmpiricalCrossSection.cuh

#ifndef _BROWNING_EMPIRICAL_CROSS_SECTION_CUH_
#define _BROWNING_EMPIRICAL_CROSS_SECTION_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace BrowningEmpiricalCrossSection
{
   class BrowningEmpiricalCrossSection
   {
   public:
      __host__ __device__ explicit BrowningEmpiricalCrossSection(const ElementT& elm);

      __host__ __device__ const ElementT& getElement() const;
      __host__ __device__ double totalCrossSection(const double energy) const;
      __host__ __device__ double randomScatteringAngle(const double energy) const;

   private:
      const ElementT& mElement;

      const double mZp17;
      const double mZp2;
      const double mZp3;
   };

   __host__ __device__ extern const BrowningEmpiricalCrossSection& getBECS(int an);

   extern void init();

   extern __global__ void initCuda();
}

#endif