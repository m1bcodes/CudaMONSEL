// package gov.nist.microanalysis.EPQLibrary.BrowningEmpiricalCrossSection.cuh

#ifndef _BROWNING_EMPIRICAL_CROSS_SECTION_CUH_
#define _BROWNING_EMPIRICAL_CROSS_SECTION_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace BrowningEmpiricalCrossSection
{
   class BrowningEmpiricalCrossSection
   {
   public:
      __host__ __device__ BrowningEmpiricalCrossSection(const ElementT& elm);

      const ElementT& getElement() const;
      double totalCrossSection(double energy) const;
      double randomScatteringAngle(double energy) const;

   private:
      const ElementT& mElement;

      const double mZp17;
      const double mZp2;
      const double mZp3;
   };

   extern const BrowningEmpiricalCrossSection& getBECS(int an);
}

#endif