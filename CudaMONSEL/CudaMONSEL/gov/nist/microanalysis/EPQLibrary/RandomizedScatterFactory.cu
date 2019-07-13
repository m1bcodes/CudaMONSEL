#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatterFactory.cuh"
#include "gov\nist\microanalysis\EPQLibrary\GasScatteringCrossSection.cuh"
#include "gov\nist\microanalysis\EPQLibrary\NISTMottScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ScreenedRutherfordScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\CzyzewskiMottScatteringAngle.cuh"

namespace RandomizedScatterFactory
{
   __host__ __device__ RandomizedScatterFactory::RandomizedScatterFactory(StringT name, const ReferenceT& ref) :
      AlgorithmClassT("Scatter factory", name, ref)
   {
   }

   const AlgorithmClassT* mImplementations[] =
   {
      &GasScatteringCrossSection::Factory,
      &NISTMottScatteringAngle::Factory,
      &ScreenedRutherfordScatteringAngle::Factory,
      &CzyzewskiMottScatteringAngle::Factory
   };

   AlgorithmClassT const * const * RandomizedScatterFactory::getAllImplementations() const
   {      
      return mImplementations;
   }

   __host__ __device__ void RandomizedScatterFactory::initializeDefaultStrategy()
   {
   }
}
