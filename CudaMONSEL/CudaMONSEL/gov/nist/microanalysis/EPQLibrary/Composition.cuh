#ifndef _COMPOSITION_CUH_
#define _COMPOSITION_CUH_

#include "..\..\..\..\Amphibian\LinkedList.cuh"
#include "..\..\..\..\Amphibian\String.cuh"

#include "Element.cuh"
#include "..\Utility\UncertainValue2.cuh"

#include <cuda_runtime.h>

namespace Composition
{
   enum Representation
   {
      UNDETERMINED, WEIGHT_PCT, STOICIOMETRY
   };

   class Composition
   {
   public:


   private:
      LinkedListKV::Node<Element::Element, UncertainValue2::UncertainValue2>* mConstituents = NULL;
      LinkedListKV::Node<Element::Element, UncertainValue2::UncertainValue2>* mConstituentsAtomic = NULL;

      UncertainValue2::UncertainValue2 mNormalization = UncertainValue2::ONE;
      UncertainValue2::UncertainValue2 mAtomicNormalization = UncertainValue2::ONE;
      String::String mName;
      Representation mOptimalRepresentation = Representation::UNDETERMINED;
      UncertainValue2::UncertainValue2 mMoleNorm = UncertainValue2::NaN;

   protected:
      __device__ void renormalize();

      int mHashCode; // = CUDART_INF_F;
   };
}
#endif
