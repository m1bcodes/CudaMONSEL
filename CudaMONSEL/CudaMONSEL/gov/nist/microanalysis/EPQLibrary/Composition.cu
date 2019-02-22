#include "Composition.cuh"

namespace Composition
{
   __device__ const long long serialVersionUID = 0x42;
   __device__ const double OUT_OF_THIS_MANY_ATOMS = 1.0;

   __device__ void Composition::renormalize()
   {
      if (LinkedListKV::Size<Element::Element, UncertainValue2::UncertainValue2>(mConstituents) > 0) {
         mNormalization = UncertainValue2::ZERO;
         auto constituentHead = mConstituents;
         while (constituentHead != NULL) {
            auto uv = constituentHead->GetValue();
            if (uv.doubleValue() > 0.0) {
               mNormalization = UncertainValue2::add(mNormalization, uv);
            }
            constituentHead = constituentHead->GetNext();
         }
         mAtomicNormalization = UncertainValue2::ZERO;

         auto constituentsAtomicHead = mConstituentsAtomic;
         while (constituentHead != NULL) {
            auto uv = constituentsAtomicHead->GetValue();
            if (uv.doubleValue() > 0.0) {
               mAtomicNormalization = UncertainValue2::add(mAtomicNormalization, uv);
            }
            constituentsAtomicHead = constituentsAtomicHead->GetNext();
         }
      }
      else {
         mNormalization = UncertainValue2::ONE;
         mAtomicNormalization = UncertainValue2::ONE;
      }
      mMoleNorm = UncertainValue2::NaN;
   }
}
