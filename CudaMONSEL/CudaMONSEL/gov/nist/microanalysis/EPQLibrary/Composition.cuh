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
      __device__ Composition();
      __device__ Composition(const Composition& comp);
      __device__ Composition(Element::Element elms[], int elmsLen, double massFracs[], int massFracsLen);
      __device__ Composition(Element::Element elm);
      __device__ Composition(Element::Element elms[], int elmsLen, double massFracs[], int massFracsLen, char* name);

      __device__ LinkedList::Node<Element::Element>* getElementSet();
      __device__ LinkedList::Node<Element::Element>* getSortedElements();
      __device__ int getElementCount();
      __device__ void addElement(int atomicNo, double massFrac);
      __device__ void addElement(int atomicNo, UncertainValue2::UncertainValue2 massFrac);
      __device__ void addElement(Element::Element elm, double massFrac);
      __device__ double weightFraction(Element::Element elm, bool normalized);
      __device__ void addElement(Element::Element elm, UncertainValue2::UncertainValue2 massFrac);
      __device__ UncertainValue2::UncertainValue2 weightFractionU(Element::Element elm, bool normalized);
      __device__ UncertainValue2::UncertainValue2 weightFractionU(Element::Element elm, bool normalized, bool positiveOnly);
      __device__ void addElementByStoiciometry(Element::Element elm, UncertainValue2::UncertainValue2 moleFrac);
      __device__ void addElementByStoiciometry(Element::Element elm, double moleFrac);
      __device__ UncertainValue2::UncertainValue2 atomicPercentU(Element::Element elm);
      __device__ UncertainValue2::UncertainValue2 atomicPercentU(Element::Element elm, bool positiveOnly);
      __device__ void defineByWeightFraction(Element::Element elms[], int elmsLen, double wgtFracs[], int wgtFracsLen);
      __device__ void defineByWeightFraction(Element::Element elms[], int elmsLen, UncertainValue2::UncertainValue2 wgtFracs[], int wgtFracsLen);
      __device__ void defineByWeightFraction(LinkedListKV::Node<Element::Element, double> map);

      template<typename T>
      __device__ Element::Element GetElementBy(T k)
      {
         return NULL;
      }

      template<>
      __device__ Element::Element GetElementBy<int>(int k)
      {
         return Element::byAtomicNumber(k);
      }

      template<>
      __device__ Element::Element GetElementBy<String::String>(String::String s)
      {
         return Element::byName(s.Get());
      }

      template<>
      __device__ Element::Element GetElementBy<Element::Element>(Element::Element e)
      {
         return e;
      }

      template<typename T>
      __device__ void defineByWeightFraction(LinkedListKV::Node<T, double>* map)
      {
         while (map != NULL) {
            double wp = map->GetValue();
            T key = map->GetKey();
            Element::Element elm = GetElementBy<T>(key);
            if (*((int*)&elm) == NULL) {
               printf("Composition::defineByWeightFraction: wrong type");
            }
            if (elm.getAtomicNumber() == Element::elmNone) {
               printf("Composition::defineByWeightFraction: bad element");
            }
            if ((*((int*)&elm) != NULL) && (elm.getAtomicNumber() == Element::elmNone)) {
               LinkedListKV::InsertHead(&mConstituents, elm, UncertainValue2::UncertainValue2(wp));
            }
         }
         recomputeStoiciometry();
         renormalize();
      }

      __device__ void defineByMoleFraction(Element::Element elms[], int elmsLen, double moleFracs[], int moleFracsLen);


   private:
      __device__ Composition readResolve();
      __device__ void recomputeStoiciometry();
      __device__ void recomputeWeightFractions();

      LinkedListKV::Node<Element::Element, UncertainValue2::UncertainValue2>* mConstituents = NULL;
      LinkedListKV::Node<Element::Element, UncertainValue2::UncertainValue2>* mConstituentsAtomic = NULL;

      UncertainValue2::UncertainValue2 mNormalization = UncertainValue2::ONE;
      UncertainValue2::UncertainValue2 mAtomicNormalization = UncertainValue2::ONE;
      String::String mName;
      Representation mOptimalRepresentation = Representation::UNDETERMINED;
      UncertainValue2::UncertainValue2 mMoleNorm = UncertainValue2::NaN;

   protected:
      __device__ void renormalize();
      __device__ void replicate(Composition comp);
      __device__ void clear();

      int mHashCode; // = CUDART_INF_F;
   };

   __device__ Composition positiveDefinite(Composition comp);
   __device__ UncertainValue2::UncertainValue2 normalize(UncertainValue2::UncertainValue2 val, UncertainValue2::UncertainValue2 norm, bool positive);


}
#endif
