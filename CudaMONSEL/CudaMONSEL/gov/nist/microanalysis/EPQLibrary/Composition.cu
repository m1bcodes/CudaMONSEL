#include "Composition.cuh"

extern __device__ float __int_as_float(int x);

namespace Composition
{
   __device__ const long long serialVersionUID = 0x42;
   __device__ const double OUT_OF_THIS_MANY_ATOMS = 1.0;

   __device__ void Composition::renormalize()
   {
      if (mConstituents != NULL) {
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

   __device__ Composition::Composition()
   {
      mHashCode = CUDART_INF_F;
      renormalize();
   }

   __device__ Composition::Composition(const Composition& comp)
   {
      replicate(comp);
   }

   __device__ Composition::Composition(Element::Element elms[], int elmsLen, double massFracs[], int massFracsLen)
   {
      if (elmsLen != massFracsLen) {
         printf("Composition::Composition: elmsLen != massFracsLen, (%d, %d)", elmsLen, massFracsLen);
      }
      for (int i = 0; i < elmsLen; ++i) {
         LinkedListKV::InsertHead(&mConstituents, elms[i], UncertainValue2::UncertainValue2(massFracs[i]));
      }
      recomputeStoiciometry();
      renormalize();
   }

   __device__ Composition::Composition(Element::Element elm)
   {
      LinkedListKV::InsertHead(&mConstituents, elm, UncertainValue2::ONE);
      recomputeStoiciometry();
      renormalize();
   }

   __device__ Composition::Composition(Element::Element elms[], int elmsLen, double massFracs[], int massFracsLen, char* name)
   {
      if (elmsLen == massFracsLen - 1) {
         double* wf = new double[elmsLen];
         double sum = 0.0;
         for (int i = 0; i < massFracsLen; ++i) {
            sum += massFracs[i];
            wf[i] = massFracs[i];
         }
         if (sum > 1.0) {
            printf("Composition::Composition: sum is greater than 1 (%lf)", sum);
         }
         wf[elmsLen - 1] = 1.0 - sum;
         massFracs = wf;
         delete[] wf;
         elmsLen = massFracsLen;
      }
      if (elmsLen != massFracsLen) {
         printf("Composition::Composition: elmsLen != massFracsLen, (%d, %d)", elmsLen, massFracsLen);
      }
      for (int i = 0; i < elmsLen; ++i) {
         if (massFracs[i] < 0.0) {
            printf("A mass fraction was less than zero while defining the material %s", name);
         }
         LinkedListKV::InsertHead(&mConstituents, elms[i], UncertainValue2::UncertainValue2(massFracs[i]));
      }
      mName = name;
      recomputeStoiciometry();
      renormalize();
   }

   __device__ Composition Composition::readResolve()
   {
      mHashCode = CUDART_INF_F;
      renormalize();
      return *this;
   }

   __device__ void Composition::replicate(Composition comp)
   {
      LinkedListKV::RemoveAll<Element::Element, UncertainValue2::UncertainValue2>(&mConstituents);
      LinkedListKV::RemoveAll<Element::Element, UncertainValue2::UncertainValue2>(&mConstituentsAtomic);
      auto constituentsHead = comp.mConstituents;
      while (constituentsHead != NULL) {
         LinkedListKV::InsertHead<Element::Element, UncertainValue2::UncertainValue2>(&mConstituents, constituentsHead->GetKey(), constituentsHead->GetValue());
         constituentsHead = constituentsHead->GetNext();
      }
      auto constituentsAtomicHead = comp.mConstituentsAtomic;
      while (constituentsAtomicHead != NULL) {
         LinkedListKV::InsertHead<Element::Element, UncertainValue2::UncertainValue2>(&mConstituentsAtomic, constituentsAtomicHead->GetKey(), constituentsAtomicHead->GetValue());
         constituentsAtomicHead = constituentsAtomicHead->GetNext();
      }
      mNormalization = comp.mNormalization;
      mAtomicNormalization = comp.mAtomicNormalization;
      mMoleNorm = comp.mMoleNorm;
      mHashCode = comp.mHashCode;
      mName = comp.mName;
   }

   __device__ LinkedList::Node<Element::Element>* Composition::getElementSet()
   {
      LinkedList::Node<Element::Element>* head = NULL;
      AdvancedLinkedList::AddAllKeys(&head, mConstituents, Element::AreEqual);
      return head;
   }

   __device__ LinkedList::Node<Element::Element>* Composition::getSortedElements()
   {
      LinkedList::Node<Element::Element>* res;
      AdvancedLinkedList::AddAllKeys(&res, mConstituents, Element::AreEqual);
      // TODO: need to implement quick qort for LinkedList
      //Collections.sort(res, new Comparator<Element>() {
      //   @Override
      //      public int compare(Element o1, Element o2) {
      //      return -Double.compare(weightFraction(o1, false), weightFraction(o2, false));
      //   }
      //});
      return res;
   }

   __device__ int Composition::getElementCount() {
      return LinkedListKV::Size(mConstituents);
   }

   __device__ void Composition::addElement(int atomicNo, double massFrac)
   {
      addElement(Element::byAtomicNumber(atomicNo), massFrac);
   }

   __device__ void Composition::addElement(int atomicNo, UncertainValue2::UncertainValue2 massFrac)
   {
      addElement(Element::byAtomicNumber(atomicNo), massFrac);
   }

   __device__ void Composition::addElement(Element::Element elm, double massFrac)
   {
      addElement(elm, UncertainValue2::UncertainValue2(massFrac));
   }

   __device__ double Composition::weightFraction(Element::Element elm, bool normalized)
   {
      UncertainValue2::UncertainValue2 d = LinkedListKV::GetValue(mConstituents, elm, Element::AreEqual);
      return *((int*)&d) != NULL ? (normalized ? normalize(d, mNormalization, true).doubleValue() : d.doubleValue()) : 0.0;
   }

   __device__ void Composition::recomputeStoiciometry()
   {
      mMoleNorm = UncertainValue2::ZERO;
      auto constituentsHead = mConstituents;
      while (constituentsHead != NULL) {
         mMoleNorm = UncertainValue2::add(mMoleNorm, UncertainValue2::multiply(1.0 / constituentsHead->GetKey().getAtomicWeight(), constituentsHead->GetValue()));
         constituentsHead = constituentsHead->GetNext();
      }

      LinkedListKV::RemoveAll<Element::Element, UncertainValue2::UncertainValue2>(&mConstituentsAtomic);
      LinkedList::Node<Element::Element>* constituentsAtomicKeysHead = NULL;
      AdvancedLinkedList::AddAllKeys(&constituentsAtomicKeysHead, mConstituentsAtomic, Element::AreEqual);
      auto constituentsAtomicKeysHeadItr = constituentsAtomicKeysHead;
      while (constituentsAtomicKeysHeadItr != NULL) {
         auto elm = constituentsAtomicKeysHeadItr->GetValue();
         UncertainValue2::UncertainValue2 moleFrac = (mMoleNorm.doubleValue() > 0.0 ? UncertainValue2::divide(LinkedListKV::GetValue<Element::Element, UncertainValue2::UncertainValue2>(mConstituents, elm, Element::AreEqual), UncertainValue2::multiply(elm.getAtomicWeight() / OUT_OF_THIS_MANY_ATOMS, mMoleNorm)) : UncertainValue2::ZERO);
         LinkedListKV::InsertHead<Element::Element, UncertainValue2::UncertainValue2>(&mConstituentsAtomic, elm, moleFrac);
         constituentsAtomicKeysHeadItr = constituentsAtomicKeysHeadItr->GetNext();
      }
      mOptimalRepresentation = Representation::WEIGHT_PCT;
      LinkedList::RemoveAll(&constituentsAtomicKeysHead);
   }

   __device__ void Composition::addElement(Element::Element elm, UncertainValue2::UncertainValue2 massFrac)
   {
      LinkedListKV::InsertHead<Element::Element, UncertainValue2::UncertainValue2>(&mConstituents, elm, massFrac);
      recomputeStoiciometry();
      renormalize();
   }

   __device__ UncertainValue2::UncertainValue2 Composition::weightFractionU(Element::Element elm, bool normalized)
   {
      return weightFractionU(elm, normalized, true);
   }

   __device__ UncertainValue2::UncertainValue2 Composition::weightFractionU(Element::Element elm, bool normalized, bool positiveOnly)
   {
      UncertainValue2::UncertainValue2 d = LinkedListKV::GetValue<Element::Element, UncertainValue2::UncertainValue2>(mConstituents, elm, Element::AreEqual);
      return *((int*)&d) != NULL ? (normalized ? normalize(d, mNormalization, positiveOnly) : d) : UncertainValue2::ZERO;
   }

   __device__ Composition positiveDefinite(Composition comp)
   {
      Composition res;
      auto elemSetHead = comp.getElementSet();
      while (elemSetHead != NULL) {
         auto elm = elemSetHead->GetValue();
         if (comp.weightFraction(elm, false) > 0.0) {
            res.addElement(elm, comp.weightFractionU(elm, false));
         }
      }
         
      return res;
   }

   __device__ void Composition::addElementByStoiciometry(Element::Element elm, UncertainValue2::UncertainValue2 moleFrac)
   {
      LinkedListKV::InsertHead<Element::Element, UncertainValue2::UncertainValue2>(&mConstituentsAtomic, elm, moleFrac);
      recomputeWeightFractions();
      renormalize();
   }

   __device__ void Composition::addElementByStoiciometry(Element::Element elm, double moleFrac)
   {
      addElementByStoiciometry(elm, UncertainValue2::UncertainValue2(moleFrac));
   }

   __device__ void Composition::clear()
   {
      LinkedListKV::RemoveAll(&mConstituents);
      LinkedListKV::RemoveAll(&mConstituentsAtomic);
      mNormalization = UncertainValue2::ONE;
      mAtomicNormalization = UncertainValue2::ONE;
      mMoleNorm = UncertainValue2::NaN;
   }

   __device__ void Composition::defineByWeightFraction(Element::Element elms[], int elmsLen, double wgtFracs[], int wgtFracsLen)
   {
      clear();
      if (elmsLen != wgtFracsLen) {
         printf("Composition::defineByWeightFraction1: elmsLen != wgtFracsLen (%d, %d)", elmsLen, wgtFracsLen);
      }
      for (int i = 0; i < elmsLen; ++i) {
         LinkedListKV::InsertHead(&mConstituents, elms[i], UncertainValue2::UncertainValue2(wgtFracs[i]));
      }
      recomputeStoiciometry();
      renormalize();
   }

   __device__ void Composition::defineByWeightFraction(Element::Element elms[], int elmsLen, UncertainValue2::UncertainValue2 wgtFracs[], int wgtFracsLen)
   {
      clear();
      if (elmsLen != wgtFracsLen) {
         printf("Composition::defineByWeightFraction2: elmsLen != wgtFracsLen (%d, %d)", elmsLen, wgtFracsLen);
      }
      for (int i = 0; i < elmsLen; ++i) {
         LinkedListKV::InsertHead(&mConstituents, elms[i], wgtFracs[i]);
      }
      recomputeStoiciometry();
      renormalize();
   }

   __device__ void Composition::defineByMoleFraction(Element::Element elms[], int elmsLen, double moleFracs[], int moleFracsLen)
   {
      clear();
      if (elmsLen != moleFracsLen) {
         printf("Composition::defineByWeightFraction: elmsLen != moleFracsLen (%d, %d)", elmsLen, moleFracsLen);
      }
      int mfSum = 0;
      for (int k = 0; k < moleFracsLen; ++k) {
         mfSum += moleFracs[k];
      }
      for (int k = 0; k < moleFracsLen; ++k) {
         moleFracs[k] /= mfSum;
      }
      for (int i = 0; i < moleFracsLen; ++i) {
         LinkedListKV::InsertHead(&mConstituentsAtomic, elms[i], UncertainValue2::UncertainValue2(moleFracs[i]));
      }
      recomputeWeightFractions();
      renormalize();
   }

   __device__ UncertainValue2::UncertainValue2 Composition::atomicPercentU(Element::Element elm) {
      return atomicPercentU(elm, true);
   }

   __device__ UncertainValue2::UncertainValue2 Composition::atomicPercentU(Element::Element elm, bool positiveOnly)
   {
      UncertainValue2::UncertainValue2 o = LinkedListKV::GetValue<Element::Element, UncertainValue2::UncertainValue2>(mConstituentsAtomic, elm, Element::AreEqual);
      return *((int*)&o) != NULL ? normalize(o, mAtomicNormalization, positiveOnly) : UncertainValue2::ZERO;
   }

   __device__ void Composition::recomputeWeightFractions()
   {
      UncertainValue2::UncertainValue2 totalWgt = UncertainValue2::ZERO;
      LinkedList::Node<Element::Element>* constituentsAtomicKeysHead, * constituentsAtomicKeysHead1, * constituentsAtomicKeysHead2;
      AdvancedLinkedList::AddAllKeys(&constituentsAtomicKeysHead, mConstituentsAtomic, Element::AreEqual);
      constituentsAtomicKeysHead1 = constituentsAtomicKeysHead;
      while (constituentsAtomicKeysHead1 != NULL) {
         auto elm = constituentsAtomicKeysHead1->GetValue();
         totalWgt = UncertainValue2::add(totalWgt, UncertainValue2::multiply(elm.getAtomicWeight(), atomicPercentU(elm)));
         constituentsAtomicKeysHead1 = constituentsAtomicKeysHead1->GetNext();
      }
         
      LinkedListKV::RemoveAll<Element::Element, UncertainValue2::UncertainValue2>(&mConstituents);
      constituentsAtomicKeysHead2 = constituentsAtomicKeysHead;
      while (constituentsAtomicKeysHead2 != NULL) {
         auto elm = constituentsAtomicKeysHead2->GetValue();
         UncertainValue2::UncertainValue2 wgtFrac = UncertainValue2::multiply(elm.getAtomicWeight(), UncertainValue2::divide(atomicPercentU(elm), totalWgt));
         LinkedListKV::InsertHead<Element::Element, UncertainValue2::UncertainValue2>(&mConstituents, elm, wgtFrac);
         constituentsAtomicKeysHead2 = constituentsAtomicKeysHead2->GetNext();
      }
      mOptimalRepresentation = Representation::STOICIOMETRY;
      LinkedList::RemoveAll(&constituentsAtomicKeysHead);
   }

   __device__ UncertainValue2::UncertainValue2 normalize(UncertainValue2::UncertainValue2 val, UncertainValue2::UncertainValue2 norm, bool positive)
   {
      UncertainValue2::UncertainValue2 res;
      if (norm.doubleValue() > 0.0) {
         res = UncertainValue2::divide(val, norm);
      }
      else {
         res = val;
      }
      return positive ? UncertainValue2::positiveDefinite(res) : res;
   }
}
