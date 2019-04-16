//#include "Composition.cuh"
//#include "gov\nist\microanalysis\Utility\CudaUtil.h"
//
//#include <curand.h>
//#include <curand_kernel.h>
//
////#include <time.h>
////#include <stdlib.h>
//
//extern __device__ float __int_as_float(int x);
//
//namespace Composition
//{
//   __device__ static const long long serialVersionUID = 0x42;
//   __device__ static const double OUT_OF_THIS_MANY_ATOMS = 1.0;
//
//   __device__ static bool doubleCmp (double& a, double& b) { return a == b; };
//
//   //__device__ void Composition::renormalize()
//   //{
//   //   if (!mConstituents.IsEmpty()) {
//   //      mNormalization = UncertainValue2::ZERO();
//   //      Map::Iterator<Element::Element, UncertainValue2::UncertainValue2> itr(mConstituents);
//   //      while (itr.HasNext()) {
//   //         auto uv = itr.GetValue();
//   //         if (uv.doubleValue() > 0.0) {
//   //            mNormalization = UncertainValue2::add(mNormalization, uv);
//   //         }
//   //         itr.Next();
//   //      }
//
//   //      mAtomicNormalization = UncertainValue2::ZERO();
//
//   //      auto constituentsAtomicHead = mConstituentsAtomic;
//   //      while (constituentsAtomicHead != NULL) {
//   //         auto uv = constituentsAtomicHead->GetValue();
//   //         if (uv.doubleValue() > 0.0) {
//   //            mAtomicNormalization = UncertainValue2::add(mAtomicNormalization, uv);
//   //         }
//   //         constituentsAtomicHead = constituentsAtomicHead->GetNext();
//   //      }
//   //   }
//   //   else {
//   //      mNormalization = UncertainValue2::ONE();
//   //      mAtomicNormalization = UncertainValue2::ONE();
//   //   }
//   //   mMoleNorm = UncertainValue2::NaN();
//   //}
//
//   __device__ Composition::Composition()
//      //mConstituents(DefaultHasher, Element::AreEqual, UncertainValue2::AreEqual),
//      //mConstituentsAtomic(NULL),
//      //mNormalization(UncertainValue2::UncertainValue2(1)),
//      //mAtomicNormalization(UncertainValue2::UncertainValue2(1)),
//      //mName(""),
//      //mOptimalRepresentation(Representation::UNDETERMINED),
//      //mMoleNorm(UncertainValue2::NaN())
//   {
//      mHashCode = CUDART_INF_F;
//      //renormalize();
//   }
//
//   //__device__ Composition::Composition(Composition& comp) : mConstituents(DefaultHasher, Element::AreEqual, UncertainValue2::AreEqual), mConstituentsAtomic(NULL)
//   //{
//   //   replicate(comp);
//   //}
//
//   //__device__ Composition::Composition(Element::Element elms[], int elmsLen, double massFracs[], int massFracsLen) : mConstituents(DefaultHasher, Element::AreEqual, UncertainValue2::AreEqual)
//   //{
//   //   if (elmsLen != massFracsLen) {
//   //      printf("Composition::Composition: elmsLen != massFracsLen, (%d, %d)", elmsLen, massFracsLen);
//   //   }
//   //   for (int i = 0; i < elmsLen; ++i) {
//   //      mConstituents.Put(elms[i], UncertainValue2::UncertainValue2(massFracs[i]));
//   //   }
//   //   recomputeStoiciometry();
//   //   renormalize();
//   //}
//
//   //__device__ Composition::Composition(Element::Element elm) : mConstituents(DefaultHasher, Element::AreEqual, UncertainValue2::AreEqual)
//   //{
//   //   mConstituents.Put(elm, UncertainValue2::ONE());
//   //   recomputeStoiciometry();
//   //   renormalize();
//   //}
//
//   //__device__ Composition::Composition(Element::Element elms[], int elmsLen, double massFracs[], int massFracsLen, char* name) : mConstituents(DefaultHasher, Element::AreEqual, UncertainValue2::AreEqual)
//   //{
//   //   if (elmsLen == massFracsLen - 1) {
//   //      double* wf = new double[elmsLen];
//   //      double sum = 0.0;
//   //      for (int i = 0; i < massFracsLen; ++i) {
//   //         sum += massFracs[i];
//   //         wf[i] = massFracs[i];
//   //      }
//   //      if (sum > 1.0) {
//   //         printf("Composition::Composition: sum is greater than 1 (%lf)", sum);
//   //      }
//   //      wf[elmsLen - 1] = 1.0 - sum;
//   //      massFracs = wf;
//   //      delete[] wf;
//   //      elmsLen = massFracsLen;
//   //   }
//   //   if (elmsLen != massFracsLen) {
//   //      printf("Composition::Composition: elmsLen != massFracsLen, (%d, %d)", elmsLen, massFracsLen);
//   //   }
//   //   for (int i = 0; i < elmsLen; ++i) {
//   //      if (massFracs[i] < 0.0) {
//   //         printf("A mass fraction was less than zero while defining the material %s", name);
//   //      }
//   //      mConstituents.Put(elms[i], UncertainValue2::UncertainValue2(massFracs[i]));
//   //   }
//   //   mName = name;
//   //   recomputeStoiciometry();
//   //   renormalize();
//   //}
//
//   //__device__ Composition::~Composition()
//   //{
//   //   LinkedListKV::RemoveAll(&mConstituentsAtomic);
//   //}
//
////
////   __device__ Composition Composition::readResolve()
////   {
////      mHashCode = CUDART_INF_F;
////      renormalize();
////      return *this;
////   }
////
//   //__device__ void Composition::replicate(Composition& comp)
//   //{
//   //   mConstituents.RemoveAll();
//   //   LinkedListKV::RemoveAll<Element::Element, UncertainValue2::UncertainValue2>(&mConstituentsAtomic);
//
//   //   Map::Iterator<Element::Element, UncertainValue2::UncertainValue2> itr(comp.mConstituents);
//   //   while (itr.HasNext()) {
//   //      mConstituents.Put(itr.GetKey(), itr.GetValue());
//   //      itr.Next();
//   //   }
//
//   //   auto constituentsAtomicHead = comp.mConstituentsAtomic;
//   //   while (constituentsAtomicHead != NULL) {
//   //      LinkedListKV::InsertHead<Element::Element, UncertainValue2::UncertainValue2>(&mConstituentsAtomic, constituentsAtomicHead->GetKey(), constituentsAtomicHead->GetValue());
//   //      constituentsAtomicHead = constituentsAtomicHead->GetNext();
//   //   }
//   //   mNormalization = comp.mNormalization;
//   //   mAtomicNormalization = comp.mAtomicNormalization;
//   //   mMoleNorm = comp.mMoleNorm;
//   //   mHashCode = comp.mHashCode;
//   //   mName = comp.mName;
//   //}
////
////   __device__ Set::Set<Element::Element> Composition::getElementSet()
////   {
////      return mConstituents.GetKeys();
////   }
////
////   __device__ Set::Set<Element::Element> Composition::getSortedElements()
////   {
////      auto keys = mConstituents.GetKeys();
////      // TODO: need to implement quick qort for LinkedList, or BST set/map
////      //Collections.sort(res, new Comparator<Element>() {
////      //   @Override
////      //      public int compare(Element o1, Element o2) {
////      //      return -Double.compare(weightFraction(o1, false), weightFraction(o2, false));
////      //   }
////      //});
////      return keys;
////   }
////
////   __device__ int Composition::getElementCount()
////   {
////      return mConstituents.Size();
////   }
////
////   __device__ void Composition::addElement(int atomicNo, double massFrac)
////   {
////      addElement(Element::byAtomicNumber(atomicNo), massFrac);
////   }
////
////   __device__ void Composition::addElement(int atomicNo, UncertainValue2::UncertainValue2 massFrac)
////   {
////      addElement(Element::byAtomicNumber(atomicNo), massFrac);
////   }
////
////   __device__ void Composition::addElement(Element::Element elm, double massFrac)
////   {
////      addElement(elm, UncertainValue2::UncertainValue2(massFrac));
////   }
////
////   __device__ double Composition::weightFraction(Element::Element elm, bool normalized)
////   {
////      UncertainValue2::UncertainValue2 d = mConstituents.GetValue(elm);
////      return *((int*)&d) != NULL ? (normalized ? normalize(d, mNormalization, true).doubleValue() : d.doubleValue()) : 0.0;
////   }
////
//   //__device__ void Composition::recomputeStoiciometry()
//   //{
//   //   mMoleNorm = UncertainValue2::ZERO();
//   //   Map::Iterator<Element::Element, UncertainValue2::UncertainValue2> itr(mConstituents);
//   //   while (itr.HasNext()) {
//   //      mMoleNorm = UncertainValue2::add(mMoleNorm, UncertainValue2::multiply(1.0 / itr.GetKey().getAtomicWeight(), itr.GetValue()));
//   //      itr.Next();
//   //   }
//
//   //   LinkedListKV::RemoveAll(&mConstituentsAtomic);
//   //   LinkedList::Node<Element::Element>* constituentsAtomicKeysHead = NULL;
//   //   AdvancedLinkedList::AddAllKeys(&constituentsAtomicKeysHead, mConstituentsAtomic, Element::AreEqual);
//   //   auto constituentsAtomicKeysHeadItr = constituentsAtomicKeysHead;
//   //   while (constituentsAtomicKeysHeadItr != NULL) {
//   //      auto elm = constituentsAtomicKeysHeadItr->GetValue();
//   //      UncertainValue2::UncertainValue2 moleFrac = (mMoleNorm.doubleValue() > 0.0 ? UncertainValue2::divide(mConstituents.GetValue(elm), UncertainValue2::multiply(elm.getAtomicWeight() / OUT_OF_THIS_MANY_ATOMS, mMoleNorm)) : UncertainValue2::ZERO());
//   //      LinkedListKV::InsertHead<Element::Element, UncertainValue2::UncertainValue2>(&mConstituentsAtomic, elm, moleFrac);
//   //      constituentsAtomicKeysHeadItr = constituentsAtomicKeysHeadItr->GetNext();
//   //   }
//   //   mOptimalRepresentation = Representation::WEIGHT_PCT;
//   //   LinkedList::RemoveAll(&constituentsAtomicKeysHead);
//   //}
////
////   __device__ void Composition::addElement(Element::Element elm, UncertainValue2::UncertainValue2 massFrac)
////   {
////      mConstituents.Put(elm, massFrac);
////      recomputeStoiciometry();
////      renormalize();
////   }
////
////   __device__ UncertainValue2::UncertainValue2 Composition::weightFractionU(Element::Element elm, bool normalized)
////   {
////      return weightFractionU(elm, normalized, true);
////   }
////
////   __device__ UncertainValue2::UncertainValue2 Composition::weightFractionU(Element::Element elm, bool normalized, bool positiveOnly)
////   {
////      UncertainValue2::UncertainValue2 d = mConstituents.GetValue(elm);
////      return *((int*)&d) != NULL ? (normalized ? normalize(d, mNormalization, positiveOnly) : d) : UncertainValue2::ZERO();
////   }
////
////   __device__ Composition positiveDefinite(Composition& comp)
////   {
////      Composition res;
////      auto elemSet = comp.getElementSet();
////      Set::Iterator<Element::Element> itr(elemSet);
////      while (itr.HasNext()) {
////         auto elm = itr.GetValue();
////         if (comp.weightFraction(elm, false) > 0.0) {
////            res.addElement(elm, comp.weightFractionU(elm, false));
////         }
////         itr.Next();
////      }
////      return res;
////   }
////
////   __device__ void Composition::addElementByStoiciometry(Element::Element elm, UncertainValue2::UncertainValue2 moleFrac)
////   {
////      LinkedListKV::InsertHead<Element::Element, UncertainValue2::UncertainValue2>(&mConstituentsAtomic, elm, moleFrac);
////      recomputeWeightFractions();
////      renormalize();
////   }
////
////   __device__ void Composition::addElementByStoiciometry(Element::Element elm, double moleFrac)
////   {
////      addElementByStoiciometry(elm, UncertainValue2::UncertainValue2(moleFrac));
////   }
////
////   __device__ void Composition::clear()
////   {
////      mConstituents.RemoveAll();
////      LinkedListKV::RemoveAll(&mConstituentsAtomic);
////      mNormalization = UncertainValue2::ONE();
////      mAtomicNormalization = UncertainValue2::ONE();
////      mMoleNorm = UncertainValue2::NaN();
////   }
////
////   __device__ void Composition::defineByWeightFraction(Element::Element elms[], int elmsLen, double wgtFracs[], int wgtFracsLen)
////   {
////      clear();
////      if (elmsLen != wgtFracsLen) {
////         printf("Composition::defineByWeightFraction1: elmsLen != wgtFracsLen (%d, %d)", elmsLen, wgtFracsLen);
////      }
////      for (int i = 0; i < elmsLen; ++i) {
////         mConstituents.Put(elms[i], UncertainValue2::UncertainValue2(wgtFracs[i]));
////      }
////      recomputeStoiciometry();
////      renormalize();
////   }
////
////   __device__ void Composition::defineByWeightFraction(Element::Element elms[], int elmsLen, UncertainValue2::UncertainValue2 wgtFracs[], int wgtFracsLen)
////   {
////      clear();
////      if (elmsLen != wgtFracsLen) {
////         printf("Composition::defineByWeightFraction2: elmsLen != wgtFracsLen (%d, %d)", elmsLen, wgtFracsLen);
////      }
////      for (int i = 0; i < elmsLen; ++i) {
////         mConstituents.Put(elms[i], wgtFracs[i]);
////      }
////      recomputeStoiciometry();
////      renormalize();
////   }
////
////   __device__ void Composition::defineByMoleFraction(Element::Element elms[], int elmsLen, double moleFracs[], int moleFracsLen)
////   {
////      clear();
////      if (elmsLen != moleFracsLen) {
////         printf("Composition::defineByWeightFraction: elmsLen != moleFracsLen (%d, %d)", elmsLen, moleFracsLen);
////      }
////      int mfSum = 0;
////      for (int k = 0; k < moleFracsLen; ++k) {
////         mfSum += moleFracs[k];
////      }
////      for (int k = 0; k < moleFracsLen; ++k) {
////         moleFracs[k] /= mfSum;
////      }
////      for (int i = 0; i < moleFracsLen; ++i) {
////         LinkedListKV::InsertHead(&mConstituentsAtomic, elms[i], UncertainValue2::UncertainValue2(moleFracs[i]));
////      }
////      recomputeWeightFractions();
////      renormalize();
////   }
////
////   __device__ UncertainValue2::UncertainValue2 Composition::atomicPercentU(Element::Element elm) {
////      return atomicPercentU(elm, true);
////   }
////
////   __device__ UncertainValue2::UncertainValue2 Composition::atomicPercentU(Element::Element elm, bool positiveOnly)
////   {
////      UncertainValue2::UncertainValue2 o = LinkedListKV::GetValue<Element::Element, UncertainValue2::UncertainValue2>(mConstituentsAtomic, elm, Element::AreEqual);
////      return *((int*)&o) != NULL ? normalize(o, mAtomicNormalization, positiveOnly) : UncertainValue2::ZERO();
////   }
////
////   __device__ void Composition::recomputeWeightFractions()
////   {
////      UncertainValue2::UncertainValue2 totalWgt = UncertainValue2::ZERO();
////      LinkedList::Node<Element::Element>* constituentsAtomicKeysHead, * constituentsAtomicKeysHead1, * constituentsAtomicKeysHead2;
////      AdvancedLinkedList::AddAllKeys(&constituentsAtomicKeysHead, mConstituentsAtomic, Element::AreEqual);
////      constituentsAtomicKeysHead1 = constituentsAtomicKeysHead;
////      while (constituentsAtomicKeysHead1 != NULL) {
////         auto elm = constituentsAtomicKeysHead1->GetValue();
////         totalWgt = UncertainValue2::add(totalWgt, UncertainValue2::multiply(elm.getAtomicWeight(), atomicPercentU(elm)));
////         constituentsAtomicKeysHead1 = constituentsAtomicKeysHead1->GetNext();
////      }
////
////      mConstituents.RemoveAll();
////      constituentsAtomicKeysHead2 = constituentsAtomicKeysHead;
////      while (constituentsAtomicKeysHead2 != NULL) {
////         auto elm = constituentsAtomicKeysHead2->GetValue();
////         UncertainValue2::UncertainValue2 wgtFrac = UncertainValue2::multiply(elm.getAtomicWeight(), UncertainValue2::divide(atomicPercentU(elm), totalWgt));
////         mConstituents.Put(elm, wgtFrac);
////         constituentsAtomicKeysHead2 = constituentsAtomicKeysHead2->GetNext();
////      }
////      mOptimalRepresentation = Representation::STOICIOMETRY;
////      LinkedList::RemoveAll(&constituentsAtomicKeysHead);
////   }
////
////   __device__ UncertainValue2::UncertainValue2 normalize(UncertainValue2::UncertainValue2& val, UncertainValue2::UncertainValue2& norm, bool positive)
////   {
////      UncertainValue2::UncertainValue2 uv;
////      if (norm.doubleValue() > 0.0) {
////         uv = UncertainValue2::divide(val, norm);
////      }
////      else {
////         uv = val;
////      }
////      return positive ? UncertainValue2::positiveDefinite(uv) : uv;
////   }
////
////   __device__ void Composition::setOptimalRepresentation(Representation opt)
////   {
////      switch (opt)
////      {
////      case UNDETERMINED:
////         break;
////      case WEIGHT_PCT:
////         recomputeStoiciometry();
////         break;
////      case STOICIOMETRY:
////         recomputeWeightFractions();
////         break;
////      }
////   }
////
////   __device__ Set::Set<Element::Element> elementSet(Composition compositions[], int len)
////   {
////      Set::Set<Element::Element> elms(DefaultHasher, Element::AreEqual);
////      for (int i = 0; i < len; ++i) {
////         elms.Add(compositions[i].getElementSet());
////      }
////      return elms;
////   }
////
////   __device__ void Composition::defineByMaterialFraction(Composition compositions[], int compLen, double matFracs[], int matFracsLen)
////   {
////      if (compLen != matFracsLen) {
////         printf("Composition::defineByMaterialFraction: lengths are different (%d, %d)", compLen , matFracsLen);
////      }
////      clear();
////      auto elms = elementSet(compositions, compLen);
////      int len = elms.Size();
////      Element::Element* newElms = new Element::Element[len];
////      UncertainValue2::UncertainValue2* frac = new UncertainValue2::UncertainValue2[len];
////
////      int ji = 0;
////      Set::Iterator<Element::Element> elmsItr(elms);
////      while (elmsItr.HasNext()) {
////         auto el = elmsItr.GetValue();
////         UncertainValue2::UncertainValue2 sum = UncertainValue2::ZERO();
////         for (int i = 0; i < compLen; ++i) {
////            sum = UncertainValue2::add(sum, UncertainValue2::multiply(matFracs[i], compositions[i].weightFractionU(el, true)));
////         }
////         frac[ji] = sum;
////         newElms[ji] = el;
////         ++ji;
////         elmsItr.Next();
////      }
////      defineByWeightFraction(newElms, len, frac, len);
////
////      elms.RemoveAll();
////      delete[] newElms;
////      delete[] frac;
////   }
////
////   __device__ void Composition::removeElement(Element::Element el)
////   {
////      if (mConstituents.ContainsKey(el)) {
////         mConstituents.Remove(el);
////         LinkedListKV::Remove(&mConstituentsAtomic, el, Element::AreEqual);
////         // Don't recomputeStoiciometry or recomputeWeightFractions
////         renormalize();
////      }
////   }
////
////   __device__ bool Composition::containsElement(Element::Element el)
////   {
////      return (mConstituents.ContainsKey(el) && (mConstituents.GetValue(el).doubleValue() > 0.0));
////   }
////
////   __device__ bool Composition::containsAll(LinkedList::Node<Element::Element>* elms)
////   {
////      while (elms != NULL) {
////         auto elm = elms->GetValue();
////         if (!containsElement(elm)) {
////            return false;
////         }
////         elms = elms->GetNext();
////      }
////      return true;
////   }
////
////   __device__ double Composition::atomicPercent(Element::Element elm)
////   {
////      return atomicPercentU(elm).doubleValue();
////   }
////
////   __device__ UncertainValue2::UncertainValue2 Composition::stoichiometryU(Element::Element elm)
////   {
////      UncertainValue2::UncertainValue2 o = LinkedListKV::GetValue(mConstituentsAtomic, elm, Element::AreEqual);
////      return *((int*)&o) != NULL ? o : UncertainValue2::ZERO();
////   }
////
////   __device__ double Composition::stoichiometry(Element::Element elm)
////   {
////      return stoichiometryU(elm).doubleValue();
////   }
////
////   __device__ double Composition::atomsPerKg(Element::Element elm, bool normalized)
////   {
////      return weightFraction(elm, normalized) / elm.getMass();
////   }
////
////   __device__ UncertainValue2::UncertainValue2 Composition::atomsPerKgU(Element::Element elm, bool normalized)
////   {
////      return UncertainValue2::multiply(1.0 / elm.getMass(), weightFractionU(elm, normalized));
////   }
////
////   __device__ UncertainValue2::UncertainValue2 Composition::weightAvgAtomicNumberU()
////   {
////      UncertainValue2::UncertainValue2 res = UncertainValue2::ZERO();
////      Map::Iterator<Element::Element, UncertainValue2::UncertainValue2> constituentsItr(mConstituents);
////      while (constituentsItr.HasNext()) {
////         Element::Element elm = constituentsItr.GetKey();
////         res = UncertainValue2::add(res, UncertainValue2::multiply(elm.getAtomicNumber(), weightFractionU(elm, true)));
////         constituentsItr.Next();
////      }
////      return res;
////   }
////
////   __device__ double Composition::weightAvgAtomicNumber()
////   {
////      return weightAvgAtomicNumberU().doubleValue();
////   }
////
////   __device__ double Composition::sumWeightFraction()
////   {
////      double sum = 0.0;
////      Map::Iterator<Element::Element, UncertainValue2::UncertainValue2> constituentsItr(mConstituents);
////      while (constituentsItr.HasNext()) {
////         auto uv = constituentsItr.GetValue();
////         if (uv.doubleValue() > 0.0) {
////            sum += uv.doubleValue();
////         }
////         constituentsItr.Next();
////      }
////      return sum;
////   }
////
////   __device__ UncertainValue2::UncertainValue2 Composition::sumWeightFractionU()
////   {
////      UncertainValue2::UncertainValue2 res = UncertainValue2::ZERO();
////      Map::Iterator<Element::Element, UncertainValue2::UncertainValue2> constituentsItr(mConstituents);
////      while (constituentsItr.HasNext()) {
////         auto val = constituentsItr.GetValue();
////         if (val.doubleValue() > 0.0) {
////            res = UncertainValue2::add(res, val);
////         }
////         constituentsItr.Next();
////      }
////      return res;
////   }
////
////   ////__device__ String::String Composition::toString()
////   ////{
////   ////   if ((*((int*)&mName) == NULL) || (mName.Length() == 0)) {
////   ////      return descriptiveString(false);
////   ////   }
////   ////   return mName;
////   ////}
////
////   ////__device__ String::String Composition::stoichiometryString()
////   ////{
////   ////   final StringBuffer sb = new StringBuffer();
////   ////   final NumberFormat nf = new HalfUpFormat("0.####");
////   ////   for (final Element elm : getElementSet()) {
////   ////      final UncertainValue2 d0 = atomicPercentU(elm);
////   ////      if (sb.length() > 1)
////   ////         sb.append(",");
////   ////      sb.append(elm.toAbbrev());
////   ////      sb.append("(");
////   ////      sb.append(d0.format(nf));
////   ////      sb.append(" atoms)");
////   ////   }
////   ////   return sb.toString();
////   ////}
////
////   ////__device__ String::String weightPercentString(bool normalize)
////   ////{
////   ////   final StringBuffer sb = new StringBuffer();
////   ////   final NumberFormat nf = new HalfUpFormat("0.0000");
////   ////   for (final Element elm : getElementSet()) {
////   ////      final UncertainValue2 d0 = weightFractionU(elm, normalize);
////   ////      if (sb.length() > 1)
////   ////         sb.append(",");
////   ////      sb.append(elm.toAbbrev());
////   ////      sb.append("(");
////   ////      sb.append(d0.format(nf));
////   ////      sb.append(" mass frac)");
////   ////   }
////   ////   if (!normalize) {
////   ////      sb.append(",\u03A3=");
////   ////      sb.append(sumWeightPercentU().format(nf));
////   ////   }
////   ////   return sb.toString();
////   ////}
////
////   ////__device__ String::String descriptiveString(bool normalize)
////   ////{
////   ////   StringBuffer sb = new StringBuffer();
////   ////   if ((mName != null) && (mName.length() > 0))
////   ////      sb.append(mName + " = ");
////   ////   sb.append("[");
////   ////   if (mOptimalRepresentation == Representation.STOICIOMETRY)
////   ////      sb.append(stoichiometryString());
////   ////   else
////   ////      sb.append(weightPercentString(normalize));
////   ////   sb.append("]");
////   ////   return sb.toString();
////   ////}
////
////   //__device__ Element::Element Composition::getNthElementByWeight(int n)
////   //{
////   //   LinkedListKV::Node<UncertainValue2::UncertainValue2, Element::Element>* tm = NULL;
////   //   auto constituentsItr = mConstituents;
////   //   while (constituentsItr != NULL) {
////   //      auto el = constituentsItr->GetKey();
////   //      UncertainValue2::UncertainValue2 wf = weightFractionU(el, true);
////   //      // Add hoc mechanism to handle the case in which multiple elements are
////   //      // present in the same weightPct.
////   //      curandState state;
////   //      curand_init(1, 0, 0, &state);
////   //      while (LinkedListKV::ContainsKey(tm, wf, UncertainValue2::AreEqual)) {
////   //         wf = UncertainValue2::add(1.0e-10 * curand_uniform(&state), wf);
////   //      }
////   //      LinkedListKV::InsertHead(&tm, wf, el);
////   //      constituentsItr = constituentsItr->GetNext();
////   //   }
////   //   int j = 0;
////   //   auto tmItr = tm;
////   //   while (tmItr != NULL) {
////   //      ++j;
////   //      if (j == LinkedListKV::Size(mConstituents) - n) {
////   //         return tmItr->GetValue();
////   //      }
////   //      tmItr = tmItr->GetNext();
////   //   }
////   //   return Element::None;
////   //}
////
////   //__device__ Element::Element Composition::getNthElementByAtomicFraction(int n)
////   //{
////   //   LinkedListKV::Node<UncertainValue2::UncertainValue2, Element::Element>* tm = NULL;
////   //   auto constituentsItr = mConstituents;
////   //   while (constituentsItr != NULL) {
////   //      auto el = constituentsItr->GetKey();
////   //      auto mf = atomicPercentU(el);
////   //      curandState state;
////   //      curand_init(1, 0, 0, &state);
////   //      while (LinkedListKV::ContainsKey(tm, mf, UncertainValue2::AreEqual)) {
////   //         mf = UncertainValue2::add(1.0e-10 * curand_uniform(&state), mf);
////   //      }
////   //      LinkedListKV::InsertHead(&tm, mf, el);
////   //      constituentsItr = constituentsItr->GetNext();
////   //   }
////   //   int j = 0;
////   //   auto tmItr = tm;
////   //   while (tmItr != NULL) {
////   //      ++j;
////   //      if (j == n) {
////   //         return tmItr->GetValue();
////   //      }
////   //      tmItr = tmItr->GetNext();
////   //   }
////   //   return Element::None;
////   //}
////
////   __device__ void Composition::setName(String::String name)
////   {
////      mName = name;
////   }
////
////   __device__ String::String Composition::getName()
////   {
////      return mName;
////   }
////
////   __device__ int Composition::compareTo(Composition& comp)
////   {
////      if (this == &comp) {
////         return 0;
////      }
////      // hashers have to be the same also
////      Map::Iterator<Element::Element, UncertainValue2::UncertainValue2> i(mConstituents);
////      Map::Iterator<Element::Element, UncertainValue2::UncertainValue2> j(comp.mConstituents);
////      while (i.HasNext() && j.HasNext()) {
////         int zi = i.GetKey().getAtomicNumber();
////         int zj = j.GetKey().getAtomicNumber();
////         if (zi < zj) {
////            return +1;
////         }
////         else if (zi > zj) {
////            return -1;
////         }
////         else {
////            UncertainValue2::UncertainValue2 ci = i.GetValue();
////            UncertainValue2::UncertainValue2 cj = j.GetValue();
////            if (ci.lessThan(cj)) {
////               return -1;
////            }
////            else if (ci.greaterThan(cj)) {
////               return +1;
////            }
////         }
////         i.Next();
////         j.Next();
////      }
////      if (!(i.HasNext())) {
////         return +1;
////      }
////      if (!(j.HasNext())) {
////         return -1;
////      }
////      return 0;
////   }
////
////   __device__ Composition Composition::asComposition()
////   {
////      Composition res;
////      res.replicate(*this);
////      return res;
////   }
////
////   __device__ Composition Composition::clone()
////   {
////      return asComposition();
////   }
////
////   __device__ UncertainValue2::UncertainValue2 Composition::differenceU(Composition& comp)
////   {
////      // assert (comp.getElementCount() == this.getElementCount());
////      UncertainValue2::UncertainValue2 delta = UncertainValue2::ZERO();
////      Set::Set<Element::Element> allElms(DefaultHasher, Element::AreEqual);
////      allElms.Add(getElementSet());
////      allElms.Add(comp.getElementSet());
////      Set::Iterator<Element::Element> allElmsItr(allElms);
////      while (allElmsItr.HasNext()) {
////         auto el = allElmsItr.GetValue();
////         delta = UncertainValue2::add(delta, UncertainValue2::sqr(UncertainValue2::subtract(comp.weightFractionU(el, false), weightFractionU(el, false))));
////         allElmsItr.Next();
////      }
////      return UncertainValue2::multiply(1.0 / allElms.Size(), delta).sqrt();
////   }
////
////   __device__ double Composition::difference(Composition& comp)
////   {
////      return differenceU(comp).doubleValue();
////   }
////
////   __device__ Representation Composition::getOptimalRepresentation()
////   {
////      return mOptimalRepresentation;
////   }
////
////   //__device__ int Composition::hashCode()
////   //{
////   //   if (mHashCode == INT_MAX) {
////   //      int result = 1;
////   //      int PRIME = 31;
////   //      result = PRIME * result + mConstituents.hashCode();
////   //      result = PRIME * result + mConstituentsAtomic.hashCode();
////   //      result = PRIME * result + ((mName == null) ? 0 : mName.hashCode());
////   //      long temp;
////   //      temp = mNormalization.hashCode();
////   //      result = PRIME * result + (int)(temp ^ (temp >> > 32));
////   //      result = PRIME * result + ((mOptimalRepresentation == null) ? 0 : mOptimalRepresentation.hashCode());
////   //      if (result == Integer.MAX_VALUE)
////   //         result = Integer.MIN_VALUE;
////   //      mHashCode = result;
////   //   }
////   //   return mHashCode;
////   //}
////
////   __device__ Map::Map<Element::Element, UncertainValue2::UncertainValue2> Composition::GetConstituents()
////   {
////      return mConstituents;
////   }
////
////   //__device__ bool Composition::AreEqualConstituents(LinkedListKV::Node<Element::Element, UncertainValue2::UncertainValue2>* a, LinkedListKV::Node<Element::Element, UncertainValue2::UncertainValue2>* b)
////   //{
////   //   //if (!LinkedListKV::AreEquivalentSets(a, b, Element::AreEqual, UncertainValue2::AreEqual)) {
////   //   //   return false;
////   //   //}
////   //   return true;
////   //}
////
////   __device__ bool Composition::equals(Composition& obj)
////   {
////      if (this == &obj) {
////         return true;
////      }
////      if (*((int*)&obj) == NULL) {
////         return false;
////      }
////      
////      //auto a = IsSet(mConstituents, Element::AreEqual, [](UnceratinValue2::UnceratinValue2, UnceratinValue2::UnceratinValue2) { return true; });
////      
////      //if (!LinkedListKV::AreEquivalentSets(GetConstituents(), obj.GetConstituents(), Element::AreEqual, UncertainValue2::AreEqual)) {
////      //   return false;
////      //}
////      if (GetConstituents() == obj.GetConstituents()) {
////         return false;
////      }
////      // TODO: fix below
////      //if (!LinkedListKV::AreEquivalentSets(mConstituentsAtomic, obj.mConstituentsAtomic, Element::AreEqual, UncertainValue2::AreEqual)) {
////      //   return false;
////      //}
////      if ((*((int*)&mName) != NULL) && (*((int*)&(obj.mName)) != NULL) && (!(mName == obj.mName))) {
////         return false;
////      }
////      if (!mNormalization.equals(obj.mNormalization)) {
////         return false;
////      }
////      if (!(mOptimalRepresentation == obj.mOptimalRepresentation)) {
////         return false;
////      }
////      return true;
////   }
////
////   __device__ bool Composition::almostEquals(Composition& other, double tol)
////   {
////      if (this == &other) {
////         return true;
////      }
////      if (*((int*)&other) == NULL) {
////         return false;
////      }
////      if (abs(mNormalization.doubleValue() - other.mNormalization.doubleValue()) > tol) {
////         return false;
////      }
////      Set::Set<Element::Element> allElms(DefaultHasher, Element::AreEqual);
////      allElms.Add(other.getElementSet());
////      allElms.Add(getElementSet());
////      Set::Iterator<Element::Element> allElmsItr(allElms);
////      while (allElmsItr.HasNext()) {
////         Element::Element elm = allElmsItr.GetValue();
////         {
////            UncertainValue2::UncertainValue2 uv1 = weightFractionU(elm, false);
////            UncertainValue2::UncertainValue2 uv2 = other.weightFractionU(elm, false);
////            if ((*((int*)&uv1) == NULL) || (*((int*)&uv2) == NULL)) {
////               return false;
////            }
////            if ((abs(uv1.doubleValue() - uv2.doubleValue()) > tol)
////               || (abs(uv1.uncertainty() - uv2.uncertainty()) > tol)) {
////               return false;
////            }
////         }
////         {
////            UncertainValue2::UncertainValue2 uv1 = this->atomicPercentU(elm);
////            UncertainValue2::UncertainValue2 uv2 = other.atomicPercentU(elm);
////            if ((*((int*)&uv1) == NULL) || (*((int*)&uv2) == NULL)) {
////               return false;
////            }
////            if ((abs(uv1.doubleValue() - uv2.doubleValue()) > tol)
////               || (abs(uv1.uncertainty() - uv2.uncertainty()) > tol)) {
////               return false;
////            }
////         }
////         allElmsItr.Next();
////      }
////      return true;
////   }
////
////   __device__ Map::Map<Element::Element, double> Composition::absoluteError(Composition& std, bool normalize)
////   {
////      Set::Set<Element::Element> elms(DefaultHasher, Element::AreEqual);
////      elms.Add(std.getElementSet());
////      elms.Add(getElementSet());
////      Map::Map<Element::Element, double> res(DefaultHasher, Element::AreEqual, [](double& a, double& b) {return a == b; });
////      Set::Iterator<Element::Element> elmsItr(elms);
////      while (elmsItr.HasNext()) {
////         auto elm = elmsItr.GetValue();
////         double u = weightFractionU(elm, normalize).doubleValue();
////         double s = std.weightFractionU(elm, normalize).doubleValue();
////         res.Put(elm, s != 0.0 ? (u - s) / s : (u == 0.0 ? 0.0 : 1.0));
////         elmsItr.Next();
////      }
////      return res;
////   }
////
////   __device__ Map::Map<Element::Element, double> Composition::relativeError(Composition& std, bool normalize)
////   {
////      Set::Set<Element::Element> elms(DefaultHasher, Element::AreEqual);
////      elms.Add(std.getElementSet());
////      elms.Add(getElementSet());
////      Map::Map<Element::Element, double> res(DefaultHasher, Element::AreEqual, [](double& a, double& b) { return a == b; });
////      Set::Iterator<Element::Element> elmsItr(elms);
////      while (elmsItr.HasNext()) {
////         auto elm = elmsItr.GetValue();
////         double u = weightFractionU(elm, normalize).doubleValue();
////         double s = std.weightFractionU(elm, normalize).doubleValue();
////         res.Put(elm, u - s);
////         elmsItr.Next();
////      }
////      return res;
////   }
////
////   __device__ bool Composition::isUncertain()
////   {
////      switch (mOptimalRepresentation) {
////      case WEIGHT_PCT:
////      case UNDETERMINED:
////         Map::Iterator<Element::Element, UncertainValue2::UncertainValue2> itr(mConstituents);
////         while (itr.HasNext()) {
////            auto v = itr.GetValue();
////            if (v.isUncertain()) {
////               return true;
////            }
////            itr.Next();
////         }
////         break;
////      case STOICIOMETRY:
////         auto j = mConstituentsAtomic;
////         while (j != NULL) {
////            auto v = j->GetValue();
////            if (v.isUncertain()) {
////               return true;
////            }
////            j = j->GetNext();
////         }
////         break;
////      }
////      return false;
////   }
////
////   __device__ UncertainValue2::UncertainValue2 Composition::meanAtomicNumberU()
////   {
////      UncertainValue2::UncertainValue2 res = UncertainValue2::ZERO();
////      auto elms = getElementSet();
////      Set::Iterator<Element::Element> itr(elms);
////      while (itr.HasNext()) {
////         auto elm = itr.GetValue();
////         res = UncertainValue2::add(res, UncertainValue2::multiply(elm.getAtomicNumber(), atomicPercentU(elm)));
////         itr.Next();
////      }
////      return res;
////   }
////
////   __device__ double Composition::meanAtomicNumber()
////   {
////      double res = 0.0;
////      auto elms = getElementSet();
////      Set::Iterator<Element::Element> itr(elms);
////      while (itr.HasNext()) {
////         auto elm = itr.GetValue();
////         res += elm.getAtomicNumber() * atomicPercent(elm);
////         itr.Next();
////      }
////      return res;
////   }
////
////   __device__ void Composition::forceNormalization()
////   {
////      UncertainValue2::UncertainValue2 norm = sumWeightFractionU();
////      Map::Map<Element::Element, UncertainValue2::UncertainValue2> newConst(DefaultHasher, Element::AreEqual, UncertainValue2::AreEqual);
////      Map::Iterator<Element::Element, UncertainValue2::UncertainValue2> itr(mConstituents);
////      while (itr.HasNext()) {
////         newConst.Put(itr.GetKey(), norm.doubleValue() > 0.0 ? UncertainValue2::divide(itr.GetValue(), norm) : UncertainValue2::ZERO());
////         itr.Next();
////      }
////      mConstituents.RemoveAll();
////      mConstituents = newConst;
////      mOptimalRepresentation = Representation::WEIGHT_PCT;
////      renormalize();
////   }
////
////   __device__ Composition parseGlass(char str[], int numlines)
////   {
////      Composition result;
////      int pos = 0;
////      int c = 0;
////      for (int n = 0; n < numlines; ++n) {
////         char line[256];
////         int k = 0;
////         while (str[c] != '\n') {
////            line[k] = str[c];
////            ++c;
////            ++k;
////         }
////         ++c;
////         if (pos == 0) {
////            if (String::StartsWith(line, "NBS GLASS K ")) {
////               result.setName("K" "NBS GLASS K");
////            }
////            else if (String::StartsWith(line, "CATIO")) {
////               pos = 1;
////            }
////         }
////         else if (pos == 1) {
////            if (String::StartsWith(line, "AVERAGE ATOMIC NUMBER")) {
////               pos = 2;
////            }
////            else {
////               char elmName[4];
////               char elmWgtPct[16];
////               int r = 0, l = 0;
////               int tabCount = 0;
////               while (true) {
////                  if (line[r] == '\t') {
////                     l = 0;
////                     ++tabCount;
////                  }
////                  if (tabCount > 6) {
////                     break;
////                  }
////                  if (tabCount == 0) {
////                     elmName[l] = line[r];
////                     ++l;
////                  }
////                  if (tabCount == 5) {
////                     elmWgtPct[l] = line[r];
////                     ++l;
////                  }
////                  ++r;
////               }
////
////               Element::Element elm = Element::byName(elmName);
////               double wgtPct = String::AToF<double>(elmWgtPct);
////               result.addElement(elm, wgtPct / 100.0);
////            }
////         }
////         else if (pos == 2) {
////            if (String::StartsWith(line, "WEIGHT PERCENT OXYGEN")) {
////               char oWgtPctStr[256];
////               int r = 0, l = 0;
////               int tabCount = 0;
////               while (true) {
////                  if (line[r] == '\t') {
////                     ++tabCount;
////                  }
////                  if (tabCount > 1) {
////                     break;
////                  }
////                  if (tabCount == 1) {
////                     oWgtPctStr[l] = line[r];
////                     ++l;
////                  }
////                  ++r;
////               }
////
////               double oWgtPct = String::AToF<double>(oWgtPctStr);
////               result.addElement(Element::O, oWgtPct / 100.0);
////               break;
////            }
////         }
////      }
////      return result;
////   }
////
////   __device__ Composition Composition::randomize(double offset, double proportional)
////   {
////      curandState state;
////      curand_init(2, 0, 0, &state);
////      Composition res;
////      auto elms = getElementSet();
////      Set::Iterator<Element::Element> itr(elms);
////      while (itr.HasNext()) {
////         Element::Element elm = itr.GetValue();
////         double w = weightFraction(elm, false);
////         double v = w + w * curand_normal(&state) * proportional + offset * curand_normal(&state);
////         v = v > 0.0 ? v : 0.0;
////         v = v < 1.1 ? v : 1.1;
////         res.addElement(elm, v);
////         itr.Next();
////      }
////      return res;
////   }
////
////   int DIM = 9;
////   __device__ long PROJECTORS[100]; // = createProjectors(2762689630628022905L);
////
////   __device__ long mIndexHashS = INT_MAX;
////   __device__ long mIndexHashL = INT_MAX;
////
////   void createProjectors(long seed)
////   {
////      srand(time(NULL));
////
////      long* hPROJECTORS = new long[100];
////      LinkedList::Node<long>* eval = NULL;
////      for (int j = 0; j < 100; ++j) {
////         long tmp;
////         do {
////            long mult = 1;
////            tmp = 0;
////            for (int i = 0; i < DIM; ++i, mult *= 10) {
////               double r = (double)rand() / (double)RAND_MAX;
////               tmp += r * 2 * mult;
////            }
////         } while (LinkedList::Exists<long>(eval, tmp, [](long& a, long& b) { return a == b; }));
////         hPROJECTORS[j] = tmp;
////         LinkedList::InsertHead(&eval, tmp);
////      }
////      LinkedList::RemoveAll(&eval);
////      checkCudaErrors(cudaMemcpyToSymbol(PROJECTORS, &hPROJECTORS, sizeof(long) * 100));
////      delete[] hPROJECTORS;
////   }
////
////   __device__ long Composition::indexHashCodeS()
////   {
////      if (mIndexHashS == INT_MAX) {
////         long res = 0;
////         auto elms = getElementSet();
////         Set::Iterator<Element::Element> itr(elms);
////         while (itr.HasNext()) {
////            Element::Element elm = itr.GetValue();
////            int v = (int)sqrt(100.0 * weightFraction(elm, false));
////            v = v > 0 ? v : 0;
////            v = v < 10 ? v : 0;
////            res += v * PROJECTORS[elm.getAtomicNumber()];
////            itr.Next();
////         }
////         mIndexHashS = res;
////      }
////      return mIndexHashS;
////   }
////
////   __device__ long Composition::indexHashCodeL()
////   {
////      if (mIndexHashL == INT_MAX) {
////         long res = 0;
////         auto elms = getElementSet();
////         Set::Iterator<Element::Element> itr(elms);
////         while (itr.HasNext()) {
////            Element::Element elm = itr.GetValue();
////            int v = (int)(10.0 * weightFraction(elm, false));
////            v = v > 0 ? v : 0;
////            v = v < 10 ? v : 0;
////            res += v * PROJECTORS[elm.getAtomicNumber()];
////            itr.Next();
////         }
////         mIndexHashL = res;
////      }
////      return mIndexHashL;
////   }
//}
