#include "Composition.cuh"

//#include <time.h>
//#include <stdlib.h>

namespace Composition
{
   static const long long serialVersionUID = 0x42;
   static const double OUT_OF_THIS_MANY_ATOMS = 1.0;
   
   void Composition::renormalize()
   {
      if (!mConstituents.empty()) {
         mNormalization = UncertainValue2::ZERO();
         for (auto e : mConstituents) {
            auto uv = e.second;
            if (uv.doubleValue() > 0.0) {
               mNormalization = UncertainValue2::add(mNormalization, uv);
            }
         }

         mAtomicNormalization = UncertainValue2::ZERO();
         for (auto e : mConstituentsAtomic) {
            auto uv = e.second;
            if (uv.doubleValue() > 0.0) {
               mAtomicNormalization = UncertainValue2::add(mAtomicNormalization, uv);
            }
         }
      }
      else {
         mNormalization = UncertainValue2::ONE();
         mAtomicNormalization = UncertainValue2::ONE();
      }
      mMoleNorm = UncertainValue2::NaN();
   }

   Composition::Composition() :
      mNormalization(UncertainValue2::UncertainValue2(1)),
      mAtomicNormalization(UncertainValue2::UncertainValue2(1)),
      mName(""),
      mOptimalRepresentation(Representation::UNDETERMINED),
      mMoleNorm(UncertainValue2::NaN())
   {
      renormalize();
   }

   Composition::Composition(const Composition& comp)
   {
      replicate(comp);
   }

   Composition::Composition(const Element::Element* elms[], int elmsLen, double massFracs[], int massFracsLen)
   {
      if (elmsLen != massFracsLen) {
         printf("Composition::Composition: elmsLen != massFracsLen, (%d, %d)", elmsLen, massFracsLen);
      }
      for (int i = 0; i < elmsLen; ++i) {
         mConstituents.insert(std::make_pair(elms[i], UncertainValue2::UncertainValue2(massFracs[i])));
      }
      recomputeStoiciometry();
      renormalize();
   }

   Composition::Composition(const Element::Element& elm)
   {
      mConstituents.insert(std::make_pair(&elm, UncertainValue2::ONE()));
      recomputeStoiciometry();
      renormalize();
   }

   Composition::Composition(const Element::Element* elms[], int elmsLen, double massFracs[], int massFracsLen, char const* name)
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
         mConstituents.insert(std::make_pair(elms[i], UncertainValue2::UncertainValue2(massFracs[i])));
      }
      mName = name;
      recomputeStoiciometry();
      renormalize();
   }

   Composition::~Composition()
   {
   }

   bool Composition::operator==(const Composition& obj) const
   {
      if (this == &obj) {
         return true;
      }

      if (!sameConstituents(obj.mConstituents)) {
         return false;
      }

      if (!sameConstituentsAtomic(obj.mConstituentsAtomic)) {
         return false;
      }
      if (!(mName == obj.mName)) {
         return false;
      }

      if (!(mNormalization == obj.mNormalization)) {
         return false;
      }
      if (!(mOptimalRepresentation == obj.mOptimalRepresentation)) {
         return false;
      }
      return true;
   }

   void Composition::operator=(const Composition& comp)
   {
      replicate(comp);
   }

   Composition Composition::readResolve()
   {
      renormalize();
      return *this;
   }

   void Composition::replicate(const Composition& comp)
   {
      if (&comp == this) return;

      mConstituents.clear();
      mConstituentsAtomic.clear();

      for (auto itr : comp.mConstituents) {
         mConstituents.insert(std::make_pair(itr.first, itr.second));
      }

      auto constituentsAtomicHead = comp.mConstituentsAtomic;
      for (auto ca : comp.mConstituentsAtomic) {
         mConstituentsAtomic.insert(std::make_pair(ca.first, ca.second));
      }
      mNormalization = comp.mNormalization;
      mAtomicNormalization = comp.mAtomicNormalization;
      mMoleNorm = comp.mMoleNorm;
      mName = comp.mName;
      mOptimalRepresentation = comp.mOptimalRepresentation;
   }

   Element::UnorderedSetT Composition::getElementSet() const
   {
      Element::UnorderedSetT elmset;
      for (auto c : mConstituents) {
         elmset.insert(c.first);
      }
      return elmset;
   }

   Element::OrderedSetT Composition::getSortedElements() const
   {
      Element::OrderedSetT elmset;
      for (auto c : mConstituents) {
         elmset.insert(c.first);
      }
      return elmset;
   }

   int Composition::getElementCount() const
   {
      return mConstituents.size();
   }

   void Composition::addElement(int atomicNo, double massFrac)
   {
      addElement(Element::byAtomicNumber(atomicNo), massFrac);
   }

   void Composition::addElement(int atomicNo, const UncertainValue2::UncertainValue2 massFrac)
   {
      addElement(Element::byAtomicNumber(atomicNo), massFrac);
   }

   void Composition::addElement(const Element::Element& elm, double massFrac)
   {
      auto uv = UncertainValue2::UncertainValue2(massFrac);
      addElement(elm, uv);
   }

   double Composition::weightFraction(const Element::Element& elm, bool normalized) const
   {
      auto itr = mConstituents.find(&elm);
      if (itr != mConstituents.end()) {
         auto d = itr->second;
         return normalized ? normalize(d, mNormalization, true).doubleValue() : d.doubleValue();
      }
      return 0.0;
   }

   void Composition::recomputeStoiciometry()
   {
      mMoleNorm = UncertainValue2::ZERO();
      for (auto itr : mConstituents) {
         mMoleNorm = UncertainValue2::add(mMoleNorm, UncertainValue2::multiply(1.0 / itr.first->getAtomicWeight(), itr.second));
      }

      mConstituentsAtomic.clear();
      Element::UnorderedSetT constituentsAtomicKeys;
      for (auto ca : mConstituentsAtomic) {
         constituentsAtomicKeys.insert(ca.first);
      }

      for (auto elm : constituentsAtomicKeys) {
         auto uv0 = UncertainValue2::ZERO();
         auto uv1 = mConstituents.find(elm)->second;
         auto uv2 = UncertainValue2::multiply(elm->getAtomicWeight() / OUT_OF_THIS_MANY_ATOMS, mMoleNorm);
         auto uv3 = UncertainValue2::divide(uv1, uv2);
         auto moleFrac = (mMoleNorm.doubleValue() > 0.0 ? uv3 : uv0);

         mConstituentsAtomic.insert(std::make_pair(elm, moleFrac));
      }
      mOptimalRepresentation = Representation::WEIGHT_PCT;
   }

   void Composition::addElement(const Element::Element& elm, const UncertainValue2::UncertainValue2& massFrac)
   {
      mConstituents.insert(std::make_pair(&elm, massFrac));
      recomputeStoiciometry();
      renormalize();
   }

   UncertainValue2::UncertainValue2 Composition::weightFractionU(const Element::Element& elm, bool normalized) const
   {
      return weightFractionU(elm, normalized, true);
   }

   UncertainValue2::UncertainValue2 Composition::weightFractionU(const Element::Element& elm, bool normalized, bool positiveOnly) const
   {
      auto itr = mConstituents.find(&elm);
      if (itr != mConstituents.end()) {
         auto d = itr->second;
         return normalized ? normalize(d, mNormalization, positiveOnly) : d;
      }
      return UncertainValue2::ZERO();
   }

   Composition positiveDefinite(const Composition& comp)
   {
      Composition res;
      auto elemSet = comp.getElementSet();
      for (auto elm : elemSet) {
         if (comp.weightFraction(*elm, false) > 0.0) {
            res.addElement(*elm, comp.weightFractionU(*elm, false));
         }
      }
      return res;
   }

   void Composition::addElementByStoiciometry(const Element::Element& elm, const UncertainValue2::UncertainValue2& moleFrac)
   {
      mConstituentsAtomic.insert(std::make_pair(&elm, moleFrac));
      recomputeWeightFractions();
      renormalize();
   }

   void Composition::addElementByStoiciometry(const Element::Element& elm, double moleFrac)
   {
      addElementByStoiciometry(elm, UncertainValue2::UncertainValue2(moleFrac));
   }

   void Composition::clear()
   {
      mConstituents.clear();
      mConstituentsAtomic.clear();
      mNormalization = UncertainValue2::ONE();
      mAtomicNormalization = UncertainValue2::ONE();
      mMoleNorm = UncertainValue2::NaN();
   }

   void Composition::defineByWeightFraction(const Element::Element* elms[], int elmsLen, double wgtFracs[], int wgtFracsLen)
   {
      clear();
      if (elmsLen != wgtFracsLen) {
         printf("Composition::defineByWeightFraction1: elmsLen != wgtFracsLen (%d, %d)", elmsLen, wgtFracsLen);
      }
      for (int i = 0; i < elmsLen; ++i) {
         mConstituents.insert(std::make_pair(elms[i], UncertainValue2::UncertainValue2(wgtFracs[i])));
      }
      recomputeStoiciometry();
      renormalize();
   }

   void Composition::defineByWeightFraction(const Element::Element* elms[], int elmsLen, const UncertainValue2::UncertainValue2 wgtFracs[], int wgtFracsLen)
   {
      clear();
      if (elmsLen != wgtFracsLen) {
         printf("Composition::defineByWeightFraction2: elmsLen != wgtFracsLen (%d, %d)", elmsLen, wgtFracsLen);
      }
      for (int i = 0; i < elmsLen; ++i) {
         mConstituents.insert(std::make_pair(elms[i], wgtFracs[i]));
      }
      recomputeStoiciometry();
      renormalize();
   }

   void Composition::defineByMoleFraction(const Element::Element* elms[], int elmsLen, double moleFracs[], int moleFracsLen)
   {
      clear();
      if (elmsLen != moleFracsLen) {
         printf("Composition::defineByWeightFraction: elmsLen != moleFracsLen (%d, %d)", elmsLen, moleFracsLen);
      }
      double mfSum = 0;
      for (int k = 0; k < moleFracsLen; ++k) {
         mfSum += moleFracs[k];
      }
      for (int k = 0; k < moleFracsLen; ++k) {
         moleFracs[k] /= mfSum;
      }
      for (int i = 0; i < moleFracsLen; ++i) {
         auto elm = elms[i];
         mConstituentsAtomic.insert(std::make_pair(elm, UncertainValue2::UncertainValue2(moleFracs[i])));
      }
      recomputeWeightFractions();
      renormalize();
   }

   UncertainValue2::UncertainValue2 Composition::atomicPercentU(const Element::Element& elm) const
   {
      return atomicPercentU(elm, true);
   }

   UncertainValue2::UncertainValue2 Composition::atomicPercentU(const Element::Element& elm, bool positiveOnly) const
   {
      ConstituentsMapT::const_iterator itr = mConstituentsAtomic.find(&elm);
      return itr == mConstituentsAtomic.cend() ? UncertainValue2::ZERO() : normalize(itr->second, mAtomicNormalization, positiveOnly);
   }

   void Composition::recomputeWeightFractions()
   {
      UncertainValue2::UncertainValue2 totalWgt = UncertainValue2::ZERO();
      Element::UnorderedSetT constituentsAtomicKeys;
      for (auto ca : mConstituentsAtomic) {
         constituentsAtomicKeys.insert(ca.first);
      }
      for (auto elm : constituentsAtomicKeys) {
         auto uv0 = atomicPercentU(*elm);
         auto uv1 = UncertainValue2::multiply(elm->getAtomicWeight(), uv0);
         totalWgt = UncertainValue2::add(totalWgt, uv1);
      }

      mConstituents.clear();
      for (auto elm : constituentsAtomicKeys) {
         UncertainValue2::UncertainValue2 wgtFrac = UncertainValue2::multiply(elm->getAtomicWeight(), UncertainValue2::divide(atomicPercentU(*elm), totalWgt));
         mConstituents.insert(std::make_pair(elm, wgtFrac));
      }
      mOptimalRepresentation = Representation::STOICIOMETRY;
   }

   UncertainValue2::UncertainValue2 normalize(const UncertainValue2::UncertainValue2& val, const UncertainValue2::UncertainValue2& norm, bool positive)
   {
      UncertainValue2::UncertainValue2 uv;
      if (norm.doubleValue() > 0.0) {
         uv = UncertainValue2::divide(val, norm);
      }
      else {
         uv = val;
      }
      return positive ? UncertainValue2::positiveDefinite(uv) : uv;
   }

   void Composition::setOptimalRepresentation(const Representation opt)
   {
      switch (opt)
      {
      case UNDETERMINED:
         break;
      case WEIGHT_PCT:
         recomputeStoiciometry();
         break;
      case STOICIOMETRY:
         recomputeWeightFractions();
         break;
      }
   }

   Element::UnorderedSetT elementSet(const Composition compositions[], int len)
   {
      Element::UnorderedSetT elms;
      for (int i = 0; i < len; ++i) {
         auto elmset = compositions[i].getElementSet();
         for (auto elm : elmset) {
            elms.insert(elm);
         }
      }
      return elms;
   }

   void Composition::defineByMaterialFraction(const Composition compositions[], int compLen, double matFracs[], int matFracsLen)
   {
      if (compLen != matFracsLen) {
         printf("Composition::defineByMaterialFraction: lengths are different (%d, %d)", compLen , matFracsLen);
      }
      clear();
      auto elms = elementSet(compositions, compLen);
      int len = elms.size();
      std::vector<const Element::Element*> newElms(len);
      std::vector<UncertainValue2::UncertainValue2> frac(len);

      int ji = 0;
      for (auto el : elms) {
         UncertainValue2::UncertainValue2 sum = UncertainValue2::ZERO();
         for (int i = 0; i < compLen; ++i) {
            auto uv = UncertainValue2::multiply(matFracs[i], compositions[i].weightFractionU(*el, true));
            sum = UncertainValue2::add(sum, uv);
         }
         frac[ji] = sum;
         newElms[ji] = el;
         ++ji;
      }
      defineByWeightFraction(newElms.data(), len, frac.data(), len);
   }

   void Composition::removeElement(const Element::Element& el)
   {
      if (mConstituents.find(&el) != mConstituents.end()) {
         mConstituents.erase(&el);
         mConstituentsAtomic.erase(&el);
         // Don't recomputeStoiciometry or recomputeWeightFractions
         renormalize();
      }
   }

   bool Composition::containsElement(const Element::Element& el) const
   {
      return (mConstituents.find(&el) != mConstituents.end() && (mConstituents.find(&el)->second.doubleValue() > 0.0));
   }

   bool Composition::containsAll(const Element::UnorderedSetT& elms) const
   {
      for (auto elm : elms) {
         if (!containsElement(*elm)) {
            return false;
         }
      }
      return true;
   }

   double Composition::atomicPercent(const Element::Element& elm) const
   {
      return atomicPercentU(elm).doubleValue();
   }

   UncertainValue2::UncertainValue2 Composition::stoichiometryU(const Element::Element& elm) const
   {
      ConstituentsMapT::const_iterator itr = mConstituentsAtomic.find(&elm);
      if (itr != mConstituentsAtomic.cend()) {
         return itr->second;
      }
      return UncertainValue2::ZERO();
   }

   double Composition::stoichiometry(const Element::Element& elm) const
   {
      return stoichiometryU(elm).doubleValue();
   }

   double Composition::atomsPerKg(Element::Element& elm, bool normalized)
   {
      return weightFraction(elm, normalized) / elm.getMass();
   }

   UncertainValue2::UncertainValue2 Composition::atomsPerKgU(const Element::Element& elm, bool normalized) const
   {
      return UncertainValue2::multiply(1.0 / elm.getMass(), weightFractionU(elm, normalized));
   }

   UncertainValue2::UncertainValue2 Composition::weightAvgAtomicNumberU() const
   {
      UncertainValue2::UncertainValue2 res = UncertainValue2::ZERO();
      for (auto itr : mConstituents) {
         auto elm = itr.first;
         auto uv0 = weightFractionU(*elm, true);
         auto uv = UncertainValue2::multiply(elm->getAtomicNumber(), uv0);
         res = UncertainValue2::add(res, uv);
      }
      return res;
   }

   double Composition::weightAvgAtomicNumber() const
   {
      return weightAvgAtomicNumberU().doubleValue();
   }

   double Composition::sumWeightFraction() const
   {
      double sum = 0.0;
      for (auto itr : mConstituents) {
         auto uv = itr.second;
         if (uv.doubleValue() > 0.0) {
            sum += uv.doubleValue();
         }
      }
      return sum;
   }

   UncertainValue2::UncertainValue2 Composition::sumWeightFractionU() const
   {
      UncertainValue2::UncertainValue2 res = UncertainValue2::ZERO();
      for (auto itr : mConstituents) {
         auto val = itr.second;
         if (val.doubleValue() > 0.0) {
            res = UncertainValue2::add(res, val);
         }
      }
      return res;
   }

   Composition::CompositionNameT Composition::toString() const
   {
      //if (mName.size() == 0) {
      //   return descriptiveString(false);
      //}
      return mName;
   }

   //String::String Composition::stoichiometryString()
   //{
   //   final StringBuffer sb = new StringBuffer();
   //   final NumberFormat nf = new HalfUpFormat("0.####");
   //   for (final Element elm : getElementSet()) {
   //      final UncertainValue2 d0 = atomicPercentU(elm);
   //      if (sb.length() > 1)
   //         sb.append(",");
   //      sb.append(elm.toAbbrev());
   //      sb.append("(");
   //      sb.append(d0.format(nf));
   //      sb.append(" atoms)");
   //   }
   //   return sb.toString();
   //}

   //String::String weightPercentString(bool normalize)
   //{
   //   final StringBuffer sb = new StringBuffer();
   //   final NumberFormat nf = new HalfUpFormat("0.0000");
   //   for (final Element elm : getElementSet()) {
   //      final UncertainValue2 d0 = weightFractionU(elm, normalize);
   //      if (sb.length() > 1)
   //         sb.append(",");
   //      sb.append(elm.toAbbrev());
   //      sb.append("(");
   //      sb.append(d0.format(nf));
   //      sb.append(" mass frac)");
   //   }
   //   if (!normalize) {
   //      sb.append(",\u03A3=");
   //      sb.append(sumWeightPercentU().format(nf));
   //   }
   //   return sb.toString();
   //}

   //String::String descriptiveString(bool normalize)
   //{
   //   StringBuffer sb = new StringBuffer();
   //   if ((mName != null) && (mName.length() > 0))
   //      sb.append(mName + " = ");
   //   sb.append("[");
   //   if (mOptimalRepresentation == Representation.STOICIOMETRY)
   //      sb.append(stoichiometryString());
   //   else
   //      sb.append(weightPercentString(normalize));
   //   sb.append("]");
   //   return sb.toString();
   //}

   //Element::Element Composition::getNthElementByWeight(int n)
   //{
   //   LinkedListKV::Node<UncertainValue2::UncertainValue2, Element::Element>* tm = NULL;
   //   auto constituentsItr = mConstituents;
   //   while (constituentsItr != NULL) {
   //      auto el = constituentsItr->GetKey();
   //      UncertainValue2::UncertainValue2 wf = weightFractionU(el, true);
   //      // Add hoc mechanism to handle the case in which multiple elements are
   //      // present in the same weightPct.
   //      curandState state;
   //      curand_init(1, 0, 0, &state);
   //      while (LinkedListKV::ContainsKey(tm, wf, UncertainValue2::AreEqual)) {
   //         wf = UncertainValue2::add(1.0e-10 * curand_uniform(&state), wf);
   //      }
   //      LinkedListKV::InsertHead(&tm, wf, el);
   //      constituentsItr = constituentsItr->GetNext();
   //   }
   //   int j = 0;
   //   auto tmItr = tm;
   //   while (tmItr != NULL) {
   //      ++j;
   //      if (j == LinkedListKV::Size(mConstituents) - n) {
   //         return tmItr->GetValue();
   //      }
   //      tmItr = tmItr->GetNext();
   //   }
   //   return Element::None;
   //}

   //Element::Element Composition::getNthElementByAtomicFraction(int n)
   //{
   //   LinkedListKV::Node<UncertainValue2::UncertainValue2, Element::Element>* tm = NULL;
   //   auto constituentsItr = mConstituents;
   //   while (constituentsItr != NULL) {
   //      auto el = constituentsItr->GetKey();
   //      auto mf = atomicPercentU(el);
   //      curandState state;
   //      curand_init(1, 0, 0, &state);
   //      while (LinkedListKV::ContainsKey(tm, mf, UncertainValue2::AreEqual)) {
   //         mf = UncertainValue2::add(1.0e-10 * curand_uniform(&state), mf);
   //      }
   //      LinkedListKV::InsertHead(&tm, mf, el);
   //      constituentsItr = constituentsItr->GetNext();
   //   }
   //   int j = 0;
   //   auto tmItr = tm;
   //   while (tmItr != NULL) {
   //      ++j;
   //      if (j == n) {
   //         return tmItr->GetValue();
   //      }
   //      tmItr = tmItr->GetNext();
   //   }
   //   return Element::None;
   //}

   void Composition::setName(const CompositionNameT& name)
   {
      mName = name;
   }

   Composition::CompositionNameT Composition::getName() const
   {
      return mName;
   }

   int Composition::compareTo(const Composition& comp) const
   {
      if (this == &comp) {
         return 0;
      }
      // hashers have to be the same also
      auto i = mConstituents.begin();
      auto j = comp.mConstituents.begin();
      while (i != mConstituents.end() && j != comp.mConstituents.end()) {
         int zi = i->first->getAtomicNumber();
         int zj = j->first->getAtomicNumber();
         if (zi < zj) {
            return +1;
         }
         else if (zi > zj) {
            return -1;
         }
         else {
            UncertainValue2::UncertainValue2 ci = i->second;
            UncertainValue2::UncertainValue2 cj = j->second;
            if (ci.lessThan(cj)) {
               return -1;
            }
            else if (ci.greaterThan(cj)) {
               return +1;
            }
         }
         ++i;
         ++j;
      }
      if (!(i == mConstituents.end())) {
         return +1;
      }
      if (!(j == comp.mConstituents.end())) {
         return -1;
      }
      return 0;
   }

   Composition Composition::asComposition() const
   {
      Composition res;
      res.replicate(*this);
      return res;
   }

   Composition Composition::clone() const
   {
      return asComposition();
   }

   UncertainValue2::UncertainValue2 Composition::differenceU(const Composition& comp) const
   {
      // assert (comp.getElementCount() == this.getElementCount());
      UncertainValue2::UncertainValue2 delta = UncertainValue2::ZERO();
      Element::UnorderedSetT allElms;
      auto s0 = getElementSet();
      auto s1 = comp.getElementSet();
      allElms.insert(s0.begin(), s0.end());
      allElms.insert(s1.begin(), s1.end());

      for (auto el : allElms) {
         auto uv0 = comp.weightFractionU(*el, false);
         auto uv1 = weightFractionU(*el, false);
         auto uv2 = UncertainValue2::subtract(uv0, uv1);
         delta = UncertainValue2::add(delta, UncertainValue2::sqr(uv2));
      }
      return UncertainValue2::multiply(1.0 / allElms.size(), delta).sqrt();
   }

   double Composition::difference(const Composition& comp) const
   {
      return differenceU(comp).doubleValue();
   }

   Representation Composition::getOptimalRepresentation() const
   {
      return mOptimalRepresentation;
   }

   unsigned int Composition::hashCode() const
   {
      unsigned int result = 1;
      int PRIME = 31;
      //result = PRIME * result + mConstituents.hashCode();
      for (auto c : mConstituents) {
         result = PRIME * result + c.first->hashCode();
         result = PRIME * result + c.second.hashCode();
      }
      //result = PRIME * result + mConstituentsAtomic.hashCode();
      for (auto ca : mConstituentsAtomic) {
         result = PRIME * result + ca.first->hashCode();
         result = PRIME * result + ca.second.hashCode();
      }
      result = PRIME * result + std::hash<std::string>()(mName);
      long temp;
      temp = mNormalization.hashCode();
      result = PRIME * result + (int)(temp ^ (temp >> 32));
      result = PRIME * result + mOptimalRepresentation;
      if (result == INT_MAX)
         result = INT_MIN;
      return result;
   }

   Composition::ConstituentsMapT& Composition::GetConstituents()
   {
      return mConstituents;
   }

   bool Composition::sameConstituents(const ConstituentsMapT& constituents) const
   {
      for (auto c : mConstituents) {
         auto itr = constituents.find(c.first);
         if (itr == constituents.end()) {
            return false;
         }
         if (!(c.second == itr->second)) {
            return false;
         }
      }

      return true;
   }

   bool Composition::sameConstituentsAtomic(const Composition::ConstituentsMapT& constituentsAtomic) const
   {
      for (auto c : mConstituentsAtomic) {
         auto itr = constituentsAtomic.find(c.first);
         if (itr == constituentsAtomic.end()) return false;
         if (c.second == itr->second) return false;
      }

      return true;
   }

   bool Composition::equals(const Composition& obj) const
   {
      return this == &obj || *this == obj;
   }

   bool Composition::almostEquals(const Composition& other, double tol) const
   {
      if (this == &other) {
         return true;
      }
      //if (*((int*)&other) == NULL) {
      //   return false;
      //}
      if (abs(mNormalization.doubleValue() - other.mNormalization.doubleValue()) > tol) {
         return false;
      }
      Element::UnorderedSetT allElms;
      auto elms0 = other.getElementSet();
      for (auto e : elms0) {
         allElms.insert(e);
      }
      auto elms1 = getElementSet();
      for (auto e : elms1) {
         allElms.insert(e);
      }
      for (auto elm : allElms) {
         {
            UncertainValue2::UncertainValue2 uv1 = weightFractionU(*elm, false);
            UncertainValue2::UncertainValue2 uv2 = other.weightFractionU(*elm, false);
            if ((*((int*)&uv1) == NULL) || (*((int*)&uv2) == NULL)) {
               return false;
            }
            if ((abs(uv1.doubleValue() - uv2.doubleValue()) > tol)
               || (abs(uv1.uncertainty() - uv2.uncertainty()) > tol)) {
               return false;
            }
         }
         {
            UncertainValue2::UncertainValue2 uv1 = this->atomicPercentU(*elm);
            UncertainValue2::UncertainValue2 uv2 = other.atomicPercentU(*elm);
            //if ((*((int*)&uv1) == NULL) || (*((int*)&uv2) == NULL)) {
            //   return false;
            //}
            if ((abs(uv1.doubleValue() - uv2.doubleValue()) > tol)
               || (abs(uv1.uncertainty() - uv2.uncertainty()) > tol)) {
               return false;
            }
         }
      }
      return true;
   }

   Composition::ErrorMapT Composition::absoluteError(const Composition& std, bool normalize) const
   {
      Element::UnorderedSetT elms;
      auto elms0 = std.getElementSet();
      for (auto e : elms0) {
         elms.insert(e);
      }
      auto elms1 = getElementSet();
      for (auto e : elms1) {
         elms.insert(e);
      }
      ErrorMapT res;
      for (auto elm : elms) {
         double u = weightFractionU(*elm, normalize).doubleValue();
         double s = std.weightFractionU(*elm, normalize).doubleValue();
         res.insert(std::make_pair(elm, s != 0.0 ? (u - s) / s : (u == 0.0 ? 0.0 : 1.0)));
      }
      return res;
   }

   Composition::ErrorMapT Composition::relativeError(const Composition& std, bool normalize) const
   {
      Element::UnorderedSetT elms; auto elms0 = std.getElementSet();
      for (auto e : elms0) {
         elms.insert(e);
      }
      auto elms1 = getElementSet();
      for (auto e : elms1) {
         elms.insert(e);
      }
      ErrorMapT res;
      for (auto elm : elms) {
         double u = weightFractionU(*elm, normalize).doubleValue();
         double s = std.weightFractionU(*elm, normalize).doubleValue();
         res.insert(std::make_pair(elm, u - s));
      }
      return res;
   }

   bool Composition::isUncertain()
   {
      switch (mOptimalRepresentation) {
      case WEIGHT_PCT:
      case UNDETERMINED:
         for (auto i : mConstituents) {
            auto v = i.second;
            if (v.isUncertain()) {
               return true;
            }
         }
         break;
      case STOICIOMETRY:
         auto j = mConstituentsAtomic;
         for (auto j : mConstituentsAtomic) {
            auto v = j.second;
            if (v.isUncertain()) {
               return true;
            }
         }
         break;
      }
      return false;
   }

   UncertainValue2::UncertainValue2 Composition::meanAtomicNumberU() const
   {
      UncertainValue2::UncertainValue2 res = UncertainValue2::ZERO();
      auto elms = getElementSet();
      for (auto elm : elms) {
         auto uv = UncertainValue2::multiply(elm->getAtomicNumber(), atomicPercentU(*elm));
         res = UncertainValue2::add(res, uv);
      }
      return res;
   }

   double Composition::meanAtomicNumber() const
   {
      double res = 0.0;
      auto elms = getElementSet();
      for (auto elm : elms) {
         res += elm->getAtomicNumber() * atomicPercent(*elm);
      }
      return res;
   }

   void Composition::forceNormalization()
   {
      UncertainValue2::UncertainValue2 norm = sumWeightFractionU();
      ConstituentsMapT newConst;
      for (auto itr : mConstituents) {
         newConst.insert(std::make_pair(itr.first, norm.doubleValue() > 0.0 ? UncertainValue2::divide(itr.second, norm) : UncertainValue2::ZERO()));
      }
      mConstituents.clear();
      mConstituents = newConst;
      mOptimalRepresentation = Representation::WEIGHT_PCT;
      renormalize();
   }

   Composition parseGlass(char str[], int numlines)
   {
      Composition result;
      int pos = 0;
      int c = 0;
      for (int n = 0; n < numlines; ++n) {
         char line[256];
         int k = 0;
         while (str[c] != '\n') {
            line[k] = str[c];
            ++c;
            ++k;
         }
         ++c;
         Composition::CompositionNameT strline(line);
         if (pos == 0) {
            if (!strline.find("NBS GLASS K ")) {
               result.setName("K" "NBS GLASS K");
            }
            else if (!strline.find("CATIO")) {
               pos = 1;
            }
         }
         else if (pos == 1) {
            if (!strline.find("AVERAGE ATOMIC NUMBER")) {
               pos = 2;
            }
            else {
               char elmName[4];
               char elmWgtPct[16];
               int r = 0, l = 0;
               int tabCount = 0;
               while (true) {
                  if (line[r] == '\t') {
                     l = 0;
                     ++tabCount;
                  }
                  if (tabCount > 6) {
                     break;
                  }
                  if (tabCount == 0) {
                     elmName[l] = line[r];
                     ++l;
                  }
                  if (tabCount == 5) {
                     elmWgtPct[l] = line[r];
                     ++l;
                  }
                  ++r;
               }

               Element::Element elm = Element::byName(elmName);
               double wgtPct = std::atof(elmWgtPct);
               result.addElement(elm, wgtPct / 100.0);
            }
         }
         else if (pos == 2) {
            if (!strline.find("WEIGHT PERCENT OXYGEN")) {
               char oWgtPctStr[256];
               int r = 0, l = 0;
               int tabCount = 0;
               while (true) {
                  if (line[r] == '\t') {
                     ++tabCount;
                  }
                  if (tabCount > 1) {
                     break;
                  }
                  if (tabCount == 1) {
                     oWgtPctStr[l] = line[r];
                     ++l;
                  }
                  ++r;
               }

               double oWgtPct = std::atof(oWgtPctStr);
               result.addElement(Element::O, oWgtPct / 100.0);
               break;
            }
         }
      }
      return result;
   }

   Composition Composition::randomize(double offset, double proportional) const
   {
      srand(0);
      Composition res;
      auto elms = getElementSet();
      for (auto elm : elms) {
         double w = weightFraction(*elm, false);
         double v = w + w * (double)rand() / RAND_MAX * proportional + offset * (double)rand() / RAND_MAX;
         v = v > 0.0 ? v : 0.0;
         v = v < 1.1 ? v : 1.1;
         res.addElement(*elm, v);
      }
      return res;
   }

   int DIM = 9;
   long PROJECTORS[100]; // = createProjectors(2762689630628022905L);

   long mIndexHashS = INT_MAX;
   long mIndexHashL = INT_MAX;

   void createProjectors(long seed)
   {
      srand(NULL);

      std::unordered_set<int> eval;
      for (int j = 0; j < 100; ++j) {
         long tmp;
         do {
            long mult = 1;
            tmp = 0;
            for (int i = 0; i < DIM; ++i, mult *= 10) {
               double r = (double)rand() / (double)RAND_MAX;
               tmp += r * 2 * mult;
            }
         } while (eval.find(tmp) != eval.end());
         PROJECTORS[j] = tmp;
      }
   }

   long Composition::indexHashCodeS() const
   {
      if (mIndexHashS == INT_MAX) {
         long res = 0;
         auto elms = getElementSet();
         for (auto elm : elms) {
            int v = (int)sqrt(100.0 * weightFraction(*elm, false));
            v = v > 0 ? v : 0;
            v = v < 10 ? v : 0;
            res += v * PROJECTORS[elm->getAtomicNumber()];
         }
         mIndexHashS = res;
      }
      return mIndexHashS;
   }

   long Composition::indexHashCodeL() const
   {
      if (mIndexHashL == INT_MAX) {
         long res = 0;
         auto elms = getElementSet();
         for (auto elm : elms) {
            int v = (int)(10.0 * weightFraction(*elm, false));
            v = v > 0 ? v : 0;
            v = v < 10 ? v : 0;
            res += v * PROJECTORS[elm->getAtomicNumber()];
         }
         mIndexHashL = res;
      }
      return mIndexHashL;
   }
}
