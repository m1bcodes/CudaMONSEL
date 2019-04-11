#ifndef _COMPOSITION_CUH_
#define _COMPOSITION_CUH_

#include "Element.cuh"
#include "..\Utility\UncertainValue2.cuh"

#include <map>
#include <set>

namespace Composition
{
   enum Representation
   {
      UNDETERMINED, WEIGHT_PCT, STOICIOMETRY
   };

   class Composition
   {
   public:
      //typedef std::unordered_map<Element::Element, double, Element::HashFcn> CompositionMap;
      //typedef std::unordered_map<Element::Element, double, Element::HashFcn>::iterator CompositionMapItr;
      typedef std::set<Element::Element> OrderedElementSetT;
      typedef std::set<Element::Element>::iterator OrderedElementSetTItr;
      typedef std::unordered_set<Element::Element, Element::HashFcn> ElementSetT;
      typedef std::unordered_set<Element::Element, Element::HashFcn>::iterator ElementSetTItr;
      typedef std::unordered_map<Element::Element, UncertainValue2::UncertainValue2, Element::HashFcn> ConstituentsMapT;
      typedef std::unordered_map<Element::Element, UncertainValue2::UncertainValue2, Element::HashFcn>::iterator ConstituentsMapTItr;
      typedef std::unordered_map<Element::Element, double, Element::HashFcn> ErrorMapT;
      typedef std::string CompositionNameT;

      Composition();
      ~Composition();
      Composition(const Composition& comp);
      Composition(Element::Element elms[], int elmsLen, double massFracs[], int massFracsLen);
      Composition(Element::Element elm);
      Composition(Element::Element elms[], int elmsLen, double massFracs[], int massFracsLen, char const * name);

      bool operator==(const Composition&) const;
      void operator=(const Composition&);

      ElementSetT getElementSet() const;
      OrderedElementSetT getSortedElements();
      int getElementCount();
      void addElement(int atomicNo, double massFrac);
      void addElement(int atomicNo, UncertainValue2::UncertainValue2 massFrac);
      void addElement(const Element::Element& elm, double massFrac);
      double weightFraction(const Element::Element& elm, bool normalized);
      void addElement(const Element::Element& elm, const UncertainValue2::UncertainValue2& massFrac);
      UncertainValue2::UncertainValue2 weightFractionU(const Element::Element&, bool normalized);
      UncertainValue2::UncertainValue2 weightFractionU(const Element::Element&, bool normalized, bool positiveOnly);
      void addElementByStoiciometry(const Element::Element&, const UncertainValue2::UncertainValue2&);
      void addElementByStoiciometry(Element::Element elm, double moleFrac);
      UncertainValue2::UncertainValue2 atomicPercentU(const Element::Element&);
      UncertainValue2::UncertainValue2 atomicPercentU(const Element::Element&, bool positiveOnly);
      void defineByWeightFraction(Element::Element elms[], int elmsLen, double wgtFracs[], int wgtFracsLen);
      void defineByWeightFraction(Element::Element elms[], int elmsLen, UncertainValue2::UncertainValue2 wgtFracs[], int wgtFracsLen);
      void defineByWeightFraction(ConstituentsMapT map);
      UncertainValue2::UncertainValue2 stoichiometryU(const Element::Element&);
      double atomsPerKg(Element::Element& elm, bool normalized);
      UncertainValue2::UncertainValue2 atomsPerKgU(Element::Element elm, bool normalized);

      template<typename T>
      Element::Element GetElementBy(T k)
      {
         return NULL;
      }

      template<>
      Element::Element GetElementBy<int>(int k)
      {
         return Element::byAtomicNumber(k);
      }

      template<>
      Element::Element GetElementBy<CompositionNameT>(CompositionNameT s)
      {
         return Element::byName(s.c_str());
      }

      template<>
      Element::Element GetElementBy<Element::Element>(Element::Element e)
      {
         return e;
      }

      //template<typename T>
      //void defineByWeightFraction(LinkedListKV::Node<T, double>* map)
      //{
      //   while (map != NULL) {
      //      double wp = map->GetValue();
      //      T key = map->GetKey();
      //      Element::Element elm = GetElementBy<T>(key);
      //      if (*((int*)&elm) == NULL) {
      //         printf("Composition::defineByWeightFraction: wrong type");
      //      }
      //      if (elm.getAtomicNumber() == Element::elmNone) {
      //         printf("Composition::defineByWeightFraction: bad element");
      //      }
      //      if ((*((int*)&elm) != NULL) && (elm.getAtomicNumber() == Element::elmNone)) {
      //         LinkedListKV::InsertHead(&mConstituents, elm, UncertainValue2::UncertainValue2(wp));
      //      }
      //   }
      //   recomputeStoiciometry();
      //   renormalize();
      //}

      void defineByMoleFraction(Element::Element elms[], int elmsLen, double moleFracs[], int moleFracsLen);
      void setOptimalRepresentation(Representation opt);
      void defineByMaterialFraction(Composition compositions[], int compLen, double matFracs[], int matFracsLen);
      void removeElement(const Element::Element&);
      bool containsElement(const Element::Element&);
      bool containsAll(const ElementSetT& elms);
      double atomicPercent(Element::Element elm);
      double stoichiometry(Element::Element elm);
      UncertainValue2::UncertainValue2 weightAvgAtomicNumberU();
      double weightAvgAtomicNumber();
      double sumWeightFraction();
      UncertainValue2::UncertainValue2 sumWeightFractionU();
      ////String::String toString();
      ////String::String stoichiometryString();
      ////String::String weightPercentString(bool normalize);
      ////String::String descriptiveString(bool normalize);
      //Element::Element getNthElementByWeight(int n);
      //Element::Element getNthElementByAtomicFraction(int n);
      void setName(CompositionNameT name);
      CompositionNameT getName();
      int compareTo(Composition& comp);
      Composition asComposition();
      Composition clone();
      UncertainValue2::UncertainValue2 differenceU(Composition& comp);
      double difference(Composition& comp);
      Representation getOptimalRepresentation();
      unsigned int hashCode();
      bool sameConstituents(const ConstituentsMapT&) const;
      bool sameConstituentsAtomic(const ConstituentsMapT&) const;
      bool equals(const Composition& other) const;
      bool almostEquals(Composition& other, double tol);
      ErrorMapT absoluteError(Composition& std, bool normalize);
      ErrorMapT relativeError(Composition& std, bool normalize);
      bool isUncertain();
      UncertainValue2::UncertainValue2 meanAtomicNumberU();
      double meanAtomicNumber();
      void forceNormalization();
      Composition randomize(double offset, double proportional);
      long indexHashCodeS();
      long indexHashCodeL();

      ConstituentsMapT& GetConstituents();

   private:
      Composition readResolve();
      void recomputeStoiciometry();
      void recomputeWeightFractions();

      ConstituentsMapT mConstituents;
      ConstituentsMapT mConstituentsAtomic;

      UncertainValue2::UncertainValue2 mNormalization;// = UncertainValue2::ONE();
      UncertainValue2::UncertainValue2 mAtomicNormalization;// = UncertainValue2::ONE();
      CompositionNameT mName;
      Representation mOptimalRepresentation;// = Representation::UNDETERMINED;
      UncertainValue2::UncertainValue2 mMoleNorm;// = UncertainValue2::NaN();

   protected:
      void renormalize();
      void replicate(const Composition& comp);
      void clear();
   };

   Composition positiveDefinite(Composition& comp);
   UncertainValue2::UncertainValue2 normalize(UncertainValue2::UncertainValue2& val, UncertainValue2::UncertainValue2& norm, bool positive);
   Composition parseGlass(char str[], int numlines);
   void createProjectors(long seed);
}
#endif
