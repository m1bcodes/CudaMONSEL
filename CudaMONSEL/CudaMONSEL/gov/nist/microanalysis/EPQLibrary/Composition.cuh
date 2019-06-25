#ifndef _COMPOSITION_CUH_
#define _COMPOSITION_CUH_

#include "Element.cuh"
#include "gov\nist\microanalysis\Utility\UncertainValue2.cuh"

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
      //typedef std::unordered_map<const Element::Element*, UncertainValue2::UncertainValue2, Element::HashFcn> ConstituentsMapT;
      typedef amp::unordered_map<const Element::Element*, UncertainValue2::UncertainValue2, Element::CompareFcn, UncertainValue2::CompareFcn, Element::HashFcn, UncertainValue2::HashFcn> ConstituentsMapT;
      //typedef std::unordered_map<const Element::Element*, double, Element::HashFcn> ErrorMapT;
      //typedef std::string StringT;
      typedef amp::string StringT;

      __host__ __device__ Composition();
      //Composition(const Composition& comp);
      Composition(const Element::Element* elms[], int elmsLen, const double massFracs[], int massFracsLen);
      Composition(const Element::Element& elm);
      Composition(const Element::Element* elms[], int elmsLen, const double massFracs[], int massFracsLen, char const * name);

      bool operator==(const Composition&) const;
      void operator=(const Composition&);

      Element::UnorderedSetT getElementSet() const;
      Element::OrderedSetT getSortedElements() const;
      __host__ __device__ int getElementCount() const;
      void addElement(int atomicNo, double massFrac);
      void addElement(int atomicNo, const UncertainValue2::UncertainValue2 massFrac);
      void addElement(const Element::Element& elm, double massFrac);
      double weightFraction(const Element::Element& elm, bool normalized) const;
      void addElement(const Element::Element& elm, const UncertainValue2::UncertainValue2& massFrac);
      UncertainValue2::UncertainValue2 weightFractionU(const Element::Element&, bool normalized) const;
      UncertainValue2::UncertainValue2 weightFractionU(const Element::Element&, bool normalized, bool positiveOnly) const;
      void addElementByStoiciometry(const Element::Element&, const UncertainValue2::UncertainValue2&);
      void addElementByStoiciometry(const Element::Element& elm, double moleFrac);
      UncertainValue2::UncertainValue2 atomicPercentU(const Element::Element&) const;
      UncertainValue2::UncertainValue2 atomicPercentU(const Element::Element&, bool positiveOnly) const;
      void defineByWeightFraction(const Element::Element* elms[], int elmsLen, double wgtFracs[], int wgtFracsLen);
      void defineByWeightFraction(const Element::Element* elms[], int elmsLen, const UncertainValue2::UncertainValue2 wgtFracs[], int wgtFracsLen);
      //void defineByWeightFraction(ConstituentsMapT map);
      UncertainValue2::UncertainValue2 stoichiometryU(const Element::Element&) const;
      double atomsPerKg(Element::Element& elm, bool normalized);
      UncertainValue2::UncertainValue2 atomsPerKgU(const Element::Element& elm, bool normalized) const;

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
      Element::Element GetElementBy<StringT>(StringT s)
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

      __host__ __device__ void defineByMoleFraction(const Element::Element* elms[], int elmsLen, const double moleFracs[], int moleFracsLen);
      void setOptimalRepresentation(const Representation opt);
      void defineByMaterialFraction(const Composition* compositions[], int compLen, const double matFracs[], int matFracsLen);
      void removeElement(const Element::Element&);
      bool containsElement(const Element::Element&) const;
      bool containsAll(const Element::UnorderedSetT& elms) const;
      double atomicPercent(const Element::Element& elm) const;
      double stoichiometry(const Element::Element& elm) const;
      UncertainValue2::UncertainValue2 weightAvgAtomicNumberU() const;
      double weightAvgAtomicNumber() const;
      double sumWeightFraction() const;
      UncertainValue2::UncertainValue2 sumWeightFractionU() const;
      const char* toString() const;
      //String::String stoichiometryString();
      //String::String weightPercentString(bool normalize);
      //String::String descriptiveString(bool normalize);
      //Element::Element getNthElementByWeight(int n);
      //Element::Element getNthElementByAtomicFraction(int n);
      __host__ __device__ void setName(const char* name);
      const char* getName() const;
      int compareTo(const Composition& comp) const;
      Composition asComposition() const;
      Composition clone() const;
      UncertainValue2::UncertainValue2 differenceU(const Composition& comp) const;
      double difference(const Composition& comp) const;
      Representation getOptimalRepresentation() const;
      unsigned int hashCode() const;
      bool sameConstituents(const ConstituentsMapT&) const;
      bool sameConstituentsAtomic(const ConstituentsMapT&) const;
      bool equals(const Composition& other) const;
      bool almostEquals(const Composition& other, double tol) const;
      //ErrorMapT absoluteError(const Composition& std, bool normalize) const;
      //ErrorMapT relativeError(const Composition& std, bool normalize) const;
      bool isUncertain();
      UncertainValue2::UncertainValue2 meanAtomicNumberU() const;
      double meanAtomicNumber() const;
      void forceNormalization();
      Composition randomize(double offset, double proportional) const;
      long indexHashCodeS() const;
      long indexHashCodeL() const;

      ConstituentsMapT& GetConstituents();

   private:
      Composition readResolve();
      void recomputeStoiciometry();
      __host__ __device__ void recomputeWeightFractions();

      ConstituentsMapT mConstituents;
      ConstituentsMapT mConstituentsAtomic;

      UncertainValue2::UncertainValue2 mNormalization;// = UncertainValue2::ONE();
      UncertainValue2::UncertainValue2 mAtomicNormalization;// = UncertainValue2::ONE();
      StringT mName;
      Representation mOptimalRepresentation;// = Representation::UNDETERMINED;
      UncertainValue2::UncertainValue2 mMoleNorm;// = UncertainValue2::NaN();

   protected:
      __host__ __device__ void renormalize();
      void replicate(const Composition& comp);
      void clear();
   };

   Composition positiveDefinite(const Composition& comp);
   UncertainValue2::UncertainValue2 normalize(const UncertainValue2::UncertainValue2& val, const UncertainValue2::UncertainValue2& norm, bool positive);
   Composition parseGlass(char str[], int numlines);
   void createProjectors(long seed);
}
#endif
