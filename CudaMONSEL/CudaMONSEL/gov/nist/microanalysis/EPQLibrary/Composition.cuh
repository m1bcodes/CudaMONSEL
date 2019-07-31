#ifndef _COMPOSITION_CUH_
#define _COMPOSITION_CUH_

#include "Element.cuh"
#include "gov\nist\microanalysis\Utility\UncertainValue2.cuh"

#include <map>
#include <set>

namespace Composition
{
   typedef UncertainValue2::data_type data_type;

   enum Representation
   {
      UNDETERMINED, WEIGHT_PCT, STOICIOMETRY
   };

   class Composition
   {
   public:
      //typedef std::unordered_map<const Element::Element*, UncertainValue2::UncertainValue2, Element::HashFcn> ConstituentsMapT;
      typedef amp::unordered_map<const Element::Element*, UncertainValue2::UncertainValue2, Element::CompareFcn, UncertainValue2::CompareFcn, Element::HashFcn, UncertainValue2::HashFcn> ConstituentsMapT;
      //typedef std::unordered_map<const Element::Element*, data_type, Element::HashFcn> ErrorMapT;
      //typedef std::string StringT;
      typedef amp::string StringT;

      __host__ __device__ Composition();
      __host__ __device__ Composition(const Composition& comp);
      Composition(const Element::Element* elms[], int elmsLen, const data_type massFracs[], int massFracsLen);
      Composition(const Element::Element& elm);
      __host__ __device__ Composition(const Element::Element* elms[], int elmsLen, const data_type massFracs[], int massFracsLen, char const * name);

      bool operator==(const Composition&) const;
      void operator=(const Composition&);

      __host__ __device__ Element::UnorderedSetT getElementSet() const;
      Element::OrderedSetT getSortedElements() const;
      __host__ __device__ int getElementCount() const;
      void addElement(int atomicNo, data_type massFrac);
      void addElement(int atomicNo, const UncertainValue2::UncertainValue2 massFrac);
      void addElement(const Element::Element& elm, data_type massFrac);
      __host__ __device__ data_type weightFraction(const Element::Element& elm, const bool normalized) const;
      void addElement(const Element::Element& elm, const UncertainValue2::UncertainValue2& massFrac);
      UncertainValue2::UncertainValue2 weightFractionU(const Element::Element&, bool normalized) const;
      UncertainValue2::UncertainValue2 weightFractionU(const Element::Element&, bool normalized, bool positiveOnly) const;
      void addElementByStoiciometry(const Element::Element&, const UncertainValue2::UncertainValue2&);
      void addElementByStoiciometry(const Element::Element& elm, data_type moleFrac);
      __host__ __device__ UncertainValue2::UncertainValue2 atomicPercentU(const Element::Element&) const;
      __host__ __device__ UncertainValue2::UncertainValue2 atomicPercentU(const Element::Element&, const bool positiveOnly) const;
      void defineByWeightFraction(const Element::Element* elms[], int elmsLen, data_type wgtFracs[], int wgtFracsLen);
      void defineByWeightFraction(const Element::Element* elms[], int elmsLen, const UncertainValue2::UncertainValue2 wgtFracs[], int wgtFracsLen);
      //void defineByWeightFraction(ConstituentsMapT map);
      UncertainValue2::UncertainValue2 stoichiometryU(const Element::Element&) const;
      data_type atomsPerKg(Element::Element& elm, bool normalized);
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
      //void defineByWeightFraction(LinkedListKV::Node<T, data_type>* map)
      //{
      //   while (map != NULL) {
      //      data_type wp = map->GetValue();
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

      __host__ __device__ void defineByMoleFraction(const Element::Element* elms[], int elmsLen, const data_type moleFracs[], int moleFracsLen);
      void setOptimalRepresentation(const Representation opt);
      void defineByMaterialFraction(const Composition* compositions[], int compLen, const data_type matFracs[], int matFracsLen);
      void removeElement(const Element::Element&);
      bool containsElement(const Element::Element&) const;
      bool containsAll(const Element::UnorderedSetT& elms) const;
      data_type atomicPercent(const Element::Element& elm) const;
      data_type stoichiometry(const Element::Element& elm) const;
      UncertainValue2::UncertainValue2 weightAvgAtomicNumberU() const;
      data_type weightAvgAtomicNumber() const;
      data_type sumWeightFraction() const;
      UncertainValue2::UncertainValue2 sumWeightFractionU() const;
      __host__ __device__ const char* toString() const;
      //String::String stoichiometryString();
      //String::String weightPercentString(bool normalize);
      //String::String descriptiveString(bool normalize);
      //Element::Element getNthElementByWeight(int n);
      //Element::Element getNthElementByAtomicFraction(int n);
      __host__ __device__ void setName(const char* name);
      __host__ __device__ const char* getName() const;
      int compareTo(const Composition& comp) const;
      Composition asComposition() const;
      Composition clone() const;
      UncertainValue2::UncertainValue2 differenceU(const Composition& comp) const;
      data_type difference(const Composition& comp) const;
      Representation getOptimalRepresentation() const;
      unsigned int hashCode() const;
      bool sameConstituents(const ConstituentsMapT&) const;
      bool sameConstituentsAtomic(const ConstituentsMapT&) const;
      bool equals(const Composition& other) const;
      bool almostEquals(const Composition& other, data_type tol) const;
      //ErrorMapT absoluteError(const Composition& std, bool normalize) const;
      //ErrorMapT relativeError(const Composition& std, bool normalize) const;
      bool isUncertain();
      UncertainValue2::UncertainValue2 meanAtomicNumberU() const;
      data_type meanAtomicNumber() const;
      void forceNormalization();
      Composition randomize(data_type offset, data_type proportional) const;
      long indexHashCodeS() const;
      long indexHashCodeL() const;

      ConstituentsMapT& GetConstituents();

      __host__ __device__ void renormalize();

   private:
      Composition readResolve();
      __host__ __device__ void recomputeStoiciometry();
      __host__ __device__ void recomputeWeightFractions();

      ConstituentsMapT mConstituents;
      ConstituentsMapT mConstituentsAtomic;

      UncertainValue2::UncertainValue2 mNormalization;// = UncertainValue2::ONE();
      UncertainValue2::UncertainValue2 mAtomicNormalization;// = UncertainValue2::ONE();
      StringT mName;
      Representation mOptimalRepresentation;// = Representation::UNDETERMINED;
      UncertainValue2::UncertainValue2 mMoleNorm;// = UncertainValue2::NaN();

   protected:
      //__host__ __device__ void renormalize();
      __host__ __device__ void replicate(const Composition& comp);
      __host__ __device__ void clear();
   };

   Composition positiveDefinite(const Composition& comp);
   __host__ __device__ UncertainValue2::UncertainValue2 normalize(const UncertainValue2::UncertainValue2& val, const UncertainValue2::UncertainValue2& norm, const bool positive);
   Composition parseGlass(char str[], int numlines);
   void createProjectors(long seed);
}
#endif
