//#ifndef _UNCERTAIN_VALUE_2_CUH_
//#define _UNCERTAIN_VALUE_2_CUH_
//
//#include <cuda_runtime.h>
//#include <math_constants.h>
//
//#include "String.cuh"
//#include "Map.cuh"
////#include "..\..\..\..\Amphibian\String.cuh"
////#include "..\..\..\..\Amphibian\Map.cuh"
//
//namespace UncertainValue2
//{
//   class Key
//   {
//   public:
//      __device__ Key(String::String src1, String::String src2);
//      __device__ bool operator==(const Key& k2) const;
//      __device__ unsigned int Hashcode();
//
//      //__device__ static bool AreEqual(Key& k1, Key& k2);
//
//   private:
//      String::String mSource1, mSource2;
//   };
//
//   struct KeyCompareFcn
//   {
//      __device__ inline bool operator() (const Key& lhs, const Key& rhs) {
//         return lhs == rhs;
//      }
//   };
//
//   struct KeyHashFcn
//   {
//      __device__ inline unsigned int operator() (Key& k) {
//         return k.Hashcode();
//      }
//   };
//
//   class Correlations
//   {
//      typedef Map::Map<Key, double, KeyCompareFcn, Comparator::DoubleCompareFcn, KeyHashFcn, Hasher::DoubleHashFcn> CorrelationMap;
//      typedef Map::Iterator<Key, double, KeyCompareFcn, Comparator::DoubleCompareFcn, KeyHashFcn, Hasher::DoubleHashFcn> CorrelationMapItr;
//
//   private:
//      CorrelationMap mCorrelations;
//
//   public:
//      __device__ Correlations();
//
//      __device__ void add(String::String& src1, String::String& src2, double corr);
//      __device__ double get(String::String& src1, String::String& src2);
//   };
//
//   class UncertainValue2
//   {
//   public:
//      typedef Map::Map<String::String, double, String::CompareFcn, Comparator::DoubleCompareFcn, String::HashFcn, Hasher::DoubleHashFcn> ComponentMap;
//      typedef Map::Iterator<String::String, double, String::CompareFcn, Comparator::DoubleCompareFcn, String::HashFcn, Hasher::DoubleHashFcn> ComponentMapItr;
//      typedef Set::Set<String::String, String::CompareFcn, String::HashFcn> KeySet;
//      typedef Set::Iterator<String::String, String::CompareFcn, String::HashFcn> KeySetItr;
//
//      __device__ UncertainValue2();
//      //__device__ ~UncertainValue2();
//      __device__ UncertainValue2(double v, char source[], double dv);
//      __device__ UncertainValue2(double v);
//      __device__ UncertainValue2(double v, double dv);
//      __device__ UncertainValue2(double v, ComponentMap& sigmas);
//      __device__ UncertainValue2(const UncertainValue2&);
//      //__device__ UncertainValue2(UncertainValue2&);
//      __device__ UncertainValue2& operator=(UncertainValue2&);
//
//      __device__ bool operator==(UncertainValue2&);
//      __device__ bool equals(UncertainValue2& uv);
//
//      __device__ void assignInitialValue(double);
//      __device__ void assignComponent(String::String name, double sigma);
//      __device__ double getComponent(String::String src);
//      __device__ ComponentMap& getComponents();
//      __device__ bool hasComponent(String::String src);
//      __device__ void renameComponent(String::String oldName, String::String newName);
//
//      __device__ double doubleValue() const;
//      __device__ bool isUncertain();
//      __device__ double uncertainty();
//      __device__ double variance();
//      __device__ double fractionalUncertainty();
//      //__device__ bool operator==(UncertainValue2&);
//
//      __device__ int compareTo(UncertainValue2& o);
//      __device__ bool lessThan(UncertainValue2& uv2);
//      __device__ bool greaterThan(UncertainValue2& uv2);
//      __device__ bool lessThanOrEqual(UncertainValue2& uv2);
//      __device__ bool greaterThanOrEqual(UncertainValue2& uv2);
//
//      __device__ UncertainValue2 sqrt();
//
//      __device__ double variance(Correlations& corr);
//      __device__ double uncertainty(Correlations& corr);
//
//      __device__ void PrintSigmas();
//
//   private:
//      ComponentMap mSigmas;
//      double mValue;
//   };
//
//   //extern __device__ const char DEFAULT[8];
//   //extern __device__ int sDefIndex; // transient
//
//   __device__ UncertainValue2 ONE();
//   __device__ UncertainValue2 NaN();
//   __device__ UncertainValue2 POSITIVE_INFINITY();
//   __device__ UncertainValue2 NEGATIVE_INFINITY();
//   __device__ UncertainValue2 ZERO();
//
//   //extern __device__ const long long serialVersionUID;
//   extern __device__ const int MAX_LEN;
//
//   __device__ UncertainValue2 add(UncertainValue2 uvs[], int uvsLen);
//   __device__ UncertainValue2 add(double a, UncertainValue2& uva, double b, UncertainValue2& uvb);
//   __device__ UncertainValue2 subtract(UncertainValue2& uva, UncertainValue2& uvb);
//   __device__ UncertainValue2 mean(UncertainValue2 uvs[], int uvsLen);
//   __device__ UncertainValue2 weightedMean(UncertainValue2 uvs[], int uvsLen);
//   __device__ UncertainValue2 uvmin(UncertainValue2 uvs[], int uvsLen);
//   __device__ UncertainValue2 uvmax(UncertainValue2 uvs[], int uvsLen);
//   __device__ UncertainValue2 add(UncertainValue2& v1, double v2);
//   __device__ UncertainValue2 add(double v1, UncertainValue2& v2);
//   __device__ UncertainValue2 add(UncertainValue2& v1, UncertainValue2& v2);
//   __device__ UncertainValue2 multiply(double v1, UncertainValue2& v2);
//   __device__ UncertainValue2 multiply(UncertainValue2& v1, UncertainValue2& v2);
//   __device__ UncertainValue2 invert(UncertainValue2& v);
//   __device__ UncertainValue2 divide(UncertainValue2& a, UncertainValue2& b);
//   __device__ UncertainValue2 divide(double a, UncertainValue2& b);
//   __device__ UncertainValue2 divide(UncertainValue2& a, double b);
//   __device__ UncertainValue2 exp(UncertainValue2& x);
//   __device__ UncertainValue2 log(UncertainValue2& v2);
//   __device__ UncertainValue2 pow(UncertainValue2& v1, double n);
//   __device__ UncertainValue2 sqrt(UncertainValue2& uv);
//   __device__ LinkedList::Node<UncertainValue2>* quadratic(UncertainValue2& a, UncertainValue2& b, UncertainValue2& c);
//   __device__ UncertainValue2 sqr(UncertainValue2& uv);
//   __device__ UncertainValue2 negate(UncertainValue2& uv);
//   __device__ UncertainValue2 atan(UncertainValue2& uv);
//   __device__ UncertainValue2 atan2(UncertainValue2& y, UncertainValue2& x);
//   __device__ UncertainValue2 positiveDefinite(UncertainValue2& uv);
//
//   //__device__ void InitializeSpecialUncertainValues();
//   __device__ bool AreEqual(UncertainValue2&, UncertainValue2&);
//}
//
//#endif