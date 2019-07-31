#ifndef _UNCERTAIN_VALUE_2_CUH_
#define _UNCERTAIN_VALUE_2_CUH_

#include <math_constants.h>

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Amphibian/String.cuh"
#include "Amphibian/unordered_map.cuh"
#include "Amphibian/unordered_set.cuh"

namespace UncertainValue2
{
   typedef amp::string StringT;
   typedef float data_type;

   class Key
   {
   public:
      Key(const StringT& src1, const StringT& src2);
      bool operator==(const Key& k2) const;
      bool operator<(const Key& k2) const;
      size_t hashCode() const;

   private:
      StringT mSource1, mSource2;
   };

   struct KeyCompareFcn
   {
      inline bool operator() (const Key& lhs, const Key& rhs) const
      {
         return lhs == rhs;
      }
   };

   struct KeyHashFcn
   {
      inline size_t operator() (const Key& k) const
      {
         return k.hashCode();
      }
   };

   class Correlations
   {
      typedef std::unordered_map<Key, data_type, KeyHashFcn> CorrelationMapT;

   private:
      CorrelationMapT mCorrelations;

   public:
      Correlations();

      void add(const StringT& src1, const StringT& src2, data_type corr);
      data_type get(const StringT& src1, const StringT& src2) const;
   };

   class UncertainValue2
   {
   public:
      //typedef std::unordered_map<StringT, data_type, amp::string_hash> ComponentMapT;
      typedef amp::unordered_map<StringT, data_type, amp::string_cmp, Comparator::DoubleCompareFcn, amp::string_hash, Hasher::DoubleHashFcn> ComponentMapT;
      //typedef std::unordered_set<StringT, amp::string_hash> KeySetT;
      typedef amp::unordered_set<StringT, amp::string_hash, amp::string_cmp> KeySetT;
      //typedef std::vector<UncertainValue2> ResultT;

      __host__ __device__ UncertainValue2();
      //~UncertainValue2();
      UncertainValue2(data_type v, const char source[], data_type dv);
      __host__ __device__ UncertainValue2(const data_type v);
      __host__ __device__ UncertainValue2(const data_type v, const data_type dv);
      __host__ __device__ UncertainValue2(data_type v, const ComponentMapT&);
      __host__ __device__ UncertainValue2(const UncertainValue2&);
      __host__ __device__ UncertainValue2& operator=(const UncertainValue2&);
      __host__ __device__ unsigned int hashCode() const;

      bool operator==(const UncertainValue2&) const;
      bool equals(const UncertainValue2& uv) const;

      __host__ __device__ void assignComponent(const StringT&, const data_type sigma);
      __host__ __device__ data_type getComponent(const StringT&) const;
      ComponentMapT& getComponents();
      __host__ __device__ const ComponentMapT& getComponents() const;
      __host__ __device__ ComponentMapT::const_iterator getComponentsItrBegin() const;
      __host__ __device__ ComponentMapT::const_iterator getComponentsItrEnd() const;
      bool hasComponent(const StringT&) const;
      void renameComponent(const StringT& oldName, const StringT& newName);

      __host__ __device__ data_type doubleValue() const;
      bool isUncertain() const;
      __host__ __device__ data_type uncertainty() const;
      __host__ __device__ data_type variance() const;
      data_type fractionalUncertainty() const;

      int compareTo(const UncertainValue2& o);
      bool lessThan(const UncertainValue2& uv2);
      bool greaterThan(const UncertainValue2& uv2);
      bool lessThanOrEqual(const UncertainValue2& uv2);
      bool greaterThanOrEqual(const UncertainValue2& uv2);

      UncertainValue2 sqrt() const;

      data_type variance(const Correlations& corr);
      data_type uncertainty(Correlations& corr);

   private:
      ComponentMapT mSigmas;
      data_type mValue;
   };

   __host__ __device__ UncertainValue2 ONE();
   __host__ __device__ UncertainValue2 NaN();
   UncertainValue2 POSITIVE_INFINITY();
   UncertainValue2 NEGATIVE_INFINITY();
   __host__ __device__ UncertainValue2 ZERO();

   UncertainValue2 add(const UncertainValue2 uvs[], int uvsLen);
   __host__ __device__ UncertainValue2 add(const data_type a, const UncertainValue2& uva, const data_type b, const UncertainValue2& uvb);
   UncertainValue2 subtract(const UncertainValue2& uva, const UncertainValue2& uvb);
   UncertainValue2 mean(const UncertainValue2 uvs[], int uvsLen);
   UncertainValue2 weightedMean(const UncertainValue2 uvs[], int uvsLen);
   UncertainValue2 uvmin(const UncertainValue2 uvs[], int uvsLen);
   UncertainValue2 uvmax(const UncertainValue2 uvs[], int uvsLen);
   UncertainValue2 add(const UncertainValue2& v1, data_type v2);
   UncertainValue2 add(data_type v1, const UncertainValue2& v2);
   __host__ __device__ UncertainValue2 add(const UncertainValue2& v1, const UncertainValue2& v2);
   __host__ __device__ UncertainValue2 multiply(data_type v1, const UncertainValue2& v2);
   UncertainValue2 multiply(const UncertainValue2& v1, const UncertainValue2& v2);
   UncertainValue2 invert(const UncertainValue2& v);
   __host__ __device__ UncertainValue2 divide(const UncertainValue2& a, const UncertainValue2& b);
   __host__ __device__ void divide(const UncertainValue2& a, const UncertainValue2& b, UncertainValue2&);
   UncertainValue2 divide(data_type a, const UncertainValue2& b);
   UncertainValue2 divide(const UncertainValue2& a, data_type b);
   UncertainValue2 exp(const UncertainValue2& x);
   UncertainValue2 log(const UncertainValue2& v2);
   UncertainValue2 pow(const UncertainValue2& v1, data_type n);
   UncertainValue2 sqrt(const UncertainValue2& uv);
   //UncertainValue2::ResultT quadratic(const UncertainValue2& a, const UncertainValue2& b, const UncertainValue2& c);
   UncertainValue2 sqr(const UncertainValue2& uv);
   UncertainValue2 negate(const UncertainValue2& uv);
   UncertainValue2 atan(const UncertainValue2& uv);
   UncertainValue2 atan2(const UncertainValue2& y, const UncertainValue2& x);
   __host__ __device__ UncertainValue2 positiveDefinite(const UncertainValue2& uv);

   struct CompareFcn
   {
      __host__ __device__ inline bool operator() (const UncertainValue2& uv0, const UncertainValue2& uv1) const
      {
         return uv0 == uv1;
      }
   };

   struct HashFcn
   {
      __host__ __device__ inline unsigned int operator() (const UncertainValue2& uv) const
      {
         return uv.hashCode();
      }
   };
}

#endif