#ifndef _UNCERTAIN_VALUE_2_CUH_
#define _UNCERTAIN_VALUE_2_CUH_

#include <math_constants.h>

#include "string"
#include "unordered_map"
#include "unordered_set"
#include "vector"

namespace UncertainValue2
{
   typedef std::string UncertainValue2StringT;

   class Key
   {
   public:
      Key(UncertainValue2StringT src1, UncertainValue2StringT src2);
      bool operator==(const Key& k2) const;
      bool operator<(const Key& k2) const;
      size_t HashCode() const;

      //static bool AreEqual(Key& k1, Key& k2);

   private:
      UncertainValue2StringT mSource1, mSource2;
   };

   struct KeyCompareFcn
   {
      inline bool operator() (const Key& lhs, const Key& rhs)
      {
         return lhs < rhs;
      }
   };

   struct KeyHashFcn
   {
      inline size_t operator() (const Key& k) const
      {
         return k.HashCode();
      }
   };

   class Correlations
   {
      typedef std::unordered_map<Key, double, KeyHashFcn> CorrelationMap;
      typedef std::unordered_map<Key, double, KeyHashFcn>::iterator CorrelationMapItr;

   private:
      CorrelationMap mCorrelations;

   public:
      Correlations();

      void add(const UncertainValue2StringT& src1, const UncertainValue2StringT& src2, double corr);
      double get(const UncertainValue2StringT& src1, const UncertainValue2StringT& src2) const;
   };

   class UncertainValue2
   {
   public:
      typedef std::unordered_map<UncertainValue2StringT, double> ComponentMapT;
      typedef std::unordered_map<UncertainValue2StringT, double>::iterator ComponentMapTItr;
      typedef std::unordered_set<UncertainValue2StringT> KeySetT;
      typedef std::unordered_set<UncertainValue2StringT> KeySetTItr;
      typedef std::vector<UncertainValue2> ResultT;

      UncertainValue2();
      //~UncertainValue2();
      UncertainValue2(double v, char source[], double dv);
      UncertainValue2(double v);
      UncertainValue2(double v, double dv);
      UncertainValue2(double v, const ComponentMapT&);
      UncertainValue2(const UncertainValue2&);
      UncertainValue2& operator=(const UncertainValue2&);
      unsigned int hashCode();

      bool operator==(const UncertainValue2&) const;
      bool equals(UncertainValue2& uv);

      void assignInitialValue(double);
      void assignComponent(UncertainValue2StringT name, double sigma);
      double getComponent(const UncertainValue2StringT& src) const;
      ComponentMapT& getComponents();
      bool hasComponent(const UncertainValue2StringT& src) const;
      void renameComponent(const UncertainValue2StringT& oldName, const UncertainValue2StringT& newName);

      double doubleValue() const;
      bool isUncertain() const;
      double uncertainty() const;
      double variance() const;
      double fractionalUncertainty() const;

      int compareTo(UncertainValue2& o);
      bool lessThan(UncertainValue2& uv2);
      bool greaterThan(UncertainValue2& uv2);
      bool lessThanOrEqual(UncertainValue2& uv2);
      bool greaterThanOrEqual(UncertainValue2& uv2);

      UncertainValue2 sqrt();

      double variance(const Correlations& corr);
      double uncertainty(Correlations& corr);

      void PrintSigmas();

   private:
      ComponentMapT mSigmas;
      double mValue;
   };

   UncertainValue2 ONE();
   UncertainValue2 NaN();
   UncertainValue2 POSITIVE_INFINITY();
   UncertainValue2 NEGATIVE_INFINITY();
   UncertainValue2 ZERO();

   UncertainValue2 add(UncertainValue2 uvs[], int uvsLen);
   UncertainValue2 add(double a, UncertainValue2& uva, double b, UncertainValue2& uvb);
   UncertainValue2 subtract(UncertainValue2& uva, UncertainValue2& uvb);
   UncertainValue2 mean(UncertainValue2 uvs[], int uvsLen);
   UncertainValue2 weightedMean(UncertainValue2 uvs[], int uvsLen);
   UncertainValue2 uvmin(UncertainValue2 uvs[], int uvsLen);
   UncertainValue2 uvmax(UncertainValue2 uvs[], int uvsLen);
   UncertainValue2 add(UncertainValue2& v1, double v2);
   UncertainValue2 add(double v1, UncertainValue2& v2);
   UncertainValue2 add(UncertainValue2& v1, UncertainValue2& v2);
   UncertainValue2 multiply(double v1, UncertainValue2& v2);
   UncertainValue2 multiply(UncertainValue2& v1, UncertainValue2& v2);
   UncertainValue2 invert(UncertainValue2& v);
   UncertainValue2 divide(UncertainValue2& a, UncertainValue2& b);
   UncertainValue2 divide(double a, UncertainValue2& b);
   UncertainValue2 divide(UncertainValue2& a, double b);
   UncertainValue2 exp(UncertainValue2& x);
   UncertainValue2 log(UncertainValue2& v2);
   UncertainValue2 pow(UncertainValue2& v1, double n);
   UncertainValue2 sqrt(UncertainValue2& uv);
   UncertainValue2::ResultT quadratic(UncertainValue2& a, UncertainValue2& b, UncertainValue2& c);
   UncertainValue2 sqr(UncertainValue2& uv);
   UncertainValue2 negate(UncertainValue2& uv);
   UncertainValue2 atan(UncertainValue2& uv);
   UncertainValue2 atan2(UncertainValue2& y, UncertainValue2& x);
   UncertainValue2 positiveDefinite(const UncertainValue2& uv);
}

#endif