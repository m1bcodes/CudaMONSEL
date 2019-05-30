#ifndef _UNCERTAIN_VALUE_2_CUH_
#define _UNCERTAIN_VALUE_2_CUH_

#include <math_constants.h>

#include "string"
#include "unordered_map"
#include "unordered_set"
#include "vector"

namespace UncertainValue2
{
   typedef std::string StringT;

   class Key
   {
   public:
      Key(const StringT& src1, const StringT& src2);
      bool operator==(const Key& k2) const;
      bool operator<(const Key& k2) const;
      size_t HashCode() const;

   private:
      StringT mSource1, mSource2;
   };

   struct KeyCompareFcn
   {
      inline bool operator() (const Key& lhs, const Key& rhs) const
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

   private:
      CorrelationMap mCorrelations;

   public:
      Correlations();

      void add(const StringT& src1, const StringT& src2, double corr);
      double get(const StringT& src1, const StringT& src2) const;
   };

   class UncertainValue2
   {
   public:
      typedef std::unordered_map<StringT, double> ComponentMapT;
      typedef std::unordered_map<StringT, double>::iterator ComponentMapTItr;
      typedef std::unordered_set<StringT> KeySetT;
      typedef std::unordered_set<StringT> KeySetTItr;
      typedef std::vector<UncertainValue2> ResultT;

      UncertainValue2();
      //~UncertainValue2();
      UncertainValue2(double v, const char source[], double dv);
      UncertainValue2(double v);
      UncertainValue2(double v, double dv);
      UncertainValue2(double v, const ComponentMapT&);
      UncertainValue2(const UncertainValue2&);
      UncertainValue2& operator=(const UncertainValue2&);
      unsigned int hashCode() const;

      bool operator==(const UncertainValue2&) const;
      bool equals(const UncertainValue2& uv) const;

      void assignComponent(const StringT&, double sigma);
      double getComponent(const StringT&) const;
      ComponentMapT& getComponents();
      const ComponentMapT& getComponentsConst() const;
      ComponentMapT::const_iterator getComponentsItrBegin() const;
      ComponentMapT::const_iterator getComponentsItrEnd() const;
      bool hasComponent(const StringT&) const;
      void renameComponent(const StringT& oldName, const StringT& newName);

      double doubleValue() const;
      bool isUncertain() const;
      double uncertainty() const;
      double variance() const;
      double fractionalUncertainty() const;

      int compareTo(const UncertainValue2& o);
      bool lessThan(const UncertainValue2& uv2);
      bool greaterThan(const UncertainValue2& uv2);
      bool lessThanOrEqual(const UncertainValue2& uv2);
      bool greaterThanOrEqual(const UncertainValue2& uv2);

      UncertainValue2 sqrt() const;

      double variance(const Correlations& corr);
      double uncertainty(Correlations& corr);

   private:
      ComponentMapT mSigmas;
      double mValue;
   };

   UncertainValue2 ONE();
   UncertainValue2 NaN();
   UncertainValue2 POSITIVE_INFINITY();
   UncertainValue2 NEGATIVE_INFINITY();
   UncertainValue2 ZERO();

   UncertainValue2 add(const UncertainValue2 uvs[], int uvsLen);
   UncertainValue2 add(double a, const UncertainValue2& uva, double b, const UncertainValue2& uvb);
   UncertainValue2 subtract(const UncertainValue2& uva, const UncertainValue2& uvb);
   UncertainValue2 mean(const UncertainValue2 uvs[], int uvsLen);
   UncertainValue2 weightedMean(const UncertainValue2 uvs[], int uvsLen);
   UncertainValue2 uvmin(const UncertainValue2 uvs[], int uvsLen);
   UncertainValue2 uvmax(const UncertainValue2 uvs[], int uvsLen);
   UncertainValue2 add(const UncertainValue2& v1, double v2);
   UncertainValue2 add(double v1, const UncertainValue2& v2);
   UncertainValue2 add(const UncertainValue2& v1, const UncertainValue2& v2);
   UncertainValue2 multiply(double v1, const UncertainValue2& v2);
   UncertainValue2 multiply(const UncertainValue2& v1, const UncertainValue2& v2);
   UncertainValue2 invert(const UncertainValue2& v);
   UncertainValue2 divide(const UncertainValue2& a, const UncertainValue2& b);
   void divide(const UncertainValue2& a, const UncertainValue2& b, UncertainValue2&);
   UncertainValue2 divide(double a, const UncertainValue2& b);
   UncertainValue2 divide(const UncertainValue2& a, double b);
   UncertainValue2 exp(const UncertainValue2& x);
   UncertainValue2 log(const UncertainValue2& v2);
   UncertainValue2 pow(const UncertainValue2& v1, double n);
   UncertainValue2 sqrt(const UncertainValue2& uv);
   UncertainValue2::ResultT quadratic(const UncertainValue2& a, const UncertainValue2& b, const UncertainValue2& c);
   UncertainValue2 sqr(const UncertainValue2& uv);
   UncertainValue2 negate(const UncertainValue2& uv);
   UncertainValue2 atan(const UncertainValue2& uv);
   UncertainValue2 atan2(const UncertainValue2& y, const UncertainValue2& x);
   UncertainValue2 positiveDefinite(const UncertainValue2& uv);
}

#endif