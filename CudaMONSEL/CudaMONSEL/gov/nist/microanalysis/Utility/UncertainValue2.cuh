#ifndef _UNCERTAIN_VALUE_2_CUH_
#define _UNCERTAIN_VALUE_2_CUH_

#include <cuda_runtime.h>
#include <math_constants.h>

#include "..\..\..\..\Amphibian\String.cuh"
#include "..\..\..\..\Amphibian\Map.cuh"

namespace UncertainValue2
{
   class Key
   {
   private:
      String::String mSource1, mSource2;

   public:
      __device__ Key(String::String src1, String::String src2);
      __device__ bool operator==(Key& k2);

      __device__ static bool AreEqual(Key& k1, Key& k2);
   };

   class Correlations
   {
   private:
      Map::Map<Key, double> mCorrelations;

   public:
      __device__ Correlations();

      __device__ void add(String::String& src1, String::String& src2, double corr);
      __device__ double get(String::String& src1, String::String& src2);
   };

   class UncertainValue2
   {
   public:
      __device__ UncertainValue2();
      //__device__ ~UncertainValue2();
      __device__ UncertainValue2(double v, char source[], double dv);
      __device__ UncertainValue2(double v);
      __device__ UncertainValue2(double v, double dv);
      __device__ UncertainValue2(double v, Map::Map<String::String, double>& sigmas);
      __device__ UncertainValue2(UncertainValue2&);
      __device__ UncertainValue2& operator=(UncertainValue2&);

      __device__ void assignInitialValue(double);
      __device__ void assignComponent(String::String name, double sigma);
      __device__ double getComponent(String::String src);
      __device__ Map::Map<String::String, double> getComponents();
      __device__ bool hasComponent(String::String src);
      __device__ void renameComponent(String::String oldName, String::String newName);

      __device__ double doubleValue();
      __device__ bool isUncertain();
      __device__ double uncertainty();
      __device__ double variance();
      __device__ double fractionalUncertainty();
      __device__ bool equals(UncertainValue2& uv);
      //__device__ bool operator==(UncertainValue2&);

      __device__ int compareTo(UncertainValue2& o);
      __device__ bool lessThan(UncertainValue2& uv2);
      __device__ bool greaterThan(UncertainValue2& uv2);
      __device__ bool lessThanOrEqual(UncertainValue2& uv2);
      __device__ bool greaterThanOrEqual(UncertainValue2& uv2);

      __device__ UncertainValue2 sqrt();

      __device__ double variance(Correlations& corr);
      __device__ double uncertainty(Correlations& corr);

      __device__ void PrintSigmas();

   private:
      Map::Map<String::String, double> mSigmas;
      double mValue;
   };

   //extern __device__ const char DEFAULT[8];
   //extern __device__ int sDefIndex; // transient

   __device__ UncertainValue2 ONE();
   __device__ UncertainValue2 NaN();
   __device__ UncertainValue2 POSITIVE_INFINITY();
   __device__ UncertainValue2 NEGATIVE_INFINITY();
   __device__ UncertainValue2 ZERO();

   //extern __device__ const long long serialVersionUID;
   extern __device__ const int MAX_LEN;

   __device__ UncertainValue2 add(UncertainValue2 uvs[], int uvsLen);
   __device__ UncertainValue2 add(double a, UncertainValue2& uva, double b, UncertainValue2& uvb);
   __device__ UncertainValue2 subtract(UncertainValue2& uva, UncertainValue2& uvb);
   __device__ UncertainValue2 mean(UncertainValue2 uvs[], int uvsLen);
   __device__ UncertainValue2 weightedMean(UncertainValue2 uvs[], int uvsLen);
   __device__ UncertainValue2 min(UncertainValue2 uvs[], int uvsLen);
   __device__ UncertainValue2 max(UncertainValue2 uvs[], int uvsLen);
   __device__ UncertainValue2 add(UncertainValue2& v1, double v2);
   __device__ UncertainValue2 add(double v1, UncertainValue2& v2);
   __device__ UncertainValue2 add(UncertainValue2& v1, UncertainValue2& v2);
   __device__ UncertainValue2 multiply(double v1, UncertainValue2& v2);
   __device__ UncertainValue2 multiply(UncertainValue2& v1, UncertainValue2& v2);
   __device__ UncertainValue2 invert(UncertainValue2& v);
   __device__ UncertainValue2 divide(UncertainValue2& a, UncertainValue2& b);
   __device__ UncertainValue2 divide(double a, UncertainValue2& b);
   __device__ UncertainValue2 divide(UncertainValue2& a, double b);
   __device__ UncertainValue2 exp(UncertainValue2& x);
   __device__ UncertainValue2 log(UncertainValue2& v2);
   __device__ UncertainValue2 pow(UncertainValue2& v1, double n);
   __device__ UncertainValue2 sqrt(UncertainValue2& uv);
   __device__ LinkedList::Node<UncertainValue2>* quadratic(UncertainValue2& a, UncertainValue2& b, UncertainValue2& c);
   __device__ UncertainValue2 sqr(UncertainValue2& uv);
   __device__ UncertainValue2 negate(UncertainValue2& uv);
   __device__ UncertainValue2 atan(UncertainValue2& uv);
   __device__ UncertainValue2 atan2(UncertainValue2& y, UncertainValue2& x);
   __device__ UncertainValue2 positiveDefinite(UncertainValue2& uv);

   //__device__ void InitializeSpecialUncertainValues();
   __device__ bool AreEqual(UncertainValue2&, UncertainValue2&);
}

#endif