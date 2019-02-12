#ifndef _UNCERTAIN_VALUE_2_CUH_
#define _UNCERTAIN_VALUE_2_CUH_

#include <cuda_runtime.h>
#include <math_constants.h>

#include "..\..\..\..\Amphibian\String.cuh"
#include "..\..\..\..\Amphibian\LinkedList.cuh"

namespace UncertainValue2
{
   class UncertainValue2
   {
   public:
      __device__ UncertainValue2(double v, char source[], double dv);
      __device__ UncertainValue2(double v);
      __device__ UncertainValue2(double v, double dv);
      __device__ UncertainValue2(double v, LinkedListKV::Node<String::String, double>* sigmas);

      __device__ void assignComponent(String::String name, double sigma);
      __device__ double getComponent(String::String src);
      __device__ LinkedListKV::Node<String::String, double> * getComponents();
      __device__ bool hasComponent(String::String src);
      __device__ void renameComponent(String::String oldName, String::String newName);

      __device__ double doubleValue();
      __device__ bool isUncertain();
      __device__ double uncertainty();
      __device__ double variance();
      __device__ double fractionalUncertainty();
      __device__ bool equals(UncertainValue2 const * uv);

      __device__ int compareTo(UncertainValue2 o);
      __device__ bool lessThan(UncertainValue2 uv2);
      __device__ bool greaterThan(UncertainValue2 uv2);
      __device__ bool lessThanOrEqual(UncertainValue2 uv2);
      __device__ bool greaterThanOrEqual(UncertainValue2 uv2);

      __device__ UncertainValue2 sqrt();
      __device__ UncertainValue2 sqrt(UncertainValue2 uv);

      //class Correlations {
      //   static private class Key {
      //      private final String mSource1;
      //      private final String mSource2;

      //      Key(String src1, String src2) {
      //         mSource1 = src1;
      //         mSource2 = src2;
      //      }

      //      @Override
      //         public int hashCode() {
      //         return mSource1.hashCode() + mSource2.hashCode();
      //      }

      //      @Override
      //         public boolean equals(Object obj) {
      //         if (obj instanceof Key) {
      //            Key k2 = (Key)obj;
      //            return (mSource1.equals(k2.mSource1) && mSource2.equals(k2.mSource2))
      //               || (mSource1.equals(k2.mSource2) && mSource2.equals(k2.mSource1));
      //         }
      //         else
      //            return false;
      //      }
      //   }

      //   private final HashMap<Key, Double> mCorrelations;

      //   public Correlations() {
      //      mCorrelations = new HashMap<Key, Double>();
      //   }

      //   /**
      //   * Adds a correlation between src1 and src2 (order does not matter.)
      //   *
      //   * @param src1
      //   * @param src2
      //   * @param corr The correlation on the range [-1.0,1.0]
      //   */
      //   public void add(String src1, String src2, double corr) {
      //      assert(corr >= -1.0) && (corr <= 1.0);
      //      mCorrelations.put(new Key(src1, src2), Math2.bound(corr, -1.0, 1.0));
      //   }

      //   /**
      //   * Returns the correlation associated with src1 and src2 (order does not
      //   * matter) or zero if one has not been specified.
      //   *
      //   * @param src1
      //   * @param src2
      //   * @return [-1.0,1.0] with 0.0 as default
      //   */
      //   public double get(String src1, String src2) {
      //      Double r = mCorrelations.get(new Correlations.Key(src1, src2));
      //      return r == null ? 0.0 : r.doubleValue();
      //   }
      //}


   private:
      double mValue;
      LinkedListKV::Node<String::String, double>* mSigmas;
   };

   extern __device__ const char DEFAULT[8];
   extern __device__ int sDefIndex; // transient

   //extern __device__ const UncertainValue2 ONE;
   //extern __device__ const UncertainValue2 NaN(CUDART_NAN);
   //extern __device__ const UncertainValue2 POSITIVE_INFINITY(CUDART_INF);
   //extern __device__ const UncertainValue2 NEGATIVE_INFINITY(-CUDART_INF);
   //extern __device__ const UncertainValue2 ZERO(0.0);

   extern __device__ const long long serialVersionUID;
   extern __device__ const int MAX_LEN;

   __device__ UncertainValue2 add(LinkedList::Node<UncertainValue2> uvs);
   __device__ UncertainValue2 add(double a, UncertainValue2 uva, double b, UncertainValue2 uvb);
   __device__ UncertainValue2 subtract(UncertainValue2 uva, UncertainValue2 uvb);
   __device__ UncertainValue2 mean(LinkedList::Node<UncertainValue2>* uvs);
   __device__ UncertainValue2 weightedMean(LinkedList::Node<UncertainValue2>* cuv);
   __device__ UncertainValue2 min(LinkedList::Node<UncertainValue2>* uvs);
   __device__ UncertainValue2 max(LinkedList::Node<UncertainValue2>* uvs);
   __device__ UncertainValue2 add(UncertainValue2 v1, double v2);
   __device__ UncertainValue2 add(double v1, UncertainValue2 v2);
   __device__ UncertainValue2 add(UncertainValue2 v1, UncertainValue2 v2);
   __device__ UncertainValue2 multiply(double v1, UncertainValue2 v2);
   __device__ UncertainValue2 invert(UncertainValue2 v);
   __device__ UncertainValue2 divide(UncertainValue2 a, UncertainValue2 b);
   __device__ UncertainValue2 divide(double a, UncertainValue2 b);
   __device__ UncertainValue2 divide(UncertainValue2 a, double b);
   __device__ UncertainValue2 exp(UncertainValue2 x);
   __device__ UncertainValue2 log(UncertainValue2 v2);
   __device__ UncertainValue2 pow(UncertainValue2 v1, double n);
   __device__ UncertainValue2 sqrt(UncertainValue2 uv);
   __device__ LinkedList::Node<UncertainValue2>* quadratic(UncertainValue2 a, UncertainValue2 b, UncertainValue2 c);
   __device__ UncertainValue2 sqr(UncertainValue2 uv);
   __device__ UncertainValue2 negate(UncertainValue2 uv);
   __device__ UncertainValue2 atan(UncertainValue2 uv);
   __device__ UncertainValue2 atan2(UncertainValue2 y, UncertainValue2 x);
   __device__ UncertainValue2 positiveDefinite(UncertainValue2 uv);
}

#endif