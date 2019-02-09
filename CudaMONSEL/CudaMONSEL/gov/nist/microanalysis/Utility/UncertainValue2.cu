#include "UncertainValue2.cuh"
#include "..\..\..\..\Amphibian\LinkedList.cuh"
#include "..\..\..\..\Amphibian\Math.cuh"

#include <stdio.h>

namespace UncertainValue2
{
   __device__ const char DEFAULT[] = "Default";
   __device__ int sDefIndex = 0;

   __device__ const long long serialVersionUID = 119495064970078787L;
   __device__ const int MAX_LEN = 11;

   __device__ UncertainValue2::UncertainValue2(double v, double dv) : mValue(v), mSigmas(NULL)
   {
      char tmpName[MAX_LEN];
      String::IToA(tmpName, atomicAdd(&sDefIndex, 1));
      UncertainValue2::UncertainValue2(v, tmpName, dv);
   }
   
   __device__ UncertainValue2::UncertainValue2(double v) : mValue(v), mSigmas(NULL)
   {
      UncertainValue2::UncertainValue2(v, 0.0);
   }
   
   __device__ UncertainValue2::UncertainValue2(double v, char source[], double dv) : mValue(v), mSigmas(NULL)
   {
      assignComponent(source, dv);
   }
   
   __device__ UncertainValue2::UncertainValue2(double v, LinkedListKV::Node<String::String, double>* sigmas) : mValue(v), mSigmas(NULL)
   {
      while (sigmas != NULL) {
         assignComponent(sigmas->GetKey(), sigmas->GetValue());
         sigmas = sigmas->GetNext();
      }
   }

   __device__ double UncertainValue2::GetValue()
   {
      return mValue;
   }

   __device__ LinkedListKV::Node<String::String, double>* UncertainValue2::GetSigmas()
   {
      return mSigmas;
   }
   
   __device__ void UncertainValue2::assignComponent(String::String name, double sigma)
   {
      if (sigma != 0.0) {
         LinkedListKV::InsertHead<String::String, double>(&mSigmas, name, sigma);
      }
      else {
         LinkedListKV::Remove<String::String, double>(&mSigmas, name, String::AreEqual);
      }
   }
   
   __device__ double UncertainValue2::getComponent(String::String src)
   {
      auto v = LinkedListKV::GetValue<String::String, double>(mSigmas, src, String::AreEqual);
      return v != NULL ? v : 0.0;
   }

   __device__ LinkedListKV::Node<String::String, double> const * UncertainValue2::getComponents() const
   {
      return mSigmas;
   }
   
   __device__ bool UncertainValue2::hasComponent(String::String src)
   {
      return LinkedListKV::GetValue<String::String, double>(mSigmas, src, String::AreEqual) != NULL;
   }

   __device__ void UncertainValue2::renameComponent(String::String oldName, String::String newName)
   {
      if (LinkedListKV::ContainsKey<String::String, double>(mSigmas, newName, String::AreEqual)) {
         printf("A component named %s already exists.", newName.Get());
      }
      double val = LinkedListKV::Remove<String::String, double>(&mSigmas, oldName, String::AreEqual);
      if (val != NULL) {
         LinkedListKV::InsertHead<String::String, double>(&mSigmas, newName, val);
      }
   }

   //UncertainValue2 add(LinkedList::Node<UncertainValue2> uvs)
   //{
   //   HashSet<String> srcs = new HashSet<String>();
   //   double sum = 0.0;
   //   for (UncertainValue2 uv2 : uvs) {
   //      srcs.addAll(uv2.mSigmas.keySet());
   //      sum += uv2.mValue;
   //   }
   //   UncertainValue2 res = new UncertainValue2(sum);
   //   for (String src : srcs) {
   //      double unc = 0.0;
   //      // This seems right but is it????
   //      for (UncertainValue2 uv2 : uvs)
   //         unc += Math.signum(uv2.mValue) * uv2.getComponent(src);
   //      res.assignComponent(src, unc);
   //   }
   //   return res;
   //}

   //UncertainValue2 add(UncertainValue2[] uvs) {
   //   return add(Arrays.asList(uvs));
   //}

   __device__ UncertainValue2 add(double a, UncertainValue2 uva, double b, UncertainValue2 uvb)
   {
      UncertainValue2 res(a * uva.GetValue() + b * uvb.GetValue());
      LinkedList::Node<String::String>** srcs = NULL;
      AdvancedLinkedList::AddAllKeys<String::String, double>(srcs, uva.GetSigmas(), String::AreEqual);
      AdvancedLinkedList::AddAllKeys<String::String, double>(srcs, uvb.GetSigmas(), String::AreEqual);
      while (srcs != NULL) {
         String::String src;
         res.assignComponent(src, a * Math::signum(uva.GetValue()) * uva.getComponent(src) + b * Math::signum(uvb.GetValue()) * uvb.getComponent(src));
      }
      return res;
   }

   //static public UncertainValue2 subtract(UncertainValue2 uva, UncertainValue2 uvb) {
   //   return add(1.0, uva, -1.0, uvb);
   //}

   //public static UncertainValue2 mean(Collection<UncertainValue2> uvs) {
   //   return divide(add(uvs), uvs.size());
   //}

   ///**
   //* <p>
   //* Computes the variance weighted mean - the maximum likelyhood estimator of
   //* the mean under the assumption that the samples are independent and
   //* normally distributed.
   //* </p>
   //*
   //* @param cuv
   //* @return UncertainValue2
   //*/
   //static public UncertainValue2 weightedMean(Collection<UncertainValue2> cuv)
   //   throws UtilException{
   //   double varSum = 0.0, sum = 0.0;
   //   for (final UncertainValue2 uv : cuv) {
   //      final double ivar = 1.0 / uv.variance();
   //      if (Double.isNaN(ivar))
   //         throw new UtilException("Unable to compute the weighted mean when one or more datapoints have zero uncertainty.");
   //      varSum += ivar;
   //      sum += ivar * uv.doubleValue();
   //   }
   //   final double iVarSum = 1.0 / varSum;
   //   return Double.isNaN(iVarSum) ? UncertainValue2.NaN : new UncertainValue2(sum / varSum, "WM", Math.sqrt(1.0 / varSum));
   //}

   //   /**
   //   * Returns the uncertain value with the smallest value. When two uncertain
   //   * values are equal by value, the one with the larger uncertainty is
   //   * returned.
   //   *
   //   * @param uvs
   //   * @return UncertainValue2
   //   */
   //   public static UncertainValue2 min(Collection<UncertainValue2> uvs) {
   //   UncertainValue2 res = null;
   //   for (UncertainValue2 uv : uvs)
   //      if (res == null)
   //         res = uv;
   //      else if (uv.doubleValue() < res.doubleValue())
   //         res = uv;
   //      else if (uv.doubleValue() == res.doubleValue())
   //         if (uv.uncertainty() > res.uncertainty())
   //            res = uv;
   //      return res;
   //}

   ///**
   //* Returns the uncertain value with the largest value. When two uncertain
   //* values are equal by value, the one with the larger uncertainty is
   //* returned.
   //*
   //* @param uvs
   //* @return UncertainValue2
   //*/
   //public static UncertainValue2 max(Collection<UncertainValue2> uvs) {
   //   UncertainValue2 res = null;
   //   for (UncertainValue2 uv : uvs)
   //      if (res == null)
   //         res = uv;
   //      else if (uv.doubleValue() > res.doubleValue())
   //         res = uv;
   //      else if (uv.doubleValue() == res.doubleValue())
   //         if (uv.uncertainty() > res.uncertainty())
   //            res = uv;
   //      return res;
   //}

   ///**
   //* Add a quantity known without error to an UncertainValue2.
   //*
   //* @param v1
   //* @param v2
   //* @return An UncertainValue2
   //*/
   //static public UncertainValue2 add(UncertainValue2 v1, double v2) {
   //   return new UncertainValue2(v1.mValue + v2, v1.mSigmas);
   //}

   ///**
   //* Add a quantity known without error to an UncertainValue2.
   //*
   //* @param v1
   //* @param v2
   //* @return An UncertainValue2
   //*/
   //static public UncertainValue2 add(double v1, UncertainValue2 v2) {
   //   return new UncertainValue2(v2.mValue + v1, v2.mSigmas);
   //}

   ///**
   //* Add a quantity known without error to an UncertainValue2.
   //*
   //* @param v1
   //* @param v2
   //* @return An UncertainValue2
   //*/
   //static public UncertainValue2 add(UncertainValue2 v1, UncertainValue2 v2) {
   //   return add(1.0, v1, 1.0, v2);
   //}

   ///**
   //* Multiply a constant times an UncertainValue2
   //*
   //* @param v1
   //* @param v2
   //* @return An UncertainValue2
   //*/
   //static public UncertainValue2 multiply(double v1, UncertainValue2 v2) {
   //   assert v2.uncertainty() >= 0.0;
   //   UncertainValue2 res = new UncertainValue2(v1 * v2.mValue);
   //   for (Map.Entry<String, Double> me : v2.mSigmas.entrySet())
   //      res.assignComponent(me.getKey(), v1 * me.getValue());
   //   return res;
   //}

   ///**
   //* Multiply two uncertain values
   //*
   //* @param v1
   //* @param v2
   //* @return An UncertainValue2
   //*/
   //static public UncertainValue2 multiply(UncertainValue2 v1, UncertainValue2 v2) {
   //   UncertainValue2 res = new UncertainValue2(v1.mValue * v2.mValue);
   //   Set<String> srcs = new TreeSet<String>();
   //   srcs.addAll(v1.mSigmas.keySet());
   //   srcs.addAll(v2.mSigmas.keySet());
   //   for (String src : srcs)
   //      res.assignComponent(src, v1.mValue * v2.getComponent(src) + v2.mValue * v1.getComponent(src));
   //   return res;
   //}

   //static public UncertainValue2 invert(UncertainValue2 v) {
   //   return UncertainValue2.divide(1.0, v);
   //}

   ///**
   //* Divide two uncertain values.
   //*
   //* @param a
   //* @param b
   //* @return An UncertainValue2
   //*/
   //static public UncertainValue2 divide(UncertainValue2 a, UncertainValue2 b) {
   //   UncertainValue2 res = new UncertainValue2(a.mValue / b.mValue);
   //   if (!Double.isNaN(res.doubleValue())) {
   //      Set<String> srcs = new TreeSet<String>();
   //      srcs.addAll(a.mSigmas.keySet());
   //      srcs.addAll(b.mSigmas.keySet());
   //      final double ua = Math.abs(1.0 / b.mValue);
   //      final double ub = Math.abs(a.mValue / (b.mValue * b.mValue));
   //      for (String src : srcs)
   //         res.assignComponent(src, ua * a.getComponent(src) + ub * b.getComponent(src));
   //   }
   //   return res;
   //}

   //static public UncertainValue2 divide(double a, UncertainValue2 b) {
   //   UncertainValue2 res = new UncertainValue2(a / b.mValue);
   //   if (!Double.isNaN(res.doubleValue())) {
   //      final double ub = Math.abs(a / (b.mValue * b.mValue));
   //      for (Map.Entry<String, Double> me : b.mSigmas.entrySet())
   //         res.assignComponent(me.getKey(), ub * me.getValue().doubleValue());
   //   }
   //   return res;
   //}

   //static public UncertainValue2 divide(UncertainValue2 a, double b) {
   //   if (!Double.isNaN(1.0 / b)) {
   //      UncertainValue2 res = new UncertainValue2(a.doubleValue() / b);
   //      final double ua = Math.abs(1.0 / b);
   //      for (Map.Entry<String, Double> me : a.mSigmas.entrySet())
   //         res.assignComponent(me.getKey(), ua * me.getValue().doubleValue());
   //      return res;
   //   }
   //   else
   //      return UncertainValue2.NaN;
   //}

   ///**
   //* Compute the exponental function of an UncertainValue2.
   //*
   //* @param x
   //* @return An UncertainValue2
   //*/

   //static public UncertainValue2 exp(UncertainValue2 x) {
   //   assert !Double.isNaN(x.mValue) : x.toString();
   //   double ex = Math.exp(x.mValue);
   //   UncertainValue2 res = new UncertainValue2(ex);
   //   for (Map.Entry<String, Double> me : x.mSigmas.entrySet())
   //      res.assignComponent(me.getKey(), ex * me.getValue().doubleValue());
   //   return res;
   //}

   ///**
   //* Compute the natural logarithm of an UncertainValue2.
   //*
   //* @param v2
   //* @return An UncertainValue2
   //*/
   //static public UncertainValue2 log(UncertainValue2 v2) {
   //   double tmp = 1.0 / v2.mValue;
   //   final double lv = Math.log(v2.mValue);
   //   if (!(Double.isNaN(tmp) || Double.isNaN(lv))) {
   //      UncertainValue2 res = new UncertainValue2(lv);
   //      for (Map.Entry<String, Double> me : v2.mSigmas.entrySet())
   //         res.assignComponent(me.getKey(), tmp * me.getValue());
   //      return res;
   //   }
   //   else
   //      return UncertainValue2.NaN;
   //}

   ///**
   //* Compute an UncertainValue2 raised to the specified power.
   //*
   //* @param v1 The value
   //* @param n The exponent
   //* @return An UncertainValue2
   //*/
   //static public UncertainValue2 pow(UncertainValue2 v1, double n) {
   //   if (v1.mValue != 0.0) {
   //      final double f = Math.pow(v1.mValue, n);
   //      final double df = n * Math.pow(v1.mValue, n - 1.0);
   //      UncertainValue2 res = new UncertainValue2(f);
   //      for (Map.Entry<String, Double> me : v1.mSigmas.entrySet())
   //         res.assignComponent(me.getKey(), me.getValue() * df);
   //      return res;
   //   }
   //   else
   //      return UncertainValue2.ZERO;
   //}

   //public UncertainValue2 sqrt() {
   //   return pow(this, 0.5);
   //}

   //public static UncertainValue2 sqrt(UncertainValue2 uv) {
   //   return pow(uv, 0.5);
   //}

   //public static UncertainValue2[] quadratic(UncertainValue2 a, UncertainValue2 b, UncertainValue2 c) {
   //   // q=-0.5*(b+signum(b)*sqrt(pow(b,2.0)-4*a*c))
   //   // return [ q/a, c/q ]
   //   UncertainValue2 r = UncertainValue2.add(1.0, UncertainValue2.pow(b, 2.0), -4.0, UncertainValue2.multiply(a, c));
   //   if (r.doubleValue() > 0.0) {
   //      UncertainValue2 q = UncertainValue2.multiply(-0.5, UncertainValue2.add(b, UncertainValue2.multiply(Math.signum(b.mValue), r.sqrt())));
   //      return new UncertainValue2[] {
   //         UncertainValue2.divide(q, a),
   //            UncertainValue2.divide(c, q)
   //      };
   //   }
   //   else
   //      return null;

   //}
   
   __device__ double UncertainValue2::doubleValue()
   {
      return mValue;
   }
   
   __device__ bool UncertainValue2::isUncertain()
   {
      return mSigmas != NULL;
   }
   
   __device__ double UncertainValue2::uncertainty()
   {
      return sqrt(variance());
   }
   
   __device__ double UncertainValue2::variance()
   {
      double sigma2 = 0.0;
      LinkedListKV::Node<String::String, double>* head = mSigmas;
      while (head != NULL) {
         double v = head->GetValue();
         sigma2 += v * v;
         head = head->GetNext();
      }
      return sigma2;
   }
   
   __device__ double UncertainValue2::fractionalUncertainty()
   {
      return (isnan(1.0 / mValue)) ? CUDART_NAN : abs(uncertainty() / mValue);
   }
   
   __device__ bool UncertainValue2::equals(UncertainValue2 const * obj)
   {
      if (this == obj) {
         return true;
      }
      if (obj == NULL) {
         return false;
      }
      UncertainValue2 other = (UncertainValue2)*obj;
      return LinkedListKV::AreEquivalentSets<String::String, double>(mSigmas, other.mSigmas, String::AreEqual, [](double a, double b) { return a == b; }) && (mValue == other.mValue);
   }
   
   __device__ int UncertainValue2::compareTo(UncertainValue2 o)
   {
      return (mValue == o.mValue) && (uncertainty() == o.uncertainty());
   }
   
   __device__ bool UncertainValue2::lessThan(UncertainValue2 uv2)
   {
      return mValue < uv2.mValue;
   }
   
   __device__ bool UncertainValue2::greaterThan(UncertainValue2 uv2)
   {
      return mValue > uv2.mValue;
   }
   
   __device__ bool UncertainValue2::lessThanOrEqual(UncertainValue2 uv2)
   {
      return mValue <= uv2.mValue;
   }
   
   __device__ bool UncertainValue2::greaterThanOrEqual(UncertainValue2 uv2)
   {
      return mValue >= uv2.mValue;
   }

   //__device__ UncertainValue2 sqr(UncertainValue2 uv)
   //{
   //   return pow(uv, 2.0);
   //}
   
   //__device__ UncertainValue2 UncertainValue2::negate(UncertainValue2 uv)
   //{
   //   UncertainValue2 ret = UncertainValue2(-uv.mValue, uv.mSigmas);
   //   return ret;
   //}
   
   //__device__ UncertainValue2 UncertainValue2::atan(UncertainValue2 uv)
   //{
   //   double f = atan(uv.doubleValue());
   //   double df = 1.0 / (1.0 + uv.doubleValue() * uv.doubleValue());
   //
   //   if (!(isnan(f) || isnan(df))) {
   //      UncertainValue2 res(f);
   //      LinkedListKV::Node<String::String, double>* sigmas = uv.mSigmas;
   //      while (sigmas != NULL) {
   //         res.assignComponent(sigmas->GetKey(), df * sigmas->GetValue());
   //         sigmas = sigmas->GetNext();
   //      }
   //      return res;
   //   }
   //   else {
   //      UncertainValue2 nan(CUDART_NAN);
   //      return nan;
   //   }
   //}

   //__device__ UncertainValue2 UncertainValue2::atan2(UncertainValue2 y, UncertainValue2 x)
   //{
   //   double f = atan2(y.doubleValue(), x.doubleValue());
   //   double df = 1.0 / (1.0 + sqr(y.doubleValue() / x.doubleValue()));
   //
   //   if (!(isnan(f) || isnan(df))) {
   //      UncertainValue2 res(f);
   //      LinkedListKV::Node<String::String, double>* sigmas;
   //      while (sigmas != NULL) {
   //         res.assignComponent(sigmas->GetKey(), df * sigmas->GetValue());
   //         sigmas = sigmas->GetNext();
   //      }
   //      return res;
   //   }
   //   else {
   //      UncertainValue2 nan(CUDART_NAN);
   //      return nan;
   //   }
   //}

   //__device__ UncertainValue2 UncertainValue2::positiveDefinite(UncertainValue2 uv)
   //{
   //   UncertainValue2 ret(0.0, uv.mSigmas);
   //   return uv.doubleValue() >= 0.0 ? uv : ret;
   //}
}