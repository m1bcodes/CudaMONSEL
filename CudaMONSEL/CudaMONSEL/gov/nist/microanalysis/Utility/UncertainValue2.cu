#include "UncertainValue2.cuh"
#include "..\..\..\..\Amphibian\LinkedList.cuh"
#include "..\..\..\..\Amphibian\Math.cuh"

#include <stdio.h>
#include <math.h>

extern __device__ double __longlong_as_double(long long int);

extern __device__ int atomicAdd(int* address, int val);

namespace UncertainValue2
{
   __device__ const char DEFAULT[] = "Default";
   __device__ int sDefIndex = 0;

   __device__ const long long serialVersionUID = 119495064970078787L;
   __device__ const int MAX_LEN = 11;

   __device__ const Hasher::pHasher DefaultHasher = Hasher::APHash;

   __device__ UncertainValue2::UncertainValue2() : mSigmas(DefaultHasher, String::AreEqual)
   {
   }

   __device__ UncertainValue2::UncertainValue2(double v, double dv) : mValue(v), mSigmas(DefaultHasher, String::AreEqual)
   {
      char tmpName[MAX_LEN];
      String::IToA(tmpName, atomicAdd(&sDefIndex, 1));
      assignComponent(tmpName, dv);
   }

   __device__ UncertainValue2::UncertainValue2(double v) : mValue(v), mSigmas(DefaultHasher, String::AreEqual)
   {
   }

   __device__ UncertainValue2::UncertainValue2(double v, char source[], double dv) : mValue(v), mSigmas(DefaultHasher, String::AreEqual)
   {
      assignComponent(source, dv);
   }

   __device__ UncertainValue2::UncertainValue2(double v, Map::Map<String::String, double> sigmas) : mValue(v), mSigmas(DefaultHasher, String::AreEqual)
   {
      if (*((int*)&sigmas) != NULL) {
         mSigmas.DeepCopy(sigmas);
      }
   }

   __device__ UncertainValue2::UncertainValue2(UncertainValue2& other) : mValue(other.doubleValue()), mSigmas(DefaultHasher, String::AreEqual)
   {
      mSigmas.DeepCopy(other.getComponents());
   }

   __device__ UncertainValue2 ONE()
   {
      return UncertainValue2(1.0);
   }

   __device__ UncertainValue2 NaN()
   {
      return UncertainValue2(CUDART_NAN);
   }

   __device__ UncertainValue2 POSITIVE_INFINITY()
   {
      return UncertainValue2(CUDART_INF);
   }

   __device__ UncertainValue2 NEGATIVE_INFINITY()
   {
      return UncertainValue2(-CUDART_INF);
   }

   __device__ UncertainValue2 ZERO()
   {
      return UncertainValue2(0.0);
   }

   __device__ UncertainValue2& UncertainValue2::operator=(UncertainValue2& other)
   {
      mValue = other.doubleValue();
      mSigmas.DeepCopy(other.getComponents());

      return *this;
   }

   __device__ void UncertainValue2::assignInitialValue(double v)
   {
      mValue = v;
   }

   __device__ void UncertainValue2::assignComponent(String::String name, double sigma)
   {
      if (sigma != 0.0) {
         mSigmas.Put(name, sigma);
      }
      else {
         mSigmas.Remove(name);
      }
   }

   __device__ double UncertainValue2::getComponent(String::String src)
   {
      auto v = mSigmas.GetValue(src);
      return v != NULL ? v : 0.0;
   }

   __device__ Map::Map<String::String, double> UncertainValue2::getComponents()
   {
      return mSigmas;
   }

   __device__ bool UncertainValue2::hasComponent(String::String src)
   {
      return getComponent(src) != 0.0;
   }

   __device__ void UncertainValue2::renameComponent(String::String oldName, String::String newName)
   {
      // if (LinkedListKV::ContainsKey<String::String, double>(mSigmas, newName, String::AreEqual)) {
      if (mSigmas.ContainsKey(newName)) {
         printf("A component named %s already exists.", newName.Get());
         return;
      }
      double val = mSigmas.Remove(oldName);
      if (val != NULL) {
         mSigmas.Put(newName, val);
      }
   }

   __device__ UncertainValue2 add(UncertainValue2 uvs[], int uvsLen)
   {
      Set::Set<String::String> keys(DefaultHasher, String::AreEqual);
      double sum = 0.0;
      for (int k = 0; k < uvsLen; ++k) {
         sum += uvs[k].doubleValue();
         keys.Add(uvs[k].getComponents().GetKeys());
      }
      UncertainValue2 res(sum);

      Set::Iterator<String::String> itr(keys);

      while (itr.HasNext()) {
         auto src = itr.GetValue();
         itr.Next();
         double unc = 0.0;
         // This seems right but is it????
         for (int k = 0; k < uvsLen; ++k) {
            auto uv = uvs[k];
            unc += uv.getComponent(src) * copysign(1.0, uv.doubleValue());
         }
         res.assignComponent(src, unc);
      }
      return res;
   }

   __device__ UncertainValue2 add(double a, UncertainValue2 uva, double b, UncertainValue2 uvb)
   {
      Set::Set<String::String> keys(DefaultHasher, String::AreEqual);
      UncertainValue2 res(a * uva.doubleValue() + b * uvb.doubleValue());
      keys.Add(uva.getComponents().GetKeys());
      keys.Add(uvb.getComponents().GetKeys());

      Set::Iterator<String::String> itr(keys);
      while (itr.HasNext()) {
         String::String src = itr.GetValue();
         res.assignComponent(src, a * copysign(1.0, uva.doubleValue()) * uva.getComponent(src) + b * copysign(1.0, uvb.doubleValue()) * uvb.getComponent(src));
         itr.Next();
      }
      return res;
   }

   __device__ UncertainValue2 subtract(UncertainValue2 uva, UncertainValue2 uvb)
   {
      return add(1.0, uva, -1.0, uvb);
   }

   __device__ UncertainValue2 mean(UncertainValue2 uvs[], int uvsLen)
   {
      return divide(add(uvs, uvsLen), (double)uvsLen);
   }

   __device__ UncertainValue2 weightedMean(UncertainValue2 cuv[], int uvsLen)
   {
      double varSum = 0.0, sum = 0.0;

      for (int k = 0; k < uvsLen; ++k) {
         auto uv = cuv[k];
         const double ivar = 1.0 / uv.variance();
         if (isnan(ivar) || isinf(ivar)) {
            printf("%s\n", "Unable to compute the weighted mean when one or more datapoints have zero uncertainty.");
            return NULL;
         }
         varSum += ivar;
         sum += ivar * uv.doubleValue();
      }
      const double iVarSum = 1.0 / varSum;
      return (isnan(iVarSum) || isinf(iVarSum)) ? NULL : UncertainValue2(sum / varSum, "WM", ::sqrt(1.0 / varSum));
   }

   __device__ UncertainValue2 min(UncertainValue2 uvs[], int uvsLen)
   {
      if (uvsLen == 0) {
         return NULL;
      }
      UncertainValue2 res = uvs[0];

      for (int k = 0; k < uvsLen; ++k) {
         auto uv = uvs[k];
         if (uv.doubleValue() < res.doubleValue()) {
            res = uv;
         }
         else if (uv.doubleValue() == res.doubleValue()) {
            if (uv.uncertainty() > res.uncertainty()) {
               res = uv;
            }
         }
      }
      return res;
   }

   __device__ UncertainValue2 max(UncertainValue2 uvs[], int uvsLen)
   {
      if (uvs == 0) {
         return NULL;
      }
      UncertainValue2 res = uvs[0];

      for (int k = 0; k < uvsLen; ++k) {
         auto uv = uvs[k];
         if (uv.doubleValue() > res.doubleValue()) {
            res = uv;
         }
         else if (uv.doubleValue() == res.doubleValue()) {
            if (uv.uncertainty() > res.uncertainty()) {
               res = uv;
            }
         }
      }
      return res;
   }

   __device__ UncertainValue2 add(UncertainValue2 v1, double v2)
   {
      return UncertainValue2(v1.doubleValue() + v2, v1.getComponents());
   }

   __device__ UncertainValue2 add(double v1, UncertainValue2 v2)
   {
      return UncertainValue2(v2.doubleValue() + v1, v2.getComponents());
   }

   __device__ UncertainValue2 add(UncertainValue2 v1, UncertainValue2 v2)
   {
      return add(1.0, v1, 1.0, v2);
   }

   __device__ UncertainValue2 multiply(double v1, UncertainValue2 v2)
   {
      if (v2.uncertainty() < 0.0) {
         printf("Error: v2.uncertainty() < 0.0");
         return NULL;
      }
      UncertainValue2 res(v1 * v2.doubleValue());

      auto m = v2.getComponents();
      Map::Iterator<String::String, double> itr(m);
      while (itr.HasNext()) {
         res.assignComponent(itr.GetKey(), v1 * itr.GetValue());
         itr.Next();
      }
      return res;
   }

   __device__ UncertainValue2 multiply(UncertainValue2 v1, UncertainValue2 v2)
   {
      Set::Set<String::String> keys(DefaultHasher, String::AreEqual);
      keys.Add(v1.getComponents().GetKeys());
      keys.Add(v2.getComponents().GetKeys());

      UncertainValue2 res(v1.doubleValue() * v2.doubleValue());
      Set::Iterator<String::String> itr(keys);
      while (itr.HasNext()) {
         auto src = itr.GetValue();
         res.assignComponent(src, v1.doubleValue() * v2.getComponent(src) + v2.doubleValue() * v1.getComponent(src));
         itr.Next();
      }

      return res;
   }

   __device__ UncertainValue2 invert(UncertainValue2 v)
   {
      return divide(1.0, v);
   }

   __device__ UncertainValue2 divide(UncertainValue2 v1, UncertainValue2 v2)
   {
      UncertainValue2 res(v1.doubleValue() / v2.doubleValue());
      if (!(isnan(res.doubleValue()) || isinf(res.doubleValue()))) {
         Set::Set<String::String> keys(DefaultHasher, String::AreEqual);
         keys.Add(v1.getComponents().GetKeys());
         keys.Add(v2.getComponents().GetKeys());

         const double ua = fabs(1.0 / v2.doubleValue());
         const double ub = fabs(v1.doubleValue() / (v2.doubleValue() * v2.doubleValue()));

         Set::Iterator<String::String> itr(keys);
         while (itr.HasNext()) {
            auto src = itr.GetValue();
            res.assignComponent(src, ua * v1.getComponent(src) + ub * v2.getComponent(src));
            itr.Next();
         }
      }
      return res;
   }

   __device__ UncertainValue2 divide(double a, UncertainValue2 b)
   {
      UncertainValue2 res(a / b.doubleValue());
      if (!(isnan(res.doubleValue()) || isinf(res.doubleValue()))) {
         const double ub = fabs(a / (b.doubleValue() * b.doubleValue()));

         auto m = b.getComponents();
         Map::Iterator<String::String, double> itr(m);
         while (itr.HasNext()) {
            res.assignComponent(itr.GetKey(), ub * itr.GetValue());
            itr.Next();
         }
      }
      return res;
   }

   __device__ UncertainValue2 divide(UncertainValue2 a, double b)
   {
      if (isnan(1.0 / b)) {
         return UncertainValue2(CUDART_NAN);
      }
      if (isinf(1.0 / b)) {
         return UncertainValue2(CUDART_INF);
      }
      UncertainValue2 res(a.doubleValue() / b);
      const double ua = fabs(1.0 / b);

      auto m = a.getComponents();
      Map::Iterator<String::String, double> itr(m);
      while (itr.HasNext()) {
         res.assignComponent(itr.GetKey(), ua * itr.GetValue());
         itr.Next();
      }
      return res;
   }

   __device__ UncertainValue2 exp(UncertainValue2 x)
   {
      if (isnan(x.doubleValue()) || isinf(x.doubleValue())) {
         printf("exp: invalid value\n");
         return NULL;
      }

      double ex = ::exp(x.doubleValue());
      UncertainValue2 res(ex);

      auto m = x.getComponents();
      Map::Iterator<String::String, double> itr(m);
      while (itr.HasNext()) {
         res.assignComponent(itr.GetKey(), ex * itr.GetValue());
         itr.Next();
      }
      return res;
   }

   __device__ UncertainValue2 log(UncertainValue2 v2)
   {
      double tmp = 1.0 / v2.doubleValue();
      const double lv = ::log(v2.doubleValue());
      if (isnan(tmp) || isnan(lv)) {
         return UncertainValue2(CUDART_NAN);
      }
      if (isinf(tmp) || isinf(lv)) {
         return UncertainValue2(CUDART_INF);
      }
      UncertainValue2 res(lv);

      auto m = v2.getComponents();
      Map::Iterator<String::String, double> itr(m);
      while (itr.HasNext()) {
         res.assignComponent(itr.GetKey(), tmp * itr.GetValue());
         itr.Next();
      }
      return res;
   }

   __device__ UncertainValue2 pow(UncertainValue2 v1, double n)
   {
      if (v1.doubleValue() == 0.0) {
         return UncertainValue2(0.0);
      }
      const double f = ::pow(v1.doubleValue(), n);
      const double df = n * ::pow(v1.doubleValue(), n - 1.0);
      UncertainValue2 res(f);

      auto m = v1.getComponents();
      Map::Iterator<String::String, double> itr(m);
      while (itr.HasNext()) {
         res.assignComponent(itr.GetKey(), itr.GetValue() * df);
         itr.Next();
      }

      return res;
   }

   __device__ UncertainValue2 UncertainValue2::sqrt()
   {
      return pow(*this, 0.5);
   }

   __device__ UncertainValue2 sqrt(UncertainValue2 uv)
   {
      return pow(uv, 0.5);
   }

   __device__ LinkedList::Node<UncertainValue2>* quadratic(UncertainValue2 a, UncertainValue2 b, UncertainValue2 c)
   {
      // q=-0.5*(b+signum(b)*sqrt(pow(b,2.0)-4*a*c))
      // return [ q/a, c/q ]
      UncertainValue2 r = add(1.0, pow(b, 2.0), -4.0, multiply(a, c));
      if (r.doubleValue() <= 0.0) {
         return NULL;
      }
      UncertainValue2 q = multiply(-0.5, add(b, multiply(copysign(1.0, b.doubleValue()), r.sqrt())));
      LinkedList::Node<UncertainValue2>* head = NULL;
      LinkedList::InsertHead(&head, divide(q, a));
      LinkedList::InsertHead(&head, divide(c, q));
      return head;
   }

   __device__ double UncertainValue2::doubleValue()
   {
      return mValue;
   }

   __device__ bool UncertainValue2::isUncertain()
   {
      return !mSigmas.IsEmpty();
   }

   __device__ double UncertainValue2::uncertainty()
   {
      return ::sqrt(variance());
   }

   __device__ double UncertainValue2::variance()
   {
      return mSigmas.Aggregate([](double a) { return a*a; });
   }

   __device__ double UncertainValue2::fractionalUncertainty()
   {
      if (isnan(1.0 / mValue)) {
         return CUDART_NAN;
      }
      if (isinf(1.0 / mValue)) {
         return CUDART_INF;
      }
      return fabs(uncertainty() / mValue);
   }

   __device__ bool UncertainValue2::equals(UncertainValue2& other)
   {
      if (*((int*)&other) == NULL) {
         return false;
      }
      if (this == &other) {
         return true;
      }
      return true;
      //return LinkedListKV::AreEquivalentSets<String::String, double>(mSigmas, other.getComponents(), String::AreEqual, [](double a, double b) { return a == b; }) && (mValue == other.doubleValue());
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

   __device__ UncertainValue2 sqr(UncertainValue2 uv)
   {
      return pow(uv, 2.0);
   }

   __device__ UncertainValue2 negate(UncertainValue2 uv)
   {
      return UncertainValue2(-uv.doubleValue(), uv.getComponents());
   }

   __device__ UncertainValue2 atan(UncertainValue2 uv)
   {
      double f = ::atan(uv.doubleValue());
      double df = 1.0 / (1.0 + uv.doubleValue() * uv.doubleValue());

      if (isnan(f)) {
         return UncertainValue2(CUDART_NAN);
      }
      if (isinf(df)) {
         return UncertainValue2(CUDART_INF);
      }
      UncertainValue2 res(f);

      auto m = uv.getComponents();
      Map::Iterator<String::String, double> itr(m);
      while (itr.HasNext()) {
         res.assignComponent(itr.GetKey(), df * itr.GetValue());
         itr.Next();
      }
      return res;
   }

   __device__ UncertainValue2 atan2(UncertainValue2 y, UncertainValue2 x)
   {
      double f = ::atan2(y.doubleValue(), x.doubleValue());
      double df = 1.0 / (1.0 + ::pow(y.doubleValue() / x.doubleValue(), 2.0));

      if (isnan(f)) {
         return UncertainValue2(CUDART_NAN);
      }
      if (isinf(df)) {
         return UncertainValue2(CUDART_INF);
      }
      UncertainValue2 res(f);

      auto m = divide(y, x).getComponents();
      Map::Iterator<String::String, double> itr(m);
      while (itr.HasNext()) {
         res.assignComponent(itr.GetKey(), df * itr.GetValue());
         itr.Next();
      }
      return res;
   }

   __device__ UncertainValue2 positiveDefinite(UncertainValue2 uv)
   {
      UncertainValue2 ret(0.0, uv.getComponents());
      return uv.doubleValue() >= 0.0 ? uv : ret;
   }

   __device__ Key::Key(String::String src1, String::String src2)
   {
      mSource1 = src1;
      mSource2 = src2;
   }

   __device__ bool Key::operator==(Key& k2)
   {
      return (mSource1 == k2.mSource1 && mSource2 == k2.mSource2) || (mSource1 == k2.mSource2 && mSource2 == k2.mSource1);
   }

   __device__ bool Key::AreEqual(Key k1, Key k2)
   {
      return k1 == k2;
   }

   __device__ Correlations::Correlations() : mCorrelations(DefaultHasher, Key::AreEqual)
   {
   }

   __device__ void Correlations::add(String::String src1, String::String src2, double corr)
   {
      if (!(corr >= -1.0) && (corr <= 1.0)) {
         printf("%s\n", "Correlations::add: invalid bound");
         return;
      }
      corr = ::fmax(corr, 1.0);
      corr = ::fmin(corr, -1.0);
      mCorrelations.Put(Key(src1, src2), corr);
   }

   __device__ double Correlations::get(String::String src1, String::String src2)
   {
      double r = mCorrelations.GetValue(Key(src1, src2));
      return r == NULL ? 0.0 : r;
   }

   __device__ double UncertainValue2::variance(Correlations corr)
   {
      Map::Map<String::String, double> sigmas = getComponents();

      Set::Set<String::String> keys(DefaultHasher, String::AreEqual);
      keys.Add(sigmas.GetKeys());

      Set::Iterator<String::String> itr1(keys);
      double res = 0.0;
      while (itr1.HasNext()) {
         String::String key = itr1.GetValue();
         auto val = sigmas.GetValue(key);
         res += val * val;
         itr1.Next();
      }

      Set::Iterator<String::String> itr2(keys);

      while (itr1.HasNext()) {
         itr2 = itr1;
         itr2.Next();
         while (itr2.HasNext()) {
            auto key1 = itr1.GetValue();
            auto key2 = itr2.GetValue();
            auto sigma1 = sigmas.GetValue(key1);
            auto sigma2 = sigmas.GetValue(key2);
            res += 2.0 * sigma1 * sigma2 * corr.get(key1, key2);
            itr2.Next();
         }
         itr1.Next();
      }
      return res;
   }

   __device__ double UncertainValue2::uncertainty(Correlations corr)
   {
      return ::sqrt(variance(corr));
   }

   __device__ void UncertainValue2::PrintSigmas()
   {
      Map::Iterator <String::String, double> itr(mSigmas);
      while (itr.HasNext()) {
         printf("%s: %lf\n", itr.GetKey().Get(), itr.GetValue());
         itr.Next();
      }
   }

   __device__ bool AreEqual(UncertainValue2 a, UncertainValue2 b)
   {
      return a.equals(b);
   }
}
