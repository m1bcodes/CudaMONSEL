//#include "UncertainValue2.cuh"
//
//#include <stdio.h>
//#include <math.h>
//
//extern __device__ double __longlong_as_double(long long int);
//
//extern __device__ int atomicAdd(int* address, int val);
//
//namespace UncertainValue2
//{
//   __device__ const char DEFAULT[] = "Default";
//   __device__ int sDefIndex = 0;
//
//   __device__ const long long serialVersionUID = 119495064970078787L;
//   __device__ const int MAX_LEN = 11;
//
//   __device__ bool doubleCmp(double& a, double& b)
//   {
//      return a == b;
//   }
//
//   //__device__ Map::Map<String::String, double>::pValCmp doubleCmp = [](double a, double b) { return a == b; };
//
//   __device__ UncertainValue2::UncertainValue2()
//   {
//   }
//
//   __device__ UncertainValue2::UncertainValue2(double v, double dv) : mValue(v)
//   {
//      char tmpName[MAX_LEN];
//      String::IToA(tmpName, atomicAdd(&sDefIndex, 1));
//      assignComponent(tmpName, dv);
//   }
//
//   __device__ UncertainValue2::UncertainValue2(double v) : mValue(v)
//   {
//   }
//
//   __device__ UncertainValue2::UncertainValue2(double v, char source[], double dv) : mValue(v)
//   {
//      assignComponent(source, dv);
//   }
//
//   __device__ UncertainValue2::UncertainValue2(double v, ComponentMap& sigmas) : mValue(v)
//   {
//      mSigmas.DeepCopy(sigmas);
//   }
//
//   __device__ UncertainValue2::UncertainValue2(const UncertainValue2& other) : mValue(other.doubleValue())
//   {
//      if (&other == this) return;
//      mSigmas.DeepCopy(other.mSigmas);
//   }
//
//   __device__ UncertainValue2 ONE()
//   {
//      return UncertainValue2(1.0);
//   }
//
//   __device__ UncertainValue2 NaN()
//   {
//      return UncertainValue2(CUDART_NAN);
//   }
//
//   __device__ UncertainValue2 POSITIVE_INFINITY()
//   {
//      return UncertainValue2(CUDART_INF);
//   }
//
//   __device__ UncertainValue2 NEGATIVE_INFINITY()
//   {
//      return UncertainValue2(-CUDART_INF);
//   }
//
//   __device__ UncertainValue2 ZERO()
//   {
//      return UncertainValue2(0.0);
//   }
//
//   __device__ UncertainValue2& UncertainValue2::operator=(UncertainValue2& other)
//   {
//      mValue = other.doubleValue();
//      mSigmas.DeepCopy(other.getComponents());
//
//      return *this;
//   }
//
//   __device__ void UncertainValue2::assignInitialValue(double v)
//   {
//      mValue = v;
//   }
//
//   __device__ void UncertainValue2::assignComponent(String::String name, double sigma)
//   {
//      if (sigma != 0.0) {
//         mSigmas.Put(name, sigma);
//      }
//      else {
//         mSigmas.Remove(name);
//      }
//   }
//
//   __device__ double UncertainValue2::getComponent(String::String src)
//   {
//      double v;
//      return mSigmas.GetValue(src, v) ? v : 0.0;
//   }
//
//   __device__ UncertainValue2::ComponentMap& UncertainValue2::getComponents()
//   {
//      return mSigmas;
//   }
//
//   __device__ bool UncertainValue2::hasComponent(String::String src)
//   {
//      return getComponent(src) != 0.0;
//   }
//
//   __device__ void UncertainValue2::renameComponent(String::String oldName, String::String newName)
//   {
//      // if (LinkedListKV::ContainsKey<String::String, double>(mSigmas, newName, String::AreEqual)) {
//      if (mSigmas.ContainsKey(newName)) {
//         printf("A component named %s already exists.", newName.Get());
//         return;
//      }
//      double val = mSigmas.Remove(oldName);
//      if (val != NULL) {
//         mSigmas.Put(newName, val);
//      }
//   }
//
//   __device__ UncertainValue2 add(UncertainValue2 uvs[], int uvsLen)
//   {
//      UncertainValue2::KeySet keys;
//      double sum = 0.0;
//      for (int k = 0; k < uvsLen; ++k) {
//         sum += uvs[k].doubleValue();
//         auto sigmasKeys = uvs[k].getComponents().GetKeys();
//         keys.Add(sigmasKeys);
//      }
//      UncertainValue2 res(sum);
//
//      UncertainValue2::KeySetItr itr(keys);
//
//      while (itr.HasNext()) {
//         auto src = itr.GetValue();
//         itr.Next();
//         double unc = 0.0;
//         // This seems right but is it????
//         for (int k = 0; k < uvsLen; ++k) {
//            auto uv = uvs[k];
//            unc += uv.getComponent(src) * copysign(1.0, uv.doubleValue());
//         }
//         res.assignComponent(src, unc);
//      }
//      return res;
//   }
//
//   __device__ UncertainValue2 add(double a, UncertainValue2& uva, double b, UncertainValue2& uvb)
//   {
//      UncertainValue2::KeySet keys;
//      UncertainValue2 res(a * uva.doubleValue() + b * uvb.doubleValue());
//      auto akeys = uva.getComponents().GetKeys();
//      auto bkeys = uvb.getComponents().GetKeys();
//      keys.Add(akeys);
//      keys.Add(bkeys);
//
//      UncertainValue2::KeySetItr itr(keys);
//      while (itr.HasNext()) {
//         String::String src = itr.GetValue();
//         res.assignComponent(src, a * copysign(1.0, uva.doubleValue()) * uva.getComponent(src) + b * copysign(1.0, uvb.doubleValue()) * uvb.getComponent(src));
//         itr.Next();
//      }
//      return res;
//   }
//
//   __device__ UncertainValue2 subtract(UncertainValue2& uva, UncertainValue2& uvb)
//   {
//      return add(1.0, uva, -1.0, uvb);
//   }
//
//   __device__ UncertainValue2 mean(UncertainValue2 uvs[], int uvsLen)
//   {
//      auto uv = add(uvs, uvsLen);
//      return divide(uv, (double)uvsLen);
//   }
//
//   __device__ UncertainValue2 weightedMean(UncertainValue2 cuv[], int uvsLen)
//   {
//      double varSum = 0.0, sum = 0.0;
//
//      for (int k = 0; k < uvsLen; ++k) {
//         auto uv = cuv[k];
//         const double ivar = 1.0 / uv.variance();
//         if (isnan(ivar) || isinf(ivar)) {
//            printf("%s\n", "Unable to compute the weighted mean when one or more datapoints have zero uncertainty.");
//            return NaN();
//         }
//         varSum += ivar;
//         sum += ivar * uv.doubleValue();
//      }
//      const double iVarSum = 1.0 / varSum;
//      if (isnan(iVarSum) || isinf(iVarSum)) {
//         printf("UncertainValue2::weightedMean: badddddd\n");
//         return NaN();
//      }
//      char str[4] = "WM";
//      return UncertainValue2(sum / varSum, str, ::sqrt(1.0 / varSum));
//   }
//
//   __device__ UncertainValue2 uvmin(UncertainValue2 uvs[], int uvsLen)
//   {
//      if (uvsLen == 0) {
//         return NULL;
//      }
//      UncertainValue2 res = uvs[0];
//
//      for (int k = 0; k < uvsLen; ++k) {
//         auto uv = uvs[k];
//         if (uv.doubleValue() < res.doubleValue()) {
//            res = uv;
//         }
//         else if (uv.doubleValue() == res.doubleValue()) {
//            if (uv.uncertainty() > res.uncertainty()) {
//               res = uv;
//            }
//         }
//      }
//      return res;
//   }
//
//   __device__ UncertainValue2 uvmax(UncertainValue2 uvs[], int uvsLen)
//   {
//      if (uvs == 0) {
//         return NULL;
//      }
//      UncertainValue2 res = uvs[0];
//
//      for (int k = 0; k < uvsLen; ++k) {
//         auto uv = uvs[k];
//         if (uv.doubleValue() > res.doubleValue()) {
//            res = uv;
//         }
//         else if (uv.doubleValue() == res.doubleValue()) {
//            if (uv.uncertainty() > res.uncertainty()) {
//               res = uv;
//            }
//         }
//      }
//      return res;
//   }
//
//   __device__ UncertainValue2 add(UncertainValue2& v1, double v2)
//   {
//      return UncertainValue2(v1.doubleValue() + v2, v1.getComponents());
//   }
//
//   __device__ UncertainValue2 add(double v1, UncertainValue2& v2)
//   {
//      return UncertainValue2(v2.doubleValue() + v1, v2.getComponents());
//   }
//
//   __device__ UncertainValue2 add(UncertainValue2& v1, UncertainValue2& v2)
//   {
//      return add(1.0, v1, 1.0, v2);
//   }
//
//   __device__ UncertainValue2 multiply(double v1, UncertainValue2& v2)
//   {
//      if (v2.uncertainty() < 0.0) {
//         printf("Error: v2.uncertainty() < 0.0");
//         return NULL;
//      }
//      UncertainValue2 res(v1 * v2.doubleValue());
//
//      auto m = v2.getComponents();
//      UncertainValue2::ComponentMapItr itr(m);
//      while (itr.HasNext()) {
//         auto k = itr.GetKey();
//         res.assignComponent(k, v1 * itr.GetValue());
//         itr.Next();
//      }
//      return res;
//   }
//
//   __device__ UncertainValue2 multiply(UncertainValue2& v1, UncertainValue2& v2)
//   {
//      UncertainValue2::KeySet keys;
//      auto v1keys = v1.getComponents().GetKeys();
//      auto v2keys = v2.getComponents().GetKeys();
//      keys.Add(v1keys);
//      keys.Add(v2keys);
//
//      UncertainValue2 res(v1.doubleValue() * v2.doubleValue());
//      UncertainValue2::KeySetItr itr(keys);
//      while (itr.HasNext()) {
//         auto src = itr.GetValue();
//         res.assignComponent(src, v1.doubleValue() * v2.getComponent(src) + v2.doubleValue() * v1.getComponent(src));
//         itr.Next();
//      }
//
//      return res;
//   }
//
//   __device__ UncertainValue2 invert(UncertainValue2& v)
//   {
//      return divide(1.0, v);
//   }
//
//   __device__ UncertainValue2 divide(UncertainValue2& v1, UncertainValue2& v2)
//   {
//      UncertainValue2 res(v1.doubleValue() / v2.doubleValue());
//      if (!(isnan(res.doubleValue()) || isinf(res.doubleValue()))) {
//         UncertainValue2::KeySet keys;
//         auto v1keys = v1.getComponents().GetKeys();
//         auto v2keys = v2.getComponents().GetKeys();
//         keys.Add(v1keys);
//         keys.Add(v2keys);
//
//         const double ua = fabs(1.0 / v2.doubleValue());
//         const double ub = fabs(v1.doubleValue() / (v2.doubleValue() * v2.doubleValue()));
//
//         UncertainValue2::KeySetItr itr(keys);
//         while (itr.HasNext()) {
//            auto src = itr.GetValue();
//            res.assignComponent(src, ua * v1.getComponent(src) + ub * v2.getComponent(src));
//            itr.Next();
//         }
//      }
//      return res;
//   }
//
//   __device__ UncertainValue2 divide(double a, UncertainValue2& b)
//   {
//      UncertainValue2 res(a / b.doubleValue());
//      if (!(isnan(res.doubleValue()) || isinf(res.doubleValue()))) {
//         const double ub = fabs(a / (b.doubleValue() * b.doubleValue()));
//
//         auto m = b.getComponents();
//         UncertainValue2::ComponentMapItr itr(m);
//         while (itr.HasNext()) {
//            res.assignComponent(itr.GetKey(), ub * itr.GetValue());
//            itr.Next();
//         }
//      }
//      return res;
//   }
//
//   __device__ UncertainValue2 divide(UncertainValue2& a, double b)
//   {
//      if (isnan(1.0 / b)) {
//         return UncertainValue2(CUDART_NAN);
//      }
//      if (isinf(1.0 / b)) {
//         return UncertainValue2(CUDART_INF);
//      }
//      UncertainValue2 res(a.doubleValue() / b);
//      const double ua = fabs(1.0 / b);
//
//      auto m = a.getComponents();
//      UncertainValue2::ComponentMapItr itr(m);
//      while (itr.HasNext()) {
//         res.assignComponent(itr.GetKey(), ua * itr.GetValue());
//         itr.Next();
//      }
//      return res;
//   }
//
//   __device__ UncertainValue2 exp(UncertainValue2& x)
//   {
//      if (isnan(x.doubleValue()) || isinf(x.doubleValue())) {
//         printf("exp: invalid value\n");
//         return NULL;
//      }
//
//      double ex = ::exp(x.doubleValue());
//      UncertainValue2 res(ex);
//
//      auto m = x.getComponents();
//      UncertainValue2::ComponentMapItr itr(m);
//      while (itr.HasNext()) {
//         res.assignComponent(itr.GetKey(), ex * itr.GetValue());
//         itr.Next();
//      }
//      return res;
//   }
//
//   __device__ UncertainValue2 log(UncertainValue2& v2)
//   {
//      double tmp = 1.0 / v2.doubleValue();
//      const double lv = ::log(v2.doubleValue());
//      if (isnan(tmp) || isnan(lv)) {
//         return UncertainValue2(CUDART_NAN);
//      }
//      if (isinf(tmp) || isinf(lv)) {
//         return UncertainValue2(CUDART_INF);
//      }
//      UncertainValue2 res(lv);
//
//      auto m = v2.getComponents();
//      UncertainValue2::ComponentMapItr itr(m);
//      while (itr.HasNext()) {
//         res.assignComponent(itr.GetKey(), tmp * itr.GetValue());
//         itr.Next();
//      }
//      return res;
//   }
//
//   __device__ UncertainValue2 pow(UncertainValue2& v1, double n)
//   {
//      if (v1.doubleValue() == 0.0) {
//         return UncertainValue2(0.0);
//      }
//      const double f = ::pow(v1.doubleValue(), n);
//      const double df = n * ::pow(v1.doubleValue(), n - 1.0);
//      UncertainValue2 res(f);
//
//      auto m = v1.getComponents();
//      UncertainValue2::ComponentMapItr itr(m);
//      while (itr.HasNext()) {
//         res.assignComponent(itr.GetKey(), itr.GetValue() * df);
//         itr.Next();
//      }
//
//      return res;
//   }
//
//   __device__ UncertainValue2 UncertainValue2::sqrt()
//   {
//      return pow(*this, 0.5);
//   }
//
//   __device__ UncertainValue2 sqrt(UncertainValue2& uv)
//   {
//      return pow(uv, 0.5);
//   }
//
//   __device__ LinkedList::Node<UncertainValue2>* quadratic(UncertainValue2& a, UncertainValue2& b, UncertainValue2& c)
//   {
//      // q=-0.5*(b+signum(b)*sqrt(pow(b,2.0)-4*a*c))
//      // return [ q/a, c/q ]
//      auto uv0 = pow(b, 2.0);
//      auto uv1 = multiply(a, c);
//      UncertainValue2 r = add(1.0, uv0, -4.0, uv1);
//      if (r.doubleValue() <= 0.0) {
//         return NULL;
//      }
//      auto uv2 = r.sqrt();
//      auto uv3 = multiply(copysign(1.0, b.doubleValue()), uv2);
//      auto uv4 = add(b, uv3);
//      UncertainValue2 q = multiply(-0.5, uv4);
//      LinkedList::Node<UncertainValue2>* head = NULL;
//      auto uv5 = divide(q, a);
//      auto uv6 = divide(c, q);
//      LinkedList::InsertHead(&head, uv5);
//      LinkedList::InsertHead(&head, uv6);
//      return head;
//   }
//
//   __device__ double UncertainValue2::doubleValue() const
//   {
//      return mValue;
//   }
//
//   __device__ bool UncertainValue2::isUncertain()
//   {
//      return !mSigmas.IsEmpty();
//   }
//
//   __device__ double UncertainValue2::uncertainty()
//   {
//      return ::sqrt(variance());
//   }
//
//   __device__ double UncertainValue2::variance()
//   {
//      return mSigmas.Aggregate([](double a) { return a*a; });
//   }
//
//   __device__ double UncertainValue2::fractionalUncertainty()
//   {
//      if (isnan(1.0 / mValue)) {
//         return CUDART_NAN;
//      }
//      if (isinf(1.0 / mValue)) {
//         return CUDART_INF;
//      }
//      return fabs(uncertainty() / mValue);
//   }
//
//   __device__ bool UncertainValue2::operator==(UncertainValue2& other)
//   {
//      if (this == &other) {
//         return true;
//      }
//      return mSigmas == other.mSigmas && mValue == other.doubleValue();
//      //return LinkedListKV::AreEquivalentSets<String::String, double>(mSigmas, other.getComponents(), String::AreEqual, [](double a, double b) { return a == b; }) && (mValue == other.doubleValue());
//   }
//
//   __device__ bool UncertainValue2::equals(UncertainValue2& other)
//   {
//      return *this == other;
//   }
//
//   __device__ int UncertainValue2::compareTo(UncertainValue2& o)
//   {
//      if (&o == this) return 0;
//      return (mValue == o.mValue) && (uncertainty() == o.uncertainty());
//   }
//
//   __device__ bool UncertainValue2::lessThan(UncertainValue2& uv2)
//   {
//      return mValue < uv2.mValue;
//   }
//
//   __device__ bool UncertainValue2::greaterThan(UncertainValue2& uv2)
//   {
//      return mValue > uv2.mValue;
//   }
//
//   __device__ bool UncertainValue2::lessThanOrEqual(UncertainValue2& uv2)
//   {
//      return mValue <= uv2.mValue;
//   }
//
//   __device__ bool UncertainValue2::greaterThanOrEqual(UncertainValue2& uv2)
//   {
//      return mValue >= uv2.mValue;
//   }
//
//   __device__ UncertainValue2 sqr(UncertainValue2& uv)
//   {
//      return pow(uv, 2.0);
//   }
//
//   __device__ UncertainValue2 negate(UncertainValue2& uv)
//   {
//      return UncertainValue2(-uv.doubleValue(), uv.getComponents());
//   }
//
//   __device__ UncertainValue2 atan(UncertainValue2& uv)
//   {
//      double f = ::atan(uv.doubleValue());
//      double df = 1.0 / (1.0 + uv.doubleValue() * uv.doubleValue());
//
//      if (isnan(f)) {
//         return UncertainValue2(CUDART_NAN);
//      }
//      if (isinf(df)) {
//         return UncertainValue2(CUDART_INF);
//      }
//      UncertainValue2 res(f);
//
//      auto m = uv.getComponents();
//      UncertainValue2::ComponentMapItr itr(m);
//      while (itr.HasNext()) {
//         res.assignComponent(itr.GetKey(), df * itr.GetValue());
//         itr.Next();
//      }
//      return res;
//   }
//
//   __device__ UncertainValue2 atan2(UncertainValue2& y, UncertainValue2& x)
//   {
//      double f = ::atan2(y.doubleValue(), x.doubleValue());
//      double df = 1.0 / (1.0 + ::pow(y.doubleValue() / x.doubleValue(), 2.0));
//
//      if (isnan(f)) {
//         return UncertainValue2(CUDART_NAN);
//      }
//      if (isinf(df)) {
//         return UncertainValue2(CUDART_INF);
//      }
//      UncertainValue2 res(f);
//
//      auto comp = divide(y, x).getComponents();
//      UncertainValue2::ComponentMapItr itr(comp);
//      while (itr.HasNext()) {
//         res.assignComponent(itr.GetKey(), df * itr.GetValue());
//         itr.Next();
//      }
//      return res;
//   }
//
//   __device__ UncertainValue2 positiveDefinite(UncertainValue2& uv)
//   {
//      UncertainValue2 ret(0.0, uv.getComponents());
//      return uv.doubleValue() >= 0.0 ? uv : ret;
//   }
//
//   __device__ Key::Key(String::String src1, String::String src2)
//   {
//      mSource1 = src1;
//      mSource2 = src2;
//   }
//
//   __device__ bool Key::operator==(const Key& k2) const
//   {
//      return (mSource1 == k2.mSource1 && mSource2 == k2.mSource2) || (mSource1 == k2.mSource2 && mSource2 == k2.mSource1);
//   }
//
//   __device__ unsigned int Key::Hashcode()
//   {
//      return mSource1.HashCode() + mSource2.HashCode();
//   }
//
//   //__device__ bool Key::AreEqual(Key& k1, Key& k2)
//   //{
//   //   if (&k1 == &k2) return true;
//   //   return k1 == k2;
//   //}
//
//   __device__ Correlations::Correlations()
//   {
//   }
//
//   __device__ void Correlations::add(String::String& src1, String::String& src2, double corr)
//   {
//      if (!(corr >= -1.0) && (corr <= 1.0)) {
//         printf("%s\n", "Correlations::add: invalid bound");
//         return;
//      }
//      corr = ::fmax(corr, 1.0);
//      corr = ::fmin(corr, -1.0);
//      Key k = Key(src1, src2);
//      mCorrelations.Put(k, corr);
//   }
//
//   __device__ double Correlations::get(String::String& src1, String::String& src2)
//   {
//      double r;
//      auto k = Key(src1, src2);
//      return mCorrelations.GetValue(k, r) ? r : 0.0;
//   }
//
//   __device__ double UncertainValue2::variance(Correlations& corr)
//   {
//      UncertainValue2::ComponentMap sigmas = getComponents();
//
//      UncertainValue2::KeySet keys;
//      auto sigkeys = sigmas.GetKeys();
//      keys.Add(sigkeys);
//
//      UncertainValue2::KeySetItr itr1(keys);
//      double res = 0.0;
//      while (itr1.HasNext()) {
//         String::String key = itr1.GetValue();
//         double val;
//         if (!sigmas.GetValue(key, val)) {
//            printf("UncertainValue2::variance: key %s not found.", key.Get());
//         }
//         res += val * val;
//         itr1.Next();
//      }
//
//      UncertainValue2::KeySetItr itr2(keys);
//
//      while (itr1.HasNext()) {
//         itr2 = itr1;
//         itr2.Next();
//         while (itr2.HasNext()) {
//            auto key1 = itr1.GetValue();
//            auto key2 = itr2.GetValue();
//            double sigma1, sigma2;
//            if (!sigmas.GetValue(key1, sigma1)) {
//               printf("UncertainValue2::variance: key1 %s not found.", key1.Get());
//            }
//            if (!sigmas.GetValue(key2, sigma2)) {
//               printf("UncertainValue2::variance: key2 %s not found.", key2.Get());
//            }
//            res += 2.0 * sigma1 * sigma2 * corr.get(key1, key2);
//            itr2.Next();
//         }
//         itr1.Next();
//      }
//      return res;
//   }
//
//   __device__ double UncertainValue2::uncertainty(Correlations& corr)
//   {
//      return ::sqrt(variance(corr));
//   }
//
//   __device__ void UncertainValue2::PrintSigmas()
//   {
//      UncertainValue2::ComponentMapItr itr(mSigmas);
//      while (itr.HasNext()) {
//         printf("%s: %lf\n", itr.GetKey().Get(), itr.GetValue());
//         itr.Next();
//      }
//   }
//
//   __device__ bool AreEqual(UncertainValue2& a, UncertainValue2& b)
//   {
//      return a == b;
//   }
//}
