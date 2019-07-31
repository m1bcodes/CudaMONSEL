#include "UncertainValue2.cuh"

#include <stdio.h>
#include <math.h>
#include <algorithm>

#include "Amphibian\Hasher.cuh"
#include "Amphibian\String.cuh"

namespace UncertainValue2
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const char DEFAULT[] = "Default";
   __device__ int sDefIndex = 0;

   __constant__ const long long serialVersionUID = 119495064970078787L;
   __constant__ const int MAX_LEN = 32;
#else
   const char DEFAULT[] = "Default";
   int sDefIndex = 0;

   const long long serialVersionUID = 119495064970078787L;
   const int MAX_LEN = 32;
#endif

   __host__ __device__ UncertainValue2::UncertainValue2()
   {
   }

   __host__ __device__ UncertainValue2::UncertainValue2(const data_type v, const data_type dv) : mValue(v)
   {
      char tmpName[MAX_LEN];
      memcpy(tmpName, DEFAULT, sizeof(DEFAULT));
      //itoa(sDefIndex++, tmpName + sizeof(DEFAULT), MAX_LEN - sizeof(DEFAULT));
      amp::IToA(sDefIndex++, tmpName + sizeof(DEFAULT), MAX_LEN - sizeof(DEFAULT)); // eg Default1, Default2 etc
      assignComponent(tmpName, dv);
   }

   __host__ __device__ UncertainValue2::UncertainValue2(const data_type v) : mValue(v)
   {
   }

   UncertainValue2::UncertainValue2(data_type v, const char source[], data_type dv) : mValue(v)
   {
      assignComponent(source, dv);
   }

   __host__ __device__ UncertainValue2::UncertainValue2(data_type v, const ComponentMapT& sigmas) : mValue(v), mSigmas(sigmas)
   {
   }

   __host__ __device__ UncertainValue2::UncertainValue2(const UncertainValue2& other) : mValue(other.doubleValue()), mSigmas(other.mSigmas)
   {
   }

   __host__ __device__ UncertainValue2 ONE()
   {
      return UncertainValue2(1.0);
   }

   __host__ __device__ UncertainValue2 NaN()
   {
      return UncertainValue2(NAN);
   }

   UncertainValue2 POSITIVE_INFINITY()
   {
      return UncertainValue2(INFINITY);
   }

   UncertainValue2 NEGATIVE_INFINITY()
   {
      return UncertainValue2(-INFINITY);
   }

   __host__ __device__ UncertainValue2 ZERO()
   {
      return UncertainValue2(0.0);
   }

   __host__ __device__ UncertainValue2& UncertainValue2::operator=(const UncertainValue2& other)
   {
      mValue = other.doubleValue();
      mSigmas = other.mSigmas;

      return *this;
   }

   __host__ __device__ unsigned int UncertainValue2::hashCode() const
   {
      // https://docs.oracle.com/javase/8/docs/api/java/util/Arrays.html#hashCode-java.lang.Object:A-
      unsigned int res = 1;
      static const unsigned int PRIME = 31; // https://docs.oracle.com/javase/8/docs/api/java/util/List.html#hashCode--
      //auto khashfcn = mSigmas.hash_function();
      amp::string_hash khashfcn;
      Hasher::DoubleHashFcn vhashfcn;
      res += PRIME * res + (long long)((long long)mValue ^ ((long long)mValue >> 32)); // https://docs.oracle.com/javase/7/docs/api/java/lang/Double.html#hashCode()
      res += PRIME * res + mSigmas.hashCode();
      return res;
   }

   __host__ __device__ void UncertainValue2::assignComponent(const StringT& name, const data_type sigma)
   {
      if (sigma != 0.0) {
         mSigmas.insert(amp::make_pair(name, sigma));
      }
      else {
         mSigmas.erase(name);
      }
   }

   __host__ __device__ data_type UncertainValue2::getComponent(const StringT& src) const
   {
      auto itr = mSigmas.find(src);
      return itr == mSigmas.end() ? 0 : (data_type)itr->second;
   }

   UncertainValue2::ComponentMapT& UncertainValue2::getComponents()
   {
      return mSigmas;
   }

   __host__ __device__ const UncertainValue2::ComponentMapT& UncertainValue2::getComponents() const
   {
      return mSigmas;
   }

   __host__ __device__ UncertainValue2::ComponentMapT::const_iterator UncertainValue2::getComponentsItrBegin() const
   {
      return mSigmas.cbegin();
   }

   __host__ __device__ UncertainValue2::ComponentMapT::const_iterator UncertainValue2::getComponentsItrEnd() const
   {
      return mSigmas.cend();
   }

   bool UncertainValue2::hasComponent(const StringT& src) const
   {
      return getComponent(src) != 0.0;
   }

   void UncertainValue2::renameComponent(const StringT& oldName, const StringT& newName)
   {
      if (mSigmas.find(newName) != mSigmas.end()) {
         printf("A component named %s already exists.", newName.c_str());
         return;
      }
      data_type val = mSigmas.erase(oldName);
      if (val != NULL) {
         mSigmas.insert(amp::make_pair(newName, val));
      }
   }

   UncertainValue2 add(const UncertainValue2 uvs[], int uvsLen)
   {
      UncertainValue2::KeySetT keys;
      data_type sum = 0.0;
      for (int k = 0; k < uvsLen; ++k) {
         sum += uvs[k].doubleValue();

         for (auto itr = uvs[k].getComponentsItrBegin(); itr != uvs[k].getComponentsItrEnd(); ++itr) {
            keys.insert(itr->first);
         }
      }
      UncertainValue2 res(sum);

      for (auto itr = keys.begin(); itr != keys.end(); ++itr) {
         auto src = *itr;
         data_type unc = 0.0;
         // This seems right but is it????
         for (int k = 0; k < uvsLen; ++k) {
            auto uv = uvs[k];
            unc += uv.getComponent(src) * copysign(1.0, uv.doubleValue());
         }
         res.assignComponent(src, unc);
      }
      return res;
   }

   __host__ __device__ UncertainValue2 add(const data_type a, const UncertainValue2& uva, const data_type b, const UncertainValue2& uvb)
   {
      UncertainValue2::KeySetT keys;
      UncertainValue2 res(a * uva.doubleValue() + b * uvb.doubleValue());

      for (auto itr = uva.getComponentsItrBegin(); itr != uva.getComponentsItrEnd(); ++itr) {
         keys.insert(itr->first);
      }

      for (auto itr = uvb.getComponentsItrBegin(); itr != uvb.getComponentsItrEnd(); ++itr) {
         keys.insert(itr->first);
      }

      for (auto itr = keys.begin(); itr != keys.end(); ++itr) {
         StringT src = *itr;
         res.assignComponent(src, a * copysign(1.0, uva.doubleValue()) * uva.getComponent(src) + b * copysign(1.0, uvb.doubleValue()) * uvb.getComponent(src));
      }
      return res;
   }

   UncertainValue2 subtract(const UncertainValue2& uva, const UncertainValue2& uvb)
   {
      return add(1.0, uva, -1.0, uvb);
   }

   UncertainValue2 mean(const UncertainValue2 uvs[], int uvsLen)
   {
      auto uv = add(uvs, uvsLen);
      return divide(uv, (data_type)uvsLen);
   }

   UncertainValue2 weightedMean(const UncertainValue2 cuv[], int uvsLen)
   {
      data_type varSum = 0.0, sum = 0.0;

      for (int k = 0; k < uvsLen; ++k) {
         auto uv = cuv[k];
         const data_type ivar = 1.0 / uv.variance();
         if (isnan(ivar) || isinf(ivar)) {
            printf("%s\n", "Unable to compute the weighted mean when one or more datapoints have zero uncertainty.");
            return NaN();
         }
         varSum += ivar;
         sum += ivar * uv.doubleValue();
      }
      const data_type iVarSum = 1.0 / varSum;
      if (isnan(iVarSum) || isinf(iVarSum)) {
         printf("UncertainValue2::weightedMean: badddddd\n");
         return NaN();
      }
      char str[4] = "WM";
      return UncertainValue2(sum / varSum, str, ::sqrt(1.0 / varSum));
   }

   UncertainValue2 uvmin(const UncertainValue2 uvs[], int uvsLen)
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

   UncertainValue2 uvmax(const UncertainValue2 uvs[], int uvsLen)
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

   UncertainValue2 add(const UncertainValue2& v1, data_type v2)
   {
      auto tmp = UncertainValue2(v2);
      return add(1.0, v1, 1.0, tmp);
   }

   UncertainValue2 add(data_type v1, const UncertainValue2& v2)
   {
      auto tmp = UncertainValue2(v1);
      return add(1.0, tmp, 1.0, v2);
   }

   __host__ __device__ UncertainValue2 add(const UncertainValue2& v1, const UncertainValue2& v2)
   {
      return add(1.0, v1, 1.0, v2);
   }

   __host__ __device__ UncertainValue2 multiply(data_type v1, const UncertainValue2& v2)
   {
      if (v2.uncertainty() < 0.0) {
         printf("UncertainValue2 multiply: v2.uncertainty() = %.10e < 0.0", v2.uncertainty());
         return NULL;
      }
      UncertainValue2 res(v1 * v2.doubleValue());

      for (auto itr = v2.getComponentsItrBegin(); itr != v2.getComponentsItrEnd(); ++itr) {
         res.assignComponent(itr->first, v1 * (data_type)itr->second);
      }
      return res;
   }

   UncertainValue2 multiply(const UncertainValue2& v1, const UncertainValue2& v2)
   {
      UncertainValue2::KeySetT keys;
      for (auto itr = v1.getComponentsItrBegin(); itr != v1.getComponentsItrEnd(); ++itr) {
         keys.insert(itr->first);
      }
      for (auto itr = v2.getComponentsItrBegin(); itr != v2.getComponentsItrEnd(); ++itr) {
         keys.insert(itr->first);
      }

      UncertainValue2 res(v1.doubleValue() * v2.doubleValue());
      for (auto src : keys) {
         //printf("v1.doubleValue: %.10e, ", v1.doubleValue());
         //printf("v2.getComponent %s: %.10e, ", src.c_str(), v2.getComponent(src));
         //printf("v2.doubleValue: %.10e, ", v2.doubleValue());
         //printf("v1.getComponent %s: %.10e, ", src.c_str(), v1.getComponent(src));
         //printf("\n");
         res.assignComponent(src, v1.doubleValue() * v2.getComponent(src) + v2.doubleValue() * v1.getComponent(src));
      }

      return res;
   }

   UncertainValue2 invert(const UncertainValue2& v)
   {
      return divide(1.0, v);
   }

   __host__ __device__ UncertainValue2 divide(const UncertainValue2& v1, const UncertainValue2& v2)
   {
      const data_type v = v1.doubleValue() / v2.doubleValue();
      if (isnan(v) || isinf(v)) {
         printf("UncertainValue2::divide: isnan(v) || isinf(v) (%.10e)\n");
         return v;
      }

      UncertainValue2::KeySetT keys;

      //std::transform(v1.getComponentsItrBegin(), v1.getComponentsItrEnd(), std::inserter(keys, keys.end()), [](LinkedListKV::Node<StringT, data_type> pair){ return pair.first; });
      //std::transform(v2.getComponentsItrBegin(), v2.getComponentsItrEnd(), std::inserter(keys, keys.end()), [](LinkedListKV::Node<StringT, data_type> pair){ return pair.first; });

      for (auto itr = v1.getComponentsItrBegin(); itr != v1.getComponentsItrEnd(); ++itr) {
         keys.insert(itr->first);
      }
      for (auto itr = v2.getComponentsItrBegin(); itr != v2.getComponentsItrEnd(); ++itr) {
         keys.insert(itr->first);
      }

      const data_type ua = ::fabs(1.0 / v2.doubleValue());
      const data_type ub = ::fabs(v1.doubleValue() / (v2.doubleValue() * v2.doubleValue()));

      UncertainValue2 res(v);
      for (auto itr = keys.begin(); itr != keys.end(); ++itr) {
         const auto &src = *itr;
         res.assignComponent(src, ua * v1.getComponent(src) + ub * v2.getComponent(src));
      }

      //std::transform(keys.begin(), keys.end(), std::inserter(res.getComponents(), res.getComponents().end()), [&](StringT name){ return amp::make_pair(name, ua * v1.getComponent(name) + ub * v2.getComponent(name)); });
      return res;
   }

   __host__ __device__ void divide(const UncertainValue2& v1, const UncertainValue2& v2, UncertainValue2& res)
   {
      const data_type v = res.doubleValue(); // v1.doubleValue() / v2.doubleValue();
      if (isnan(v) || isinf(v)) {
         printf("UncertainValue2::divide: isnan(v) || isinf(v) (%.10e)\n");
      }

      UncertainValue2::KeySetT keys;

      //std::transform(v1.getComponentsItrBegin(), v1.getComponentsItrEnd(), std::inserter(keys, keys.end()), [](LinkedListKV::Node<StringT, data_type> pair){ return pair.first; });
      //std::transform(v2.getComponentsItrBegin(), v2.getComponentsItrEnd(), std::inserter(keys, keys.end()), [](LinkedListKV::Node<StringT, data_type> pair){ return pair.first; });

      for (auto itr = v1.getComponentsItrBegin(); itr != v1.getComponentsItrEnd(); ++itr) {
         keys.insert(itr->first);
      }
      for (auto itr = v2.getComponentsItrBegin(); itr != v2.getComponentsItrEnd(); ++itr) {
         keys.insert(itr->first);
      }

      const data_type ua = fabs(1.0 / v2.doubleValue());
      const data_type ub = fabs(v1.doubleValue() / (v2.doubleValue() * v2.doubleValue()));

      for (auto itr = keys.begin(); itr != keys.end(); ++itr) {
         const auto &src = *itr;
         res.assignComponent(src, ua * v1.getComponent(src) + ub * v2.getComponent(src));
      }
      //std::transform(keys.begin(), keys.end(), std::inserter(res.getComponents(), res.getComponents().end()), [&](StringT name){ return amp::make_pair(name, ua * v1.getComponent(name) + ub * v2.getComponent(name)); });
   }

   UncertainValue2 divide(data_type a, const UncertainValue2& b)
   {
      UncertainValue2 res(a / b.doubleValue());
      if (!(isnan(res.doubleValue()) || isinf(res.doubleValue()))) {
         const data_type ub = fabs(a / (b.doubleValue() * b.doubleValue()));

         for (auto itr = b.getComponentsItrBegin(); itr != b.getComponentsItrEnd(); ++itr) {
            res.assignComponent(itr->first, ub * (data_type)itr->second);
         }
      }
      return res;
   }

   UncertainValue2 divide(const UncertainValue2& a, data_type b)
   {
      if (isnan(1.0 / b)) {
         return UncertainValue2(NAN);
      }
      if (isinf(1.0 / b)) {
         return UncertainValue2(INFINITY);
      }
      UncertainValue2 res(a.doubleValue() / b);
      const data_type ua = fabs(1.0 / b);

      for (auto itr = a.getComponentsItrBegin(); itr != a.getComponentsItrEnd(); ++itr) {
         res.assignComponent(itr->first, ua * (data_type)itr->second);
      }
      return res;
   }

   UncertainValue2 exp(const UncertainValue2& x)
   {
      if (isnan(x.doubleValue()) || isinf(x.doubleValue())) {
         printf("exp: invalid value\n");
         return NULL;
      }

      data_type ex = ::exp(x.doubleValue());
      UncertainValue2 res(ex);

      for (auto itr = x.getComponentsItrBegin(); itr != x.getComponentsItrEnd(); ++itr) {
         res.assignComponent(itr->first, ex * (data_type)itr->second);
      }
      return res;
   }

   UncertainValue2 log(const UncertainValue2& v2)
   {
      data_type tmp = 1.0 / v2.doubleValue();
      const data_type lv = ::log(v2.doubleValue());
      if (isnan(tmp) || isnan(lv)) {
         return UncertainValue2(NAN);
      }
      if (isinf(tmp) || isinf(lv)) {
         return UncertainValue2(INFINITY);
      }
      UncertainValue2 res(lv);

      for (auto itr = v2.getComponentsItrBegin(); itr != v2.getComponentsItrEnd(); ++itr) {
         res.assignComponent(itr->first, tmp * (data_type)itr->second);
      }

      return res;
   }

   UncertainValue2 pow(const UncertainValue2& v1, data_type n)
   {
      if (v1.doubleValue() == 0.0) {
         return UncertainValue2(0.0);
      }
      const data_type f = ::pow(v1.doubleValue(), n);
      const data_type df = n * ::pow(v1.doubleValue(), n - 1.0);
      UncertainValue2 res(f);

      for (auto itr = v1.getComponentsItrBegin(); itr != v1.getComponentsItrEnd(); ++itr) {
         res.assignComponent(itr->first, df * (data_type)itr->second);
      }

      return res;
   }

   UncertainValue2 UncertainValue2::sqrt() const
   {
      return pow(*this, 0.5);
   }

   UncertainValue2 sqrt(const UncertainValue2& uv)
   {
      return pow(uv, 0.5);
   }

   //UncertainValue2::ResultT quadratic(const UncertainValue2& a, const UncertainValue2& b, const UncertainValue2& c)
   //{
   //   // q=-0.5*(b+signum(b)*sqrt(pow(b,2.0)-4*a*c))
   //   // return [ q/a, c/q ]
   //   auto uv0 = pow(b, 2.0);
   //   auto uv1 = multiply(a, c);
   //   UncertainValue2 r = add(1.0, uv0, -4.0, uv1);
   //   if (r.doubleValue() <= 0.0) {
   //      return UncertainValue2::ResultT();
   //   }
   //   auto uv2 = r.sqrt();
   //   auto uv3 = multiply(copysign(1.0, b.doubleValue()), uv2);
   //   auto uv4 = add(b, uv3);
   //   UncertainValue2 q = multiply(-0.5, uv4);
   //   auto uv5 = divide(q, a);
   //   auto uv6 = divide(c, q);
   //   UncertainValue2::ResultT head;
   //   head.push_back(uv5);
   //   head.push_back(uv6);
   //   return head;
   //}

   __host__ __device__ data_type UncertainValue2::doubleValue() const
   {
      return mValue;
   }

   bool UncertainValue2::isUncertain() const
   {
      return !mSigmas.empty();
   }

   __host__ __device__ data_type UncertainValue2::uncertainty() const
   {
      return ::sqrt(variance());
   }

   __host__ __device__ data_type UncertainValue2::variance() const
   {
      data_type sigma2 = 0.0;
      for (auto &s : mSigmas) {
         //printf("%.10e, ", (data_type)s.second);
         sigma2 += (data_type)s.second * (data_type)s.second;
      }
      //printf("\n");
      return sigma2;
   }

   data_type UncertainValue2::fractionalUncertainty() const
   {
      if (isnan(1.0 / mValue)) {
         return NAN;
      }
      if (isinf(1.0 / mValue)) {
         return INFINITY;
      }
      return ::fabs(uncertainty() / mValue);
   }

   bool UncertainValue2::operator==(const UncertainValue2& other) const
   {
      if (this == &other) {
         return true;
      }

      for (auto s : other.mSigmas) {
         auto itr = mSigmas.find(s.first);
         if (itr == mSigmas.end()) return false;
         if ((data_type)itr->second != (data_type)s.second) return false;
      }

      return mValue == other.doubleValue();
   }

   bool UncertainValue2::equals(const UncertainValue2& other) const
   {
      return *this == other;
   }

   int UncertainValue2::compareTo(const UncertainValue2& o)
   {
      if (&o == this) return 0;
      return (mValue == o.mValue) && (uncertainty() == o.uncertainty());
   }

   bool UncertainValue2::lessThan(const UncertainValue2& uv2)
   {
      return mValue < uv2.mValue;
   }

   bool UncertainValue2::greaterThan(const UncertainValue2& uv2)
   {
      return mValue > uv2.mValue;
   }

   bool UncertainValue2::lessThanOrEqual(const UncertainValue2& uv2)
   {
      return mValue <= uv2.mValue;
   }

   bool UncertainValue2::greaterThanOrEqual(const UncertainValue2& uv2)
   {
      return mValue >= uv2.mValue;
   }

   UncertainValue2 sqr(const UncertainValue2& uv)
   {
      return pow(uv, 2.0);
   }

   UncertainValue2 negate(const UncertainValue2& uv)
   {
      return add(-1.0, uv);
   }

   UncertainValue2 atan(const UncertainValue2& uv)
   {
      data_type f = ::atan(uv.doubleValue());
      data_type df = 1.0 / (1.0 + uv.doubleValue() * uv.doubleValue());

      if (isnan(f)) {
         return UncertainValue2(NAN);
      }
      if (isinf(df)) {
         return UncertainValue2(INFINITY);
      }
      UncertainValue2 res(f);

      for (auto itr = uv.getComponentsItrBegin(); itr != uv.getComponentsItrEnd(); ++itr) {
         res.assignComponent(itr->first, df * (data_type)itr->second);
      }

      return res;
   }

   UncertainValue2 atan2(const UncertainValue2& y, const UncertainValue2& x)
   {
      data_type f = ::atan2(y.doubleValue(), x.doubleValue());
      data_type df = 1.0 / (1.0 + ::pow(y.doubleValue() / x.doubleValue(), 2.0));

      if (isnan(f)) {
         return UncertainValue2(NAN);
      }
      if (isinf(df)) {
         return UncertainValue2(INFINITY);
      }
      UncertainValue2 res(f);

      auto m = divide(y, x).getComponents();
      for (auto itr = m.begin(); itr != m.end(); ++itr) {
         res.assignComponent(itr->first, df * (data_type)itr->second);
      }
      return res;
   }

   __host__ __device__ UncertainValue2 positiveDefinite(const UncertainValue2& uv)
   {
      return uv.doubleValue() >= 0.0 ? uv : UncertainValue2(0, uv.getComponents());
   }

   Key::Key(const StringT& src1, const StringT& src2)
   {
      mSource1 = src1;
      mSource2 = src2;
   }

   bool Key::operator==(const Key& k2) const
   {
      return (mSource1 == k2.mSource1 && mSource2 == k2.mSource2) || (mSource1 == k2.mSource2 && mSource2 == k2.mSource1);
   }

   bool Key::operator<(const Key& k2) const
   {
      if (mSource1 != k2.mSource1) {
         return mSource1 < k2.mSource1;
      }
      return mSource2 < k2.mSource2;
   }

   size_t Key::hashCode() const
   {
      //unsigned int s = 0;
      //for (auto ch : mSource1) {
      //   s += ch;
      //}
      //for (auto ch : mSource2) {
      //   s += ch;
      //}
      //return s;

      return mSource1.hashCode() + mSource2.hashCode();;
   }

   Correlations::Correlations()
   {
   }

   void Correlations::add(const StringT& src1, const StringT& src2, data_type corr)
   {
      if (!(corr >= -1.0) && (corr <= 1.0)) {
         printf("%s\n", "Correlations::add: invalid bound");
         return;
      }
      corr = ::fmax(corr, 1.0);
      corr = ::fmin(corr, -1.0);
      Key k = Key(src1, src2);
      mCorrelations.insert(std::pair<Key, data_type>(k, corr));
   }

   data_type Correlations::get(const StringT& src1, const StringT& src2) const
   {
      auto k = Key(src1, src2);
      return mCorrelations.at(k);
   }

   data_type UncertainValue2::variance(const Correlations& corr)
   {
      UncertainValue2::ComponentMapT sigmas = getComponents();

      UncertainValue2::KeySetT keys;
      for (auto s : sigmas) {
         keys.insert(s.first);
      }

      UncertainValue2::KeySetT itr1(keys);
      data_type res = 0.0;
      for (auto s : keys) {
         data_type val = sigmas[s];
         if (!val) {
            printf("UncertainValue2::variance: key %s not found.", s.c_str());
         }
         res += val * val;
      }

      for (auto itr1 = keys.begin(); itr1 != keys.end(); ++itr1) {
         for (auto itr2 = itr1; itr2 != keys.end(); ++itr2) {
            auto key1 = *itr1, key2 = *itr2;
            data_type sigma1 = sigmas[key1], sigma2 = sigmas[key2];
            if (!sigma1) {
               printf("UncertainValue2::variance: key1 %s not found.", key1.c_str());
            }
            if (!sigma2) {
               printf("UncertainValue2::variance: key2 %s not found.", key2.c_str());
            }
            res += 2.0 * sigma1 * sigma2 * corr.get(key1, key2);
         }
      }
      return res;
   }

   data_type UncertainValue2::uncertainty(Correlations& corr)
   {
      return ::sqrt(variance(corr));
   }
}
