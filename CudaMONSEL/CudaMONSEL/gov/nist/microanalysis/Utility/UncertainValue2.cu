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

   __device__ UncertainValue2::UncertainValue2(double v, double dv) : mValue(v), mSigmas(NULL)
   {
      char tmpName[MAX_LEN];
      String::IToA(tmpName, atomicAdd(&sDefIndex, 1));
      assignComponent(tmpName, dv);
   }

   __device__ UncertainValue2::UncertainValue2(double v) : mValue(v), mSigmas(NULL)
   {
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

   __device__ UncertainValue2::UncertainValue2(UncertainValue2& other) : mValue(other.doubleValue()), mSigmas(NULL)
   {
      LinkedListKV::DeepCopy<String::String, double>(&mSigmas, other.getComponents());
   }

   __device__ UncertainValue2& UncertainValue2::operator=(UncertainValue2& other)
   {
      mValue = other.doubleValue();
      mSigmas = NULL;
      LinkedListKV::DeepCopy<String::String, double>(&mSigmas, other.getComponents());

      return *this;
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

   __device__ LinkedListKV::Node<String::String, double> * UncertainValue2::getComponents()
   {
      return mSigmas;
   }

   __device__ bool UncertainValue2::hasComponent(String::String src)
   {
      return getComponent(src) != 0.0;
   }

   __device__ void UncertainValue2::renameComponent(String::String oldName, String::String newName)
   {
      if (LinkedListKV::ContainsKey<String::String, double>(mSigmas, newName, String::AreEqual)) {
         printf("A component named %s already exists.", newName.Get());
         return;
      }
      double val = LinkedListKV::Remove<String::String, double>(&mSigmas, oldName, String::AreEqual);
      if (val != NULL) {
         LinkedListKV::InsertHead<String::String, double>(&mSigmas, newName, val);
      }
   }

   __device__ UncertainValue2 add(LinkedList::Node<UncertainValue2>* uvs)
   {
      LinkedList::Node<String::String>* srcs;
      double sum = 0.0;
      auto uvItrHead = uvs;
      while (uvItrHead != NULL) {
         AdvancedLinkedList::AddAllKeys<String::String, double>(&srcs, uvItrHead->GetValue().getComponents(), String::AreEqual);
         sum += uvItrHead->GetValue().doubleValue();
         uvItrHead = uvItrHead->GetNext();
      }
      UncertainValue2 res(sum);
      while (srcs != NULL) {
         auto src = srcs->GetValue();
         srcs = srcs->GetNext();
         double unc = 0.0;
         // This seems right but is it????
         auto uvItrHead = uvs;
         while (uvItrHead != NULL) {
            auto uv = uvItrHead->GetValue();
            unc += copysign(uv.getComponent(src), (uv.doubleValue()));
            uvItrHead = uvItrHead->GetNext();
         }
         res.assignComponent(src, unc);
      }
      return res;
   }

   //__device__ UncertainValue2 add(UncertainValue2* uvs, size_t n)
   //{
   //   LinkedList::Node<UncertainValue2>* head = NULL;
   //   for (int k = 0; k < n; ++k) {
   //      LinkedList::InsertHead(&head, uvs[k]);
   //   }
   //   return add(head);
   //}

   __device__ UncertainValue2 add(double a, UncertainValue2 uva, double b, UncertainValue2 uvb)
   {
      UncertainValue2 res(a * uva.doubleValue() + b * uvb.doubleValue());
      LinkedList::Node<String::String>* srcs = NULL;
      AdvancedLinkedList::AddAllKeys<String::String, double>(&srcs, uva.getComponents(), String::AreEqual);
      AdvancedLinkedList::AddAllKeys<String::String, double>(&srcs, uvb.getComponents(), String::AreEqual);
      while (srcs != NULL) {
         String::String src = srcs->GetValue();
         res.assignComponent(src, a * copysign(1.0, uva.doubleValue()) * uva.getComponent(src) + b * copysign(1.0, uvb.doubleValue()) * uvb.getComponent(src));
         srcs = srcs->GetNext();
      }
      return res;
   }

   __device__ UncertainValue2 subtract(UncertainValue2 uva, UncertainValue2 uvb)
   {
      return add(1.0, uva, -1.0, uvb);
   }

   __device__ UncertainValue2 mean(LinkedList::Node<UncertainValue2>* uvs)
   {
      return divide(add(uvs), (double)LinkedList::Size<UncertainValue2>(uvs));
   }

   __device__ UncertainValue2 weightedMean(LinkedList::Node<UncertainValue2>* cuv)
   {
      double varSum = 0.0, sum = 0.0;
      while (cuv != NULL) {
         auto uv = cuv->GetValue();
         const double ivar = 1.0 / uv.variance();
         if (isnan(ivar) || isinf(ivar)) {
            printf("%s\n", "Unable to compute the weighted mean when one or more datapoints have zero uncertainty.");
            return NULL;
         }
         varSum += ivar;
         sum += ivar * uv.doubleValue();
         cuv = cuv->GetNext();
      }
      const double iVarSum = 1.0 / varSum;
         return (isnan(iVarSum) || isinf(iVarSum)) ? NULL : UncertainValue2(sum / varSum, "WM", ::sqrt(1.0 / varSum));
   }

   __device__ UncertainValue2 min(LinkedList::Node<UncertainValue2>* uvs)
   {
      if (uvs == NULL) {
         return NULL;
      }
      UncertainValue2 res = uvs->GetValue();

      while (uvs != NULL) {
         auto uv = uvs->GetValue();
         if (uv.doubleValue() < res.doubleValue()) {
            res = uv;
         }
         else if (uv.doubleValue() == res.doubleValue()) {
            if (uv.uncertainty() > res.uncertainty()) {
               res = uv;
            }
         }
         uvs = uvs->GetNext();
      }
      return res;
   }

   __device__ UncertainValue2 max(LinkedList::Node<UncertainValue2>* uvs)
   {
      if (uvs == NULL) {
         return NULL;
      }
      UncertainValue2 res = uvs->GetValue();
      while (uvs != NULL) {
         auto uv = uvs->GetValue();
         if (uv.doubleValue() > res.doubleValue()) {
            res = uv;
         }
         else if (uv.doubleValue() == res.doubleValue()) {
            if (uv.uncertainty() > res.uncertainty()) {
               res = uv;
            }
         }
         uvs = uvs->GetNext();
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
      auto srcs = v2.getComponents();
      while (srcs != NULL) {
         res.assignComponent(srcs->GetKey(), v1 * srcs->GetValue());
         srcs = srcs->GetNext();
      }
      return res;
   }

   __device__ UncertainValue2 multiply(UncertainValue2 v1, UncertainValue2 v2)
   {
      UncertainValue2 res(v1.doubleValue() * v2.doubleValue());
      LinkedList::Node<String::String>* srcs = NULL;
      AdvancedLinkedList::AddAllKeys<String::String, double>(&srcs, v1.getComponents(), String::AreEqual);
      AdvancedLinkedList::AddAllKeys<String::String, double>(&srcs, v2.getComponents(), String::AreEqual);
      
      while (srcs != NULL) {
         auto src = srcs->GetValue();
         //printf("%s: ", src.Get());
         //printf("%lf\n", v1.doubleValue() * v2.getComponent(src) + v2.doubleValue() * v1.getComponent(src));
         res.assignComponent(src, v1.doubleValue() * v2.getComponent(src) + v2.doubleValue() * v1.getComponent(src));
         srcs = srcs->GetNext();
      }
      return res;
   }

   __device__ UncertainValue2 invert(UncertainValue2 v)
   {
      return divide(1.0, v);
   }

   __device__ UncertainValue2 divide(UncertainValue2 a, UncertainValue2 b)
   {
      UncertainValue2 res(a.doubleValue() / b.doubleValue());
      if (!(isnan(res.doubleValue()) || isinf(res.doubleValue()))) {
         LinkedList::Node<String::String>* srcs = NULL;
         AdvancedLinkedList::AddAllKeys(&srcs, a.getComponents(), String::AreEqual);
         AdvancedLinkedList::AddAllKeys(&srcs, b.getComponents(), String::AreEqual);
         const double ua = fabs(1.0 / b.doubleValue());
         const double ub = fabs(a.doubleValue() / (b.doubleValue() * b.doubleValue()));

         while (srcs != NULL) {
            auto src = srcs->GetValue();
            res.assignComponent(src, ua * a.getComponent(src) + ub * b.getComponent(src));
            srcs = srcs->GetNext();
         }
      }
      return res;
   }

   __device__ UncertainValue2 divide(double a, UncertainValue2 b)
   {
      UncertainValue2 res(a / b.doubleValue());
      if (!(isnan(res.doubleValue()) || isinf(res.doubleValue()))) {
         const double ub = fabs(a / (b.doubleValue() * b.doubleValue()));
         auto bSigmas = b.getComponents();
         while (bSigmas != NULL) {
            res.assignComponent(bSigmas->GetKey(), ub * bSigmas->GetValue());
            bSigmas = bSigmas->GetNext();
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
      auto sigmaItrHead = a.getComponents();
      while (sigmaItrHead != NULL) {
         res.assignComponent(sigmaItrHead->GetKey(), ua * sigmaItrHead->GetValue());
         sigmaItrHead = sigmaItrHead->GetNext();
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
      auto sigmas = x.getComponents();
      while (sigmas != NULL) {
         res.assignComponent(sigmas->GetKey(), ex * sigmas->GetValue());
         sigmas = sigmas->GetNext();
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
      auto sigmas = v2.getComponents();
      while (sigmas != NULL) {
         res.assignComponent(sigmas->GetKey(), tmp * sigmas->GetValue());
         sigmas = sigmas->GetNext();
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
      auto v1sigmas = v1.getComponents();
      while (v1sigmas != NULL) {
         res.assignComponent(v1sigmas->GetKey(), v1sigmas->GetValue() * df);
         v1sigmas = v1sigmas->GetNext();
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
      return mSigmas != NULL;
   }
   
   __device__ double UncertainValue2::uncertainty()
   {
      return ::sqrt(variance());
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
      if (isnan(1.0 / mValue)) {
         return CUDART_NAN;
      }
      if (isinf(1.0 / mValue)) {
         return CUDART_INF;
      }
      return fabs(uncertainty() / mValue);
   }
   
   __device__ bool UncertainValue2::equals(UncertainValue2* other)
   {
      if (other == NULL) {
         return false;
      }
      if (this == other) {
         return true;
      }
      return LinkedListKV::AreEquivalentSets<String::String, double>(mSigmas, other->getComponents(), String::AreEqual, [](double a, double b) { return a == b; }) && (mValue == other->doubleValue());
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
      LinkedListKV::Node<String::String, double>* sigmas = uv.getComponents();
      while (sigmas != NULL) {
         res.assignComponent(sigmas->GetKey(), df * sigmas->GetValue());
         sigmas = sigmas->GetNext();
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
      LinkedListKV::Node<String::String, double>* sigmas = divide(y, x).getComponents();
      while (sigmas != NULL) {
         res.assignComponent(sigmas->GetKey(), df * sigmas->GetValue());
         sigmas = sigmas->GetNext();
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

   __device__ Correlations::Correlations() : mCorrelations(NULL)
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
      LinkedListKV::InsertHead<Key, double>(&mCorrelations, Key(src1, src2), corr);
   }
   
   __device__ double Correlations::get(String::String src1, String::String src2)
   {
      double r = LinkedListKV::GetValue<Key, double>(mCorrelations, Key(src1, src2), Key::AreEqual);
      return r == NULL ? 0.0 : r;
   }

   __device__ double UncertainValue2::variance(Correlations corr)
   {
      LinkedList::Node<String::String>* keys = NULL;
      LinkedListKV::Node<String::String, double>* sigmas = getComponents();
      AdvancedLinkedList::AddAllKeys(&keys, sigmas, String::AreEqual);
      double res = 0.0;
      auto tmpKeys = keys;
      while (tmpKeys != NULL) {
         String::String key = tmpKeys->GetValue();
         auto val = LinkedListKV::GetValue(sigmas, key, String::AreEqual);
         res += val * val;
         tmpKeys = tmpKeys->GetNext();
      }
      auto tmpKeys1 = keys, tmpKeys2 = tmpKeys1->GetNext();
      while (tmpKeys1 != NULL) {
         while (tmpKeys2 != NULL) {
            auto key1 = tmpKeys1->GetValue();
            auto key2 = tmpKeys2->GetValue();
            auto sigma1 = LinkedListKV::GetValue(sigmas, key1, String::AreEqual);
            auto sigma2 = LinkedListKV::GetValue(sigmas, key2, String::AreEqual);
            res += 2.0 * sigma1 * sigma2 * corr.get(key1, key2);
            tmpKeys2 = tmpKeys2->GetNext();
         }
         tmpKeys1 = tmpKeys1->GetNext();
         tmpKeys2 = tmpKeys1->GetNext();
      }
            
      return res;
   }

   __device__ double UncertainValue2::uncertainty(Correlations corr)
   {
      return ::sqrt(variance(corr));
   }

   __device__ void UncertainValue2::PrintSigmas()
   {
      auto sigmas = mSigmas;
      while (sigmas != NULL) {
         printf("%s: %lf\n", sigmas->GetKey().Get(), sigmas->GetValue());
         sigmas = sigmas->GetNext();
      }
   }
}
