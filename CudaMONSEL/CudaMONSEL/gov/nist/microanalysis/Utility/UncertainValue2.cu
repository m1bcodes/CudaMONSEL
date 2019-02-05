#include "UncertainValue2.cuh"

#include <math_constants.h>

const char UncertainValue2::DEFAULT[] = "Default";

const long long UncertainValue2::serialVersionUID = 119495064970078787L;

__device__ int sDefIndex = 0; // transient

const UncertainValue2 UncertainValue2::ONE(1.0);
const UncertainValue2 UncertainValue2::ZERO(0.0);
//const UncertainValue2 UncertainValue2::NaN(CUDART_NAN);
//const UncertainValue2 UncertainValue2::POSITIVE_INFINITY(CUDART_INF);
//const UncertainValue2 UncertainValue2::NEGATIVE_INFINITY(-CUDART_INF);

__device__ UncertainValue2::UncertainValue2(double v, double dv) : mValue(v)
{
   char tmpName[MAX_LEN];
   String::IToA(tmpName, atomicAdd(&sDefIndex, 1));
   UncertainValue2::UncertainValue2(v, tmpName, dv);
}

__device__ UncertainValue2::UncertainValue2(double v) : mValue(v)
{
   UncertainValue2::UncertainValue2(v, 0.0);
}

__device__ UncertainValue2::UncertainValue2(double v, char source[], double dv) : mValue(v)
{
   assignComponent(source, dv);
}

__device__ UncertainValue2::UncertainValue2(double v, Node<String, double>* sigmas) : mValue(v)
{
   while (sigmas != NULL) {
      assignComponent(sigmas->GetKey(), sigmas->GetValue());
      sigmas = sigmas->GetNext();
   }
}

__device__ void UncertainValue2::assignComponent(String name, double sigma)
{
   if (sigma != 0.0) {
      Node<String, double>::InsertHead(&mSigmas, name, sigma);
   }
   else {
      Node<String, double>::Remove(&mSigmas, name, String::AreEqual);
   }
}

__device__ double UncertainValue2::getComponent(String src)
{
   double v = Node<String, double>::GetNode(mSigmas, src, String::AreEqual);
   return v != NULL ? v : 0.0;
}

__device__ bool UncertainValue2::hasComponent(String src)
{
   return Node<String, double>::GetNode(mSigmas, src, String::AreEqual) != NULL;
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
   return sqrt(variance());
}

__device__ double UncertainValue2::variance()
{
   double sigma2 = 0.0;
   Node<String, double>* head = mSigmas;
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

//int hashCode()
//{
//   return Objects.hash(mValue, mSigmas);
//}

__device__ bool UncertainValue2::equals(UncertainValue2 const * obj)
{
   if (this == obj) {
      return true;
   }
   if (obj == NULL) {
      return false;
   }
   UncertainValue2 other = (UncertainValue2)*obj;
   return Node<String, double>::AreEquivalentSets(mSigmas, other.mSigmas, String::AreEqual, [](double a, double b) { return a == b; }) && (mValue == other.mValue);
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

__device__ UncertainValue2 UncertainValue2::sqr(UncertainValue2 uv)
{
   return pow(uv, 2.0);
}

__device__ UncertainValue2 UncertainValue2::negate(UncertainValue2 uv)
{
   UncertainValue2 ret = UncertainValue2(-uv.mValue, uv.mSigmas);
   return ret;
}

__device__ UncertainValue2 UncertainValue2::atan(UncertainValue2 uv)
{
   double f = atan(uv.doubleValue());
   double df = 1.0 / (1.0 + uv.doubleValue() * uv.doubleValue());

   if (!(isnan(f) || isnan(df))) {
      UncertainValue2 res(f);
      Node<String, double>* sigmas = uv.mSigmas;
      while (sigmas != NULL) {
         res.assignComponent(sigmas->GetKey(), df * sigmas->GetValue());
         sigmas = sigmas->GetNext();
      }
      return res;
   }
   else {
      UncertainValue2 nan(CUDART_NAN);
      return nan;
   }
}

__device__ UncertainValue2 UncertainValue2::atan2(UncertainValue2 y, UncertainValue2 x)
{
   double f = atan2(y.doubleValue(), x.doubleValue());
   double df = 1.0 / (1.0 + sqr(y.doubleValue() / x.doubleValue()));

   if (!(isnan(f) || isnan(df))) {
      UncertainValue2 res(f);
      Node<String, double>* sigmas;
      while (sigmas != NULL) {
         res.assignComponent(sigmas->getKey(), df * sigmas->getValue());
         sigmas = sigmas->GetNext();
      }
      return res;
   }
   else {
      UncertainValue2 nan(CUDART_NAN);
      return nan;
   }
}

__device__ UncertainValue2 UncertainValue2::positiveDefinite(UncertainValue2 uv)
{
   UncertainValue2 ret(0.0, uv.mSigmas);
   return uv.doubleValue() >= 0.0 ? uv : ret;
}
