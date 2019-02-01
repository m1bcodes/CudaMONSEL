#include "UncertainValue2.cuh"
#include <math_constants.h>

const char UncertainValue2::DEFAULT[] = "Default";

__device__ int sDefIndex = 0; // transient
const long long UncertainValue2::serialVersionUID = 119495064970078787L;

const UncertainValue2 UncertainValue2::ONE(1.0);
const UncertainValue2 UncertainValue2::ZERO(0.0);
//const UncertainValue2 UncertainValue2::NaN(Double.NaN);
//const UncertainValue2 UncertainValue2::POSITIVE_INFINITY(Double.POSITIVE_INFINITY);
//const UncertainValue2 UncertainValue2::NEGATIVE_INFINITY(Double.NEGATIVE_INFINITY);

UncertainValue2::UncertainValue2(double v, double dv) : mValue(v)
{
   char tmpName[MAX_LEN];
   String::IToA(tmpName, sDefIndex);
   atomicAdd(sDefIndex, 1);
   UncertainValue2::UncertainValue2(v, tmpName, dv);
}

UncertainValue2::UncertainValue2(double v) : mValue(v)
{
   UncertainValue2::UncertainValue2(v, 0.0);
}

UncertainValue2::UncertainValue2(double v, String source, double dv) : mValue(v)
{
   assignComponent(source, dv);
}

UncertainValue2::UncertainValue2(double v, thrust::device_vector<Sigma> sigmas) : mValue(v)
{
   for (thrust::device_vector<struct Sigma>::iterator iter = sigmas.begin(); iter < sigmas.end(); ++iter) {
      struct Sigma curr = *iter;
      assignComponent(curr.name, curr.val);
   }
}

void UncertainValue2::assignComponent(String name, double sigma)
{
   if (sigma != 0.0) {
      Sigma s = { name, sigma };
      mSigmas.push_back(s);
   }
   else {
      for (thrust::device_vector<struct Sigma>::iterator iter = mSigmas.begin(); iter < mSigmas.end(); ++iter) {
         struct Sigma curr = *iter;
         if (String::AreEqual(name, curr.name) && sigma == curr.val) {
            mSigmas.erase(iter);
            break;
         }
      }
   }
}

double UncertainValue2::doubleValue()
{
   return mValue;
}

bool UncertainValue2::isUncertain()
{
   return !mSigmas.empty();
}

double UncertainValue2::uncertainty() {
   double varnc = variance();
   return sqrt(varnc);
}

double UncertainValue2::variance() {
   double sigmaSOS = 0.0;
   for (thrust::device_vector<struct Sigma>::iterator iter = mSigmas.begin(); iter < mSigmas.end(); ++iter) {
      struct Sigma s = *iter;
      sigmaSOS += s.val * s.val;
   }
   return sigmaSOS;
}

double UncertainValue2::fractionalUncertainty() {
   return (isnan(1.0 / mValue)) ? CUDART_NAN : fabs(uncertainty() / mValue);
}

//int hashCode() {
//   return Objects.hash(mValue, mSigmas);
//}

bool UncertainValue2::equals(UncertainValue2 * const obj) {
   if (this == obj) {
      return true;
   }
   if (obj == NULL) {
      return false;
   }
   UncertainValue2 other = (UncertainValue2)*obj;
   return Node<String, double>::AreEquivalentSets(mSigmas, other.mSigmas, String::AreEqual, [](double a, double b) { return a == b; }) && (mValue == other.mValue);
}
