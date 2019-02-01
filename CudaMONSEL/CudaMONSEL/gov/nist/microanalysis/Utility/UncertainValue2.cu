#include "UncertainValue2.cuh"
#include "..\..\..\..\Amphibian\Math.cuh"

const char UncertainValue2::DEFAULT[] = "Default";

//int UncertainValue2::sDefIndex = 0;
////const long UncertainValue2::serialVersionUID = 119495064970078787L;

//const UncertainValue2 UncertainValue2::ONE(1.0);
//const UncertainValue2 UncertainValue2::ZERO(0.0);
//const UncertainValue2 UncertainValue2::NaN(Double.NaN);
//const UncertainValue2 UncertainValue2::POSITIVE_INFINITY(Double.POSITIVE_INFINITY);
//const UncertainValue2 UncertainValue2::NEGATIVE_INFINITY(Double.NEGATIVE_INFINITY);

//__host__ __device__ UncertainValue2::UncertainValue2(double v, double dv) : mValue(v)
//{
//   char tmpName[MAX_LEN];
//   String::IToA(tmpName, ++sDefIndex);
//   UncertainValue2::UncertainValue2(v, tmpName, dv);
//}

//__host__ __device__ UncertainValue2::UncertainValue2(double v) : mValue(v)
//{
//   UncertainValue2::UncertainValue2(v, 0.0);
//}

__host__ __device__ UncertainValue2::UncertainValue2(double v, char source[], double dv) : mValue(v)
{
   assignComponent(source, dv);
}

__host__ __device__ UncertainValue2::UncertainValue2(double v, Node<String, double>* sigmas) : mValue(v)
{
   while (sigmas != NULL) {
      assignComponent(sigmas->GetKey(), sigmas->GetValue());
      sigmas = sigmas->GetNext();
   }
}

__host__ __device__ void UncertainValue2::assignComponent(String name, double sigma)
{
   if (sigma != 0.0) {
      Node<String, double>::InsertHead(&mSigmas, name, sigma);
   }
   else {
      Node<String, double>::Remove(&mSigmas, name, String::AreEqual);
   }
}
//
//double UncertainValue2::doubleValue() {
//   return mValue;
//}
//
//bool UncertainValue2::isUncertain() {
//   return mSigmas != NULL;
//}
//
//double UncertainValue2::uncertainty() {
//   double varnc = variance();
//   return Math::sqrt(varnc, varnc / 1000000);
//}
//
//double UncertainValue2::variance() {
//   double sigma2 = 0.0;
//   Node<String, double>* head = mSigmas;
//   while(head != NULL) {
//      double v = head->GetValue();
//      sigma2 += v * v;
//      head = head->GetNext();
//   }
//   return sigma2;
//}
//
//double UncertainValue2::fractionalUncertainty() {
//   return (mValue == 0 || isNaN()) ? 0 : Math::abs(uncertainty() / mValue);
//}
//
//bool UncertainValue2::isNaN() {
//   return mNotANumber;
//}
//
////int hashCode() {
////   return Objects.hash(mValue, mSigmas);
////}
//
//bool UncertainValue2::equals(UncertainValue2 * const obj) {
//   if (this == obj) {
//      return true;
//   }
//   if (obj == NULL) {
//      return false;
//   }
//   UncertainValue2 other = (UncertainValue2)*obj;
//   return Node<String, double>::AreEquivalentSets(mSigmas, other.mSigmas, String::AreEqual, [](double a, double b) { return a == b; }) && (mValue == other.mValue);
//}
