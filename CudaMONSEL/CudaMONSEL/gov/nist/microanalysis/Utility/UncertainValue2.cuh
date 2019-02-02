#ifndef _UNCERTAIN_VALUE_2_CUH_
#define _UNCERTAIN_VALUE_2_CUH_

#include "..\..\..\..\Amphibian\String.cuh"

#include <thrust\device_vector.h>

class UncertainValue2
{
public:
   typedef struct Sigma {
      String name;
      double val;
   };

   static const char DEFAULT[];

   static const UncertainValue2 ONE;
   static const UncertainValue2 ZERO;
   //static const UncertainValue2 NaN;
   //static const UncertainValue2 POSITIVE_INFINITY;
   //static const UncertainValue2 NEGATIVE_INFINITY;

   __device__ UncertainValue2(double v, char source[], double dv);
   __device__ UncertainValue2(double v);
   __device__ UncertainValue2(double v, double dv);
   __device__ UncertainValue2(double v, Node<String, double>* sigmas);

   __device__ void assignComponent(String name, double sigma);
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

   __device__ static UncertainValue2 sqr(UncertainValue2 uv);
   __device__ static UncertainValue2 negate(UncertainValue2 uv);
   __device__ static UncertainValue2 atan(UncertainValue2 uv);
   __device__ static UncertainValue2 atan2(UncertainValue2 y, UncertainValue2 x);
   __device__ static UncertainValue2 positiveDefinite(UncertainValue2 uv);

private:
   static const int MAX_LEN = 11;
   static const long long serialVersionUID;

   const double mValue;
   thrust::device_vector<Sigma> mSigmas;

   //bool mNotANumber;
   //bool mPosInfinity, mNegInfinity;
};

#endif