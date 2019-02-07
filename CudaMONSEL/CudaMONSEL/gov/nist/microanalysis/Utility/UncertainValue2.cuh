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
      __device__ UncertainValue2(double v, LinkedListKV::Node<String, double>* sigmas);

      __device__ void assignComponent(String name, double sigma);
      __device__ double getComponent(String src);
      __device__ LinkedListKV::Node<String, double> const * getComponents() const;
      __device__ bool hasComponent(String src);
      __device__ void renameComponent(String oldName, String newName);

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

      //__device__ static UncertainValue2 sqr(UncertainValue2 uv);
      //__device__ static UncertainValue2 negate(UncertainValue2 uv);
      //__device__ static UncertainValue2 atan(UncertainValue2 uv);
      //__device__ static UncertainValue2 atan2(UncertainValue2 y, UncertainValue2 x);
      //__device__ static UncertainValue2 positiveDefinite(UncertainValue2 uv);

   private:
      double mValue;
      LinkedListKV::Node<String, double>* mSigmas;
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

   //__device__ add(Node<UncertainValue2> uvs;
}

#endif