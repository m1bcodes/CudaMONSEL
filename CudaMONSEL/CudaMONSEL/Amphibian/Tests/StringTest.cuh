#ifndef _STRING_TEST_CUH_
#define _STRING_TEST_CUH_

#include "Amphibian\String.cuh"

namespace StringTest
{
   __host__ __device__ void assertEqual(bool expr);
   __host__ __device__ void EmptyTest();
   __host__ __device__ void TestOne();
   __host__ __device__ void AtoITest();
   __host__ __device__ void AtoFTest();
   __host__ __device__ void ItoATest();
   __host__ __device__ void findTest();
   __host__ __device__ void addTest();
};

#endif