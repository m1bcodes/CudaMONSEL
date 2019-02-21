#include "StringTest.cuh"

#include <stdio.h>

#include <cuda_runtime.h>

namespace StringTest
{
   __host__ __device__ void assertEqual(bool expr)
   {
      if (!expr) {
         printf("not equal");
      }
   }

   __host__ __device__ void AtoITest()
   {
      assertEqual(String::AToI("0") == 0);
      assertEqual(String::AToI("1") == 1);
      assertEqual(String::AToI("12") == 12);
      assertEqual(String::AToI("123") == 123);
      assertEqual(String::AToI("-123") == -123);
      assertEqual(String::AToI("-2147483647") == -2147483647);
      printf("StringTest::AtoITest() completed");
   }

   __host__ __device__ void ItoATest()
   {
      char num[16];
      String::IToA(num, 0, 16);
      assertEqual(num == "0");
      String::IToA(num, 1, 16);
      assertEqual(num == "1");
      String::IToA(num, 12, 16);
      assertEqual(num == "12");
      String::IToA(num, 123, 16);
      assertEqual(num == "123");
      String::IToA(num, -123, 16);
      assertEqual(num == "-123");
      String::IToA(num, -2147483647, 16);
      assertEqual(num == "-2147483647");
      printf("StringTest::ItoATest() completed");
   }
}