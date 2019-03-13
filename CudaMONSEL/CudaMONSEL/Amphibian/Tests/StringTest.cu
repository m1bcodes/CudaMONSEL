#include "StringTest.cuh"

#include <stdio.h>

#include <cuda_runtime.h>

namespace StringTest
{
   __host__ __device__ void assertTrue(bool expr)
   {
      if (!expr) {
         printf("not equal\n");
      }
   }

   __host__ __device__ void EmptyTest()
   {
      String::String s;
      char* a = s.Get();
      for (int k = 0; k < 32; ++k) {
         //int l = 0;
         //l |= a[k];
         //printf("%d ", a[k]);
         assertTrue(a[k] == '\0');
      }
      
      printf("StringTest::EmptyTest() completed\n");
   }

   __host__ __device__ void AtoITest()
   {
      assertTrue(String::AToI("0") == 0);
      assertTrue(String::AToI("1") == 1);
      assertTrue(String::AToI("12") == 12);
      assertTrue(String::AToI("123") == 123);
      assertTrue(String::AToI("-123") == -123);
      assertTrue(String::AToI("-2147483647") == -2147483647);
      printf("StringTest::AtoITest() completed\n");
   }

   __host__ __device__ void ItoATest()
   {
      char num[16] = { NULL };
      String::IToA(num, 0, 16);
      assertTrue(String::AreEqual(num, "0"));
      String::IToA(num, 1, 16);
      assertTrue(String::AreEqual(num, "1"));
      String::IToA(num, 12, 16);
      assertTrue(String::AreEqual(num, "12"));
      String::IToA(num, 123, 16);
      assertTrue(String::AreEqual(num, "123"));
      String::IToA(num, -1, 16);
      assertTrue(String::AreEqual(num, "-1"));
      String::IToA(num, -12, 16);
      assertTrue(String::AreEqual(num, "-12"));
      String::IToA(num, -123, 16);
      assertTrue(String::AreEqual(num, "-123"));
      String::IToA(num, -2147483647, 16);
      assertTrue(String::AreEqual(num, "-2147483647"));
      printf("StringTest::ItoATest() completed\n");
   }

   __host__ __device__ void AtoFTest()
   {
      assertTrue(String::AToF<float>("0") == 0.f);
      assertTrue(String::AToF<float>("1") == 1.f);
      assertTrue(String::AToF<float>("12") == 12.f);
      assertTrue(String::AToF<float>("123") == 123.f);
      assertTrue(String::AToF<float>("-123") == -123.f);
      assertTrue(String::AToF<float>("-2147483") == -2147483.f);
      assertTrue(String::AToF<float>("1.0") == 1.f);
      assertTrue(String::AToF<float>("1.1") == 1.1f);
      assertTrue(String::AToF<float>("10.2") == 10.2f);
      assertTrue(String::AToF<float>("10.21") == 10.21f);
      assertTrue(String::AToF<float>("101.212") == 101.212f);
      assertTrue(String::AToF<float>("101.2127") == 101.2127f);
      assertTrue(String::AToF<double>("2147483647") == 2147483647);
      assertTrue(String::AToF<double>("-2147483647") == -2147483647);

      printf("StringTest::AtoFTest() completed\n");
   }
}