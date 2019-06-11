#include "Amphibian/Tests/StringTest.cuh"

#include <stdio.h>
#include <string.h>

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
      amp::string s;
      const char* a = s.c_str();
      for (int k = 0; k < 32; ++k) {
         //int l = 0;
         //l |= a[k];
         //printf("%d ", a[k]);
         assertTrue(a[k] == '\0');
      }

      printf("StringTest::EmptyTest() completed\n");
   }

   __host__ __device__ void TestOne()
   {
      amp::string a("a");
      amp::string a1 = a;
      amp::string a2(a);
      assertTrue(a == a1);
      assertTrue(a == a2);
      assertTrue(a1 == a2);

      printf("StringTest::TestOne() completed\n");
   }

   __host__ __device__ void AtoITest()
   {
      assertTrue(amp::AToI("0") == 0);
      assertTrue(amp::AToI("1") == 1);
      assertTrue(amp::AToI("12") == 12);
      assertTrue(amp::AToI("123") == 123);
      assertTrue(amp::AToI("-123") == -123);
      assertTrue(amp::AToI("-2147483647") == -2147483647);

      printf("StringTest::AtoITest() completed\n");
   }

   __host__ __device__ void ItoATest()
   {
      const int MAX_LEN = 16;
      char num[MAX_LEN];
      memset(num, NULL_CHAR, MAX_LEN);

      amp::IToA(0, num, MAX_LEN);
      assertTrue(amp::equal(num, "0"));
      amp::IToA(1, num, MAX_LEN);
      assertTrue(amp::equal(num, "1"));
      amp::IToA(12, num, MAX_LEN);
      assertTrue(amp::equal(num, "12"));
      amp::IToA(123, num, MAX_LEN);
      assertTrue(amp::equal(num, "123"));
      amp::IToA(-1, num, MAX_LEN);
      assertTrue(amp::equal(num, "-1"));
      amp::IToA(-12, num, MAX_LEN);
      assertTrue(amp::equal(num, "-12"));
      amp::IToA(-123, num, MAX_LEN);
      assertTrue(amp::equal(num, "-123"));
      amp::IToA(-2147483647, num, MAX_LEN);
      assertTrue(amp::equal(num, "-2147483647"));

      printf("StringTest::ItoATest() completed\n");
   }

   __host__ __device__ void AtoFTest()
   {
      assertTrue(amp::AToF<float>("0") == 0.f);
      assertTrue(amp::AToF<float>("1") == 1.f);
      assertTrue(amp::AToF<float>("12") == 12.f);
      assertTrue(amp::AToF<float>("123") == 123.f);
      assertTrue(amp::AToF<float>("-123") == -123.f);
      assertTrue(amp::AToF<float>("-2147483") == -2147483.f);
      assertTrue(amp::AToF<float>("1.0") == 1.f);
      assertTrue(amp::AToF<float>("1.1") == 1.1f);
      assertTrue(amp::AToF<float>("10.2") == 10.2f);
      assertTrue(amp::AToF<float>("10.21") == 10.21f);
      assertTrue(amp::AToF<float>("101.212") == 101.212f);
      assertTrue(amp::AToF<float>("101.2127") == 101.2127f);
      assertTrue(amp::AToF<double>("2147483647") == 2147483647);
      assertTrue(amp::AToF<double>("-2147483647") == -2147483647);

      printf("StringTest::AtoFTest() completed\n");
   }
}
