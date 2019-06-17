#include "Amphibian\Tests\VectorTest.cuh"

#include "Amphibian\vector.cuh"

namespace VectorTest
{
   struct IntCmp
   {
      __host__ __device__ inline bool operator() (const int a, const int b) {
         return a == b;
      }
   };

   __host__ __device__ void testOne()
   {
      int *ptr = new int(9);
      printf("%d\n", *ptr);

      amp::vector<int> v;
      v.push_back(0);
      v.push_back(1);
      v[0] = 5;

      for (auto i : v) {
         printf("%d\n", i);
      }

      printf("size: %d\n", v.size());
      v.clear();
      printf("cleared size: %d\n", v.size());

      printf("VectorTest::testOne() completed.\n");
   }

   __host__ __device__ void testTwo()
   {
      amp::vector<amp::vector<int>> v;
      amp::vector<int> u0;
      amp::vector<int> u1;
      u0.push_back(0);
      u0.push_back(1);
      u0.push_back(2);

      u1.push_back(10);
      u1.push_back(11);
      u1.push_back(12);

      v.push_back(u0);
      v.push_back(u1);

      for (auto i : v) {
         for (auto j : i) {
            printf("%d, ", j);
         }
         printf("\n");
      }

      amp::vector<amp::vector<int>> v2(v);
      v.clear();
      for (auto i : v2) {
         for (auto j : i) {
            printf("%d, ", j);
         }
         printf("\n");
      }

      printf("VectorTest::testTwo() completed.\n");
   }
}