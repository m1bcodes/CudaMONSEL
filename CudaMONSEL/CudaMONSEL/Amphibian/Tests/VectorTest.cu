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
      if (v.operator==(v2)) printf("good!\n");
      v.clear();
      for (auto i : v) {
         for (auto j : i) {
            printf("%d, ", j);
         }
         printf("\n");
      }
      for (auto i : v2) {
         for (auto j : i) {
            printf("%d, ", j);
         }
         printf("\n");
      }

      int a[] = { 0, 1, 2 };
      amp::vector<int> v3(a, a + 3);
      for (auto i : v3) {
         printf("%d, ", i);
      }
      printf("\n");

      auto itr = amp::find(v3.begin(), v3.end(), 1);
      v3.erase(itr);
      for (auto i : v3) {
         printf("%d, ", i);
      }
      printf("\n");

      printf("VectorTest::testTwo() completed.\n");
   }

   __host__ __device__ void testThree()
   {
      amp::vector<int> v;
      for (int i = 0; i < 100000; ++i) {
         v.push_back(i);
         //printf("%d ", v[0]);
         auto itr = amp::find(v.begin(), v.end(), i);
         //v.clear();
         v.erase(itr);
      }
      //printf("\n");
      printf("VectorTest::testThree() completed.\n");
   }
}