#include "Amphibian\Tests\VectorTest.cuh"

#include "Amphibian\vector.cuh"

namespace VectorTest
{
   struct IntCmp
   {
      __host__ __device__ inline bool operator() (const int a, const int b)
      {
         return a == b;
      }
   };

   __host__ __device__ void assertTrue(bool cond)
   {
      if (!cond) printf("bad!\n");
   }

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
      assertTrue(v.operator==(v2));
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

      amp::vector<amp::vector<int>> u(100);
      for (int i = 0; i < 100; ++i) {
         amp::vector<int> tmp(100);
         for (int j = 0; j < 100; ++j) {
            tmp.push_back(j);
         }
         u.push_back(tmp);
      }

      printf("VectorTest::testThree() completed.\n");
   }

   __host__ __device__ void testFour()
   {
      amp::vector<int> v0(5, 8);
      for (auto i : v0) {
         printf("%d, ", i);
      }
      printf("\n");

      amp::vector<amp::vector<int>> v(5, amp::vector<int>(5, 8));
      for (auto u : v) {
         for (auto i : u) {
            printf("%d, ", i);
         }
         printf("\n");
      }

      v.resize(7, amp::vector<int>(9, 1));
      for (auto u : v) {
         for (auto i : u) {
            printf("%d, ", i);
         }
         printf("\n");
      }

      const int tmp[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
      v[0].assign(tmp, tmp + 11);
      for (auto i : v[0]) {
         printf("%d, ", i);
      }
      printf("\n");

      printf("VectorTest::testFour() completed.\n");
   }
}