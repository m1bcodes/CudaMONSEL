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

      amp::vector<int, IntCmp> v;
      v.push_back(0);
      v.push_back(1);
      v[0] = 5;

      for (auto i : v) {
         printf("%d\n", i);
      }

      printf("size: %d\n", v.size());
      v.clear();
      printf("cleared size: %d\n", v.size());

      printf("VectorTest::InsertTest() completed.\n");
   }
}