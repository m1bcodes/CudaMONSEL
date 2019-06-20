#include "Amphibian\Tests\StackTest.cuh"

#include "Amphibian\stack.cuh"

namespace StackTest
{
   __host__ __device__ static void assertTrue(bool expr)
   {
      if (!expr) printf("wrong!\n");
   }

   __host__ __device__ void testOne()
   {
      amp::stack<int> s;

      s.push(0);
      s.push(1);
      s.push(2);

      while (!s.empty()) {
         printf("%d ", s.top());
         s.pop();
      }
      printf("\n");

      printf("StackTest::testOne() completed.\n");
   }
}