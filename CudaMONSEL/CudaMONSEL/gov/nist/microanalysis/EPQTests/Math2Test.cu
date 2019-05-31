#include "gov\nist\microanalysis\EPQTests\Math2Test.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

#include <stdio.h>"

namespace Math2Test
{
   void testRandom1()
   {
      printf("%.10e\n", Math2::random());
      printf("%.10e\n", Math2::random());
      printf("%.10e\n", Math2::random());
      printf("%.10e\n", Math2::random());
      printf("%.10e\n", Math2::random());

      printf("Math2Test::testRandom1() completed.\n");
   }

   void testRandom2()
   {
      printf("%.10e\n", Math2::random());
      printf("%.10e\n", Math2::random());
      printf("%.10e\n", Math2::random());
      printf("%.10e\n", Math2::random());
      printf("%.10e\n", Math2::random());

      printf("Math2Test::testRandom2() completed.\n");
   }

   __device__ void testRandom1CUDA(curandState& state)
   {
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));

      printf("Math2Test::testRandom1CUDA() completed.\n");
   }

   __device__ void testRandom2CUDA(curandState& state)
   {
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));
      printf("%.10e\n", Math2::random(state));

      printf("Math2Test::testRandom2CUDA() completed.\n");
   }
}