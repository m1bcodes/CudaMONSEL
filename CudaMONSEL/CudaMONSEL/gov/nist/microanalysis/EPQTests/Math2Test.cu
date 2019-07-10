#include "gov\nist\microanalysis\EPQTests\Math2Test.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

#include "Amphibian\random.cuh"

#include <stdio.h>"

namespace Math2Test
{
   void testRandom1()
   {
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());

      printf("Math2Test::testRandom1() completed.\n");
   }

   void testRandom2()
   {
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());

      printf("Math2Test::testRandom2() completed.\n");
   }

   __device__ void testRandom1Cuda()
   {
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());

      printf("Math2Test::testRandom1CUDA() completed.\n");
   }

   __device__ void testRandom2Cuda()
   {
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());
      printf("%.10e\n", Random::random());

      printf("Math2Test::testRandom2CUDA() completed.\n");
   }
}