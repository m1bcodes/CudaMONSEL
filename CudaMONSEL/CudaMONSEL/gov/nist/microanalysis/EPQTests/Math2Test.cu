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
}