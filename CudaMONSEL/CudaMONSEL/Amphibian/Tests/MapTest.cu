#include "MapTest.cuh"

#include "..\Map.cuh"
#include "..\Hasher.cuh"

namespace MapTest
{
   __device__ void TestA()
   {
      Hasher::pHasher hasher = Hasher::APHash;
      Map::Map<int, int> map(hasher, [](int a, int b) { return a == b; });
   }
}