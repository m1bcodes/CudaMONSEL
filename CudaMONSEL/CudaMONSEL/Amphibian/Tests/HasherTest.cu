#include "Amphibian/Tests/HasherTest.cuh"

#include "Amphibian/Hasher.cuh"

#include <stdio.h>

namespace HasherTest
{
   __host__ __device__ void Equals(unsigned int hash, unsigned int target)
   {
      if (hash != target) {
         printf("Hash not equal: (%u, %u)\n", hash, target);
      }
   }

   __host__ __device__ void TestOne()
   {
      char* key = "abcdefghijklmnopqrstuvwxyz1234567890";

      Equals(Hasher::RSHash(key, 36), 4097835502);
      Equals(Hasher::JSHash(key, 36), 1651003062);
      Equals(Hasher::PJWHash(key, 36), 126631744);
      Equals(Hasher::ELFHash(key, 36), 126631744);
      Equals(Hasher::BKDRHash(key, 36), 3153586616);
      Equals(Hasher::SDBMHash(key, 36), 3449571336);
      Equals(Hasher::DJBHash(key, 36), 729241521);
      Equals(Hasher::DEKHash(key, 36), 2923964919);
      Equals(Hasher::BPHash(key, 36), 1726880944);
      Equals(Hasher::FNVHash(key, 36), 3243095106);
      Equals(Hasher::APHash(key, 36), 882643939);

      printf("HasherTest::TestOne() completed.\n");
   }
}
