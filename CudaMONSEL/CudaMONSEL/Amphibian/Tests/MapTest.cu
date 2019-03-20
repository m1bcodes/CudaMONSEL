#include "MapTest.cuh"

#include "..\String.cuh"

#include <math.h>

namespace MapTest
{
   __device__ void AssertEqual(int a, int b)
   {
      if (a != b) {
         printf("not equal: (%d, %d)", a, b);
      }
   }

   __device__ MapTest::MapTest() : DefaultHasher(Hasher::APHash)
   {
   }

   __device__ void MapTest::TestInteger()
   {
      Map::Map<int, int> m1 = CreateMapA<int, int>([](int& a, int& b) { return a == b; }, [](int& a, int& b) { return a == b; }, 0, 1);
      m1.Put(1, 1);

      AssertEqual(m1.Size(), 2);

      Map::Map<int, int> m2(m1);
      AssertEqual(m2.Size(), 2);
      m2.Put(2, 3);
      AssertEqual(m2.Size(), 3);
      AssertEqual(m1.Size(), 2);

      int c1 = 0;
      Map::Iterator<int, int> itr1(m1);
      while (itr1.HasNext()) {
         //printf("(%d, %d) ", itr1.GetKey(), itr1.GetValue());
         c1++;
         itr1.Next();
      }
      AssertEqual(c1, 2);

      int c2 = 0;
      Map::Iterator<int, int> itr2(m2);
      while (itr2.HasNext()) {
         //printf("(%d, %d) ", itr2.GetKey(), itr2.GetValue());
         c2++;
         itr2.Next();
      }
      AssertEqual(c2, 3);

      auto m3 = m2;
      AssertEqual(m3.Size(), 3);
      m3.Put(9, 1);
      AssertEqual(m3.Size(), 4);
      auto h2 = m2.HashCode();
      auto h3 = m3.HashCode();
      if (h2 == h3) {
         printf("HashCodes are different: %d, %d\n", h2, h3);
      }

      printf("MapTest::TestInteger() completed\n");
   }

   __device__ Map::Map<String::String, double> makeMap()
   {
      Map::Map<String::String, double> map(Hasher::APHash, String::AreEqual, [](double& a, double& b) {return a == b; });
      map.Put("A", 1);
      map.Put("B", 2);
      map.Put("C", 3);
      map.Put("D", 4);
      map.Put("E", 5);
      return map;
   }

   __device__ void MapTest::TestString()
   {
      Map::Map<String::String, double> m1(DefaultHasher, String::AreEqual, [](double& a, double& b) { return a == b; });
      m1.Put("A", 1);
      m1.Put("B", 2);

      AssertEqual(m1.Size(), 2);

      Map::Map<String::String, double> m2(m1);
      AssertEqual(m2.Size(), 2);
      m2.Put("C", 3);
      AssertEqual(m2.Size(), 3);
      AssertEqual(m1.Size(), 2);

      int c1 = 0;
      Map::Iterator<String::String, double> itr1(m1);
      while (itr1.HasNext()) {
         //printf("(%s, %lf) ", itr1.GetKey().Get(), itr1.GetValue());
         c1++;
         itr1.Next();
      }
      AssertEqual(c1, 2);
      //printf("\n");

      int c2 = 0;
      Map::Iterator<String::String, double> itr2(m2);
      while (itr2.HasNext()) {
         //printf("(%s, %lf) ", itr2.GetKey().Get(), itr2.GetValue());
         c2++;
         itr2.Next();
      }
      AssertEqual(c2, 3);
      //printf("\n");

      Map::Map<String::String, double> m3 = m2;
      AssertEqual(m3.Size(), 3);
      m3.Put("D", 4);
      AssertEqual(m3.Size(), 4);

      int c3 = 0;
      Map::Iterator<String::String, double> itr3(m3);
      while (itr3.HasNext()) {
         //printf("(%s, %lf) ", itr3.GetKey().Get(), itr3.GetValue());
         c3++;
         itr3.Next();
      }
      AssertEqual(m3.Size(), 4);
      //printf("\n");

      Map::Map<String::String, double> map4(DefaultHasher, String::AreEqual, [](double& a, double& b) { return a == b; });
      map4.Put("V1", ::sqrt(0.05));
      map4.Put("V2", ::sqrt(0.04));

      int c4 = 0;
      Map::Iterator<String::String, double> itr4(map4);
      while (itr4.HasNext()) {
         //printf("(%s, %lf) ", itr4.GetKey().Get(), itr4.GetValue());
         c4++;
         itr4.Next();
      }
      AssertEqual(map4.Size(), 2);
      AssertEqual(c4, 2);

      auto m5 = map4;
      if (!(m5 == map4)) {
         printf("maps sizes are different: %d, %d\n", map4.Size(), m5.Size());
      }
      if (!(m5 == map4)) {
         printf("maps are different\n");
      }
      if (m5.HashCode() != map4.HashCode()) {
         printf("maps hashcodes are different\n");
      }
      
      Map::Iterator<String::String, double> itr5(m5);
      while (itr5.HasNext()) {
         //printf("(%s, %lf) ", itr5.GetKey().Get(), itr5.GetValue());
         itr5.Next();
      }

      auto m6 = makeMap();
      Map::Iterator<String::String, double> itr6(m6);
      while (itr6.HasNext()) {
         //printf("(%s, %lf) ", itr6.GetKey().Get(), itr6.GetValue());
         itr6.Next();
      }
      //printf("\n");

      printf("MapTest::TestString() completed\n");
   }
}