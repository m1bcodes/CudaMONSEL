#include "MapTest.cuh"

#include "..\String.cuh"

#include <math.h>

namespace MapTest
{
   __device__ void AssertEqual(int a, int b)
   {
      if (a != b) {
         printf("not equal: (%d, %d)\n", a, b);
      }
   }

   __device__ MapTest::MapTest() : DefaultHasher(Hasher::APHash)
   {
   }

   __device__ void MapTest::TestInteger()
   {
      Map::Map<int, int> m1 = CreateMapA<int, int>([](int& a, int& b) { return a == b; }, [](int& a, int& b) { return a == b; }, 0, 1);
      int k0 = 1, v0 = 1;
      m1.Put(k0, v0);

      AssertEqual(m1.Size(), 2);

      Map::Map<int, int> m2(m1);
      AssertEqual(m2.Size(), 2);
      int k1 = 2, v1 = 3;
      m2.Put(k1, v1);
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
      int k2 = 2, v2 = 3;
      m3.Put(k2, v2);
      AssertEqual(m3.Size(), 3);
      auto h2 = m2.HashCode();
      auto h3 = m3.HashCode();
      if (h2 != h3) {
         printf("HashCodes are different: %d, %d\n", h2, h3);
      }

      unsigned int tmp1 = 1;
      auto h1 = m1.Hash((char*)&tmp1, sizeof(unsigned int));
      unsigned int tmp2 = 1;
      h2 = DefaultHasher((char*)&tmp2, sizeof(unsigned int));
      if (h1 != h2) {
         printf("value: (%d, %u), (%d, %u)\n", tmp1, h1, tmp2, h2);
      }

      printf("MapTest::TestInteger() completed\n");
   }

   __device__ Map::Map<String::String, double> makeMap()
   {
      Map::Map<String::String, double> map(Hasher::APHash, String::AreEqual, [](double& a, double& b) {return a == b; });
      double v0 = 1, v1 = 2, v2 = 3, v3 = 4, v4 = 5;
      map.Put(String::String("A"), v0);
      map.Put(String::String("B"), v1);
      map.Put(String::String("C"), v2);
      map.Put(String::String("D"), v3);
      map.Put(String::String("E"), v4);
      return map;
   }

   __device__ void MapTest::TestString()
   {
      Map::Map<String::String, double> m1(DefaultHasher, String::AreEqual, [](double& a, double& b) { return a == b; });
      double v0 = 1, v1 = 2, v2 = 3, v3 = 4;
      m1.Put(String::String("A"), v0);
      m1.Put(String::String("B"), v1);

      AssertEqual(m1.Size(), 2);

      Map::Map<String::String, double> m2(m1);
      AssertEqual(m2.Size(), 2);
      m2.Put(String::String("C"), v2);
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
      m3.Put(String::String("D"), v3);
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
      double v5 = ::sqrt(0.05), v6 = ::sqrt(0.04);
      map4.Put(String::String("V1"), v5);
      map4.Put(String::String("V2"), v6);

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
         printf("maps are different\n");
      }
      if (!(m5.Size() == map4.Size())) {
         printf("maps sizes are different: %d, %d\n", map4.Size(), m5.Size());
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

      unsigned int tmp1 = 1;
      auto h1 = m1.Hash((char*)&tmp1, sizeof(unsigned int));
      unsigned int tmp2 = 1;
      auto h2 = DefaultHasher((char*)&tmp2, sizeof(unsigned int));
      if (h1 != h2) {
         printf("value: (%d, %u), (%d, %u)\n", tmp1, h1, tmp2, h2);
      }

      printf("MapTest::TestString() completed\n");
   }

   __device__ bool doubleCmp(double& a, double& b) { return a == b; }
   
   class TestClass
   {
   public:
      __device__ TestClass() : m(Hasher::APHash, String::AreEqual, doubleCmp) {}
      //__device__ TestClass(int) : m(Hasher::APHash, String::AreEqual, doubleCmp) {}

      __device__  static bool AreEqualTestClass(TestClass& a, TestClass& b)
      {
         if (&a == &b) return true;
         typedef bool(*pAreEqualStrings)(String::String&, String::String&);
         pAreEqualStrings fcnptr = String::AreEqual;
         //printf("%p, %p\n", a.m.GetHasher(), Hasher::APHash);
         //printf("%p, %p\n", a.m.GetKeyCmp(), fcnptr);
         //printf("%p, %p\n", a.m.GetValCmp(), doubleCmp);

         unsigned int tmp1 = 1;
         //auto ahasher = a.m.GetHasher();
         auto h1 = Hasher::APHash((char*)&tmp1, sizeof(unsigned int));
         unsigned int tmp2 = 1;
         auto h2 = Hasher::APHash((char*)&tmp2, sizeof(unsigned int));
         printf("value: (%d, %u), (%d, %u)\n", tmp1, h1, tmp2, h2);

         printf("a size: %d\n", a.m.Size());
         printf("b size: %d\n", b.m.Size());
         Map::Iterator<String::String, double> itra(a.m);
         printf("A:\n");
         while (itra.HasNext()) {
            auto k = itra.GetKey();
            auto v = itra.GetValue();
            printf("%s %lf\n", k.Get(), v);
            //auto bm = b.m;
            //Map::Map<String::String, double> bm(b.m);
            //if (!b.m.ContainsKey(k)) return false;

            //if (v != b.m.GetValue(k)) return false;
            itra.Next();
         }
         printf("B:\n");
         Map::Iterator<String::String, double> itrb(b.m);
         while (itrb.HasNext()) {
            auto k = itrb.GetKey();
            auto v = itrb.GetValue();
            printf("%s %lf\n", k.Get(), v);
            itrb.Next();
         }
         printf("===\n");
         return true;
      }

      Map::Map<String::String, double> m;
   //private:
   };

   __device__ void MapTest::TestMapOfMap()
   {
      Map::Map<String::String, TestClass> m1(DefaultHasher, String::AreEqual, TestClass::AreEqualTestClass);
      TestClass a;
      TestClass b;

      unsigned int tmp1 = 1;
      auto h1 = a.m.Hash((char*)&tmp1, sizeof(unsigned int));
      unsigned int tmp2 = 1;
      auto h2 = DefaultHasher((char*)&tmp2, sizeof(unsigned int));
      if (h1 != h2) {
         printf("value: (%d, %u), (%d, %u)\n", tmp1, h1, tmp2, h2);
      }

      double v1 = 1.2, v2 = 2.1, v3 = 3.3;
      a.m.Put(String::String("A1"), v1);
      a.m.Put(String::String("A2"), v2);
      a.m.Put(String::String("A3"), v3);
      b.m.Put(String::String("B1"), v1);
      m1.Put(String::String("A"), a);
      m1.Put(String::String("B"), b);
      TestClass::AreEqualTestClass(a, b);
      printf("MapTest::TestMapOfMap() completed\n");
   }
}