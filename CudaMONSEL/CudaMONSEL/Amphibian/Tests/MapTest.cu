#include "Amphibian/Tests/MapTest.cuh"

#include "Amphibian/String.cuh"

#include <math.h>

namespace MapTest
{
   __host__ __device__ static void AssertEqual(int a, int b)
   {
      if (a != b) {
         printf("not equal: (%d, %d)\n", a, b);
      }
   }

   __host__ __device__ MapTest::MapTest()
   {
   }

   typedef amp::unordered_map<int, int, Comparator::IntCompareFcn, Comparator::IntCompareFcn, Hasher::IntHashFcn, Hasher::IntHashFcn> IntTestType;
   typedef amp::unordered_map<int, int, Comparator::IntCompareFcn, Comparator::IntCompareFcn, Hasher::IntHashFcn, Hasher::IntHashFcn>::iterator IntTestTypeItr;

   __host__ __device__ void MapTest::testInteger()
   {
      int k = 0, v = 1;
      IntTestType m1 = CreateMapA<int, int, Comparator::IntCompareFcn, Comparator::IntCompareFcn, Hasher::IntHashFcn, Hasher::IntHashFcn>(k, v);
      int k0 = 1, v0 = 1;
      m1.Put(k0, v0);

      AssertEqual(m1.size(), 2);

      IntTestType m2(m1);
      AssertEqual(m2.size(), 2);
      int k1 = 2, v1 = 3;
      //m2.Put(k1, v1);
      m2.insert(amp::make_pair<int, int>(k1, v1));
      AssertEqual(m2.size(), 3);
      AssertEqual(m1.size(), 2);

      int c1 = 0;
      IntTestTypeItr itr1(m1);
      while (itr1.HasNext()) {
         printf("(%d, %d) ", itr1.GetKey(), itr1.GetValue());
         c1++;
         itr1.Next();
      }
      printf("\n");

      for (auto itr = m1.begin(); itr != m1.end(); ++itr) {
         printf("(%d, %d) ", itr.GetKey(), itr.GetValue());
      }
      printf("\n");
      AssertEqual(c1, 2);

      int c2 = 0;
      IntTestTypeItr itr2(m2);
      while (itr2.HasNext()) {
         printf("(%d, %d) ", itr2.GetKey(), itr2.GetValue());
         c2++;
         itr2.Next();
      }
      printf("\n");
      for (auto itr = m2.begin(); itr != m2.end(); ++itr) {
         printf("(%d, %d) ", itr.GetKey(), itr.GetValue());
      }
      printf("\n");

      AssertEqual(c2, 3);

      IntTestType m3 = m2;
      AssertEqual(m3.size(), 3);
      int k2 = 2, v2 = 3;
      m3.Put(k2, v2);
      AssertEqual(m3.size(), 3);
      unsigned int h2 = m2.hashCode();
      unsigned int h3 = m3.hashCode();
      if (h2 != h3) {
         printf("HashCodes are different: %d, %d\n", h2, h3);
      }

      printf("MapTest::testInteger() completed\n");
   }

   typedef amp::unordered_map<amp::string, double, amp::string_cmp, Comparator::DoubleCompareFcn, amp::string_hash, Hasher::DoubleHashFcn> StringTestType;
   typedef amp::unordered_map<amp::string, double, amp::string_cmp, Comparator::DoubleCompareFcn, amp::string_hash, Hasher::DoubleHashFcn>::iterator StringTestTypeItr;

   __host__ __device__ StringTestType makeMap()
   {
      StringTestType map;
      double v0 = 1, v1 = 2, v2 = 3, v3 = 4, v4 = 5;
      amp::string A("A"), B("B"), C("C"), D("D"), E("E");
      map.Put(A, v0);
      map.Put(B, v1);
      map.Put(C, v2);
      map.Put(D, v3);
      map.Put(E, v4);
      return map;
   }

   __host__ __device__ void MapTest::testString()
   {
      StringTestType m1;
      double v0 = 1, v1 = 2, v2 = 3, v3 = 4;
      amp::string A("A"), B("B"), C("C"), D("D"), E("E");
      m1.Put(A, v0);
      m1.Put(B, v1);

      AssertEqual(m1.size(), 2);

      StringTestType m2(m1);
      AssertEqual(m2.size(), 2);
      m2.Put(C, v2);
      AssertEqual(m2.size(), 3);
      AssertEqual(m1.size(), 2);

      int c1 = 0;
      StringTestTypeItr itr1(m1);
      while (itr1.HasNext()) {
         //printf("(%s, %lf) ", itr1.GetKey().Get(), itr1.GetValue());
         c1++;
         itr1.Next();
      }
      AssertEqual(c1, 2);
      //printf("\n");

      int c2 = 0;
      StringTestTypeItr itr2(m2);
      while (itr2.HasNext()) {
         //printf("(%s, %lf) ", itr2.GetKey().Get(), itr2.GetValue());
         c2++;
         itr2.Next();
      }
      AssertEqual(c2, 3);
      //printf("\n");

      StringTestType m3 = m2;
      AssertEqual(m3.size(), 3);
      m3.Put(D, v3);
      AssertEqual(m3.size(), 4);

      int c3 = 0;
      StringTestTypeItr itr3(m3);
      while (itr3.HasNext()) {
         //printf("(%s, %lf) ", itr3.GetKey().Get(), itr3.GetValue());
         c3++;
         itr3.Next();
      }
      AssertEqual(m3.size(), 4);
      //printf("\n");

      StringTestType map4;
      double v5 = ::sqrt(0.05), v6 = ::sqrt(0.04);
      amp::string V1("V1"), V2("V2");
      map4.Put(V1, v5);
      map4.Put(V2, v6);

      int c4 = 0;
      StringTestTypeItr itr4(map4);
      while (itr4.HasNext()) {
         //printf("(%s, %lf) ", itr4.GetKey().Get(), itr4.GetValue());
         c4++;
         itr4.Next();
      }
      AssertEqual(map4.size(), 2);
      AssertEqual(c4, 2);

      StringTestType m5 = map4;
      if (!(m5 == map4)) {
         printf("maps are different\n");
      }
      if (!(m5.size() == map4.size())) {
         printf("maps sizes are different: %d, %d\n", map4.size(), m5.size());
      }
      if (m5.hashCode() != map4.hashCode()) {
         printf("maps hashcodes are different\n");
      }

      StringTestTypeItr itr5(m5);
      while (itr5.HasNext()) {
         //printf("(%s, %lf) ", itr5.GetKey().Get(), itr5.GetValue());
         itr5.Next();
      }

      StringTestType m6 = makeMap();
      StringTestTypeItr itr6(m6);
      while (itr6.HasNext()) {
         //printf("(%s, %lf) ", itr6.GetKey().Get(), itr6.GetValue());
         printf("(%s, %lf) ", ((amp::string&)(itr6->first)).c_str(), (double)(itr6->second));
         itr6.Next();
      }
      printf("\n");
      for (auto p : m6) {
         printf("(%s, %lf) ", ((amp::string&)(p.first)).c_str(), (double)(p.second));
      }
      printf("\n");

      printf("MapTest::testString() completed\n");
   }

   class TestClass
   {
   public:
      __host__ __device__ TestClass() {}
      __host__ __device__ StringTestType& GetMap() { return m; }

   private:
      StringTestType m;
   };

   struct TestClassCompare
   {
      __host__ __device__ inline bool operator() (TestClass& lhs, TestClass& rhs)
      {
         return lhs.GetMap() == rhs.GetMap();
      }
   };

   struct TestClassHashFcn
   {
      __host__ __device__ inline unsigned int operator() (TestClass& t)
      {
         return t.GetMap().hashCode();
      }
   };

   typedef amp::unordered_map<amp::string, TestClass, amp::string_cmp, TestClassCompare, amp::string_hash, TestClassHashFcn> MapTestT;
   typedef amp::unordered_map<amp::string, TestClass, amp::string_cmp, TestClassCompare, amp::string_hash, TestClassHashFcn>::iterator MapTestTItr;

   __host__ __device__ void MapTest::testMapOfMap()
   {
      MapTestT m1;
      TestClass a;
      TestClass b;

      double v1 = 1.2, v2 = 1.2, v3 = 3.3;
      amp::string s0("A1"), s1("A2"), s2("A3"), s3("B1");
      amp::string s4("A"), s5("B");
      a.GetMap().Put(s0, v1);
      a.GetMap().Put(s1, v2);
      a.GetMap().Put(s2, v3);
      b.GetMap().Put(s3, v1);
      m1.Put(s4, a);
      m1.Put(s5, b);

      printf("MapTest::testMapOfMap() completed\n");
   }

   __host__ __device__ void MapTest::testAggregate()
   {
      StringTestType t;
      double v0 = 0.1, v1 = 0.5, v2 = 100;
      amp::string s0("A"), s1("B"), s2("C");
      t.Put(s0, v0);
      t.Put(s1, v1);
      t.Put(s2, v2);

      double ret1 = t.Aggregate([](double a) {return a*a; });
      double ret2 = t.Aggregate([](double a) {return a*a; });
      double ret3 = t.Aggregate([](double a) {return a*a; });
      double ans = 0.1*0.1 + 0.5*0.5 + 100 * 100;
      if (ret1 != ans) {
         printf("wrong : %lf, %lf", ret1, ans);
      }
      if (ret2 != ans) {
         printf("wrong : %lf, %lf", ret2, ans);
      }
      if (ret3 != ans) {
         printf("wrong : %lf, %lf", ret3, ans);
      }
      printf("MapTest::testAggregate() completed\n");
   }
}
