#include "UncertainValue2Test.cuh"

#include <stdio.h>
#include <math.h>

extern __device__ double __longlong_as_double(long long int);

namespace UncertainValue2Test
{
   UncertainValue2::UncertainValue2 makeA2a()
   {
      UncertainValue2::UncertainValue2 res(1.24);
      res.assignComponent("V1", ::sqrt(0.05));
      res.assignComponent("V2", ::sqrt(0.04));
      //printf("=============\n");
      //res.uncertainty();
      //printf("=============\n");
      return res;
   }

   UncertainValue2::UncertainValue2 makeB2a()
   {
      UncertainValue2::UncertainValue2 res(8.82);
      res.assignComponent("V1", ::sqrt(0.84));
      res.assignComponent("V2", ::sqrt(0.6));
      return res;
   }

   UncertainValue2::UncertainValue2 makeC2a()
   {
      UncertainValue2::UncertainValue2 res(-9.3);
      // 2.1 * 2.1 = 4.2 + 0.21 = 4.41
      res.assignComponent("V1", ::sqrt(3.0));
      res.assignComponent("V3", ::sqrt(1.41));
      return res;
   }

   //UncertainValue2::UncertainValue2 mA(1.24, "A", 0.3);
   //UncertainValue2::UncertainValue2 mA2a = makeA2a();

   //UncertainValue2::UncertainValue2 mB(8.82, "B", 1.2);
   //UncertainValue2::UncertainValue2 mB2a = makeB2a();

   //UncertainValue2::UncertainValue2 mC(-9.3, "C", 2.1);
   //UncertainValue2::UncertainValue2 mC2a = makeC2a();

   void assertEquals(double v1, double v2, double delta)
   {
      bool b = false;
      b = ::fabs(v1 - v2) < delta;
      if (!b) {
         b = v1 == v2;
      }
      if (!b) {
         b = true;
         for (int k = 0; k < sizeof(double); ++k) {
            char * c1 = ((char *)&v1) + k;
            char * c2 = ((char *)&v2) + k;
            b = !(*c1 ^ *c2);
            if (!b) {
               break;
            }
         }
      }
      if (!b) {
         printf("values are different: %lf, %lf\n", v1, v2);
      }
   }

   void assertEquals(UncertainValue2::UncertainValue2 uv, UncertainValue2::UncertainValue2 uv2, double delta)
   {
      bool b = false;
      b = ::fabs(uv.doubleValue() - uv2.doubleValue()) < delta;
      if (!b) {
         printf("values are different: %lf, %lf\n", uv.doubleValue(), uv2.doubleValue());
      }
      double uvunc = uv.uncertainty(), uv2unc = uv2.uncertainty();
      b = ::fabs(uvunc - uv2unc) < delta;
      if (!b) {
         printf("sigmas are different: %lf, %lf\n", uvunc, uv2unc);
      }
   }

   UncertainValue2Test::UncertainValue2Test()
   {
   }

   void UncertainValue2Test::testSpecialValues()
   {
      auto one = UncertainValue2::ONE();
      auto oneU = one.uncertainty();
      assertEquals(one.doubleValue(), 1.0, 1e-10);
      assertEquals(oneU, 0.0, 1e-10);

      auto nan = UncertainValue2::NaN();
      auto nanU = nan.uncertainty();
      assertEquals(nan.doubleValue(), NAN, 1e-10);
      assertEquals(nanU, 0.0, 1e-10);

      auto ni = UncertainValue2::NEGATIVE_INFINITY();
      auto niU = ni.uncertainty();
      assertEquals(ni.doubleValue(), -INFINITY, 1e-10);
      assertEquals(niU, 0.0, 1e-10);

      auto pi = UncertainValue2::POSITIVE_INFINITY();
      auto piU = pi.uncertainty();
      assertEquals(pi.doubleValue(), INFINITY, 1e-10);
      assertEquals(piU, 0.0, 1e-10);

      auto zero = UncertainValue2::ZERO();
      auto zeroU = zero.uncertainty();
      assertEquals(zero.doubleValue(), 0.0, 1e-10);
      assertEquals(zeroU, 0.0, 1e-10);

      printf("%s completed.\n", "UncertainValue2Test::testSpecialValues()");
   }

   void UncertainValue2Test::testA()
   {
      UncertainValue2::UncertainValue2 mA(1.24, "A", 0.3);
      UncertainValue2::UncertainValue2 mA2a = makeA2a();

      assertEquals(mA, mA2a, 1.0e-10);
      assertEquals(mA.uncertainty(), mA2a.uncertainty(), 1.0e-8);
      assertEquals(mA, mA2a, 1.0e-8);

      printf("%s completed.\n", "UncertainValue2Test::testA()");
   }

   void UncertainValue2Test::testB()
   {
      UncertainValue2::UncertainValue2 mB(8.82, "B", 1.2);
      UncertainValue2::UncertainValue2 mB2a = makeB2a();

      assertEquals(mB, mB2a, 1.0e-10);
      assertEquals(mB.uncertainty(), mB2a.uncertainty(), 1.0e-8);
      assertEquals(mB, mB2a, 1.0e-8);
      printf("%s completed.\n", "UncertainValue2Test::testB()");
   }

   void UncertainValue2Test::testC()
   {
      UncertainValue2::UncertainValue2 mC(-9.3, "C", 2.1);
      UncertainValue2::UncertainValue2 mC2a = makeC2a();

      assertEquals(mC, mC2a, 1.0e-10);
      assertEquals(mC.uncertainty(), mC2a.uncertainty(), 1.0e-8);
      assertEquals(mC, mC2a, 1.0e-8);
      printf("%s completed.\n", "UncertainValue2Test::testC()");
   }

   void UncertainValue2Test::testAB()
   {
      UncertainValue2::UncertainValue2 mA(1.24, "A", 0.3);
      UncertainValue2::UncertainValue2 mA2a(makeA2a());
      UncertainValue2::UncertainValue2 mB(8.82, "B", 1.2);
      UncertainValue2::UncertainValue2 mB2a(makeB2a());

      UncertainValue2::UncertainValue2 mA2a1(makeA2a());
      UncertainValue2::UncertainValue2 mA2a2(makeA2a());
      assertEquals(mA2a1.hashCode(), mA2a2.hashCode(), 0.0001);

      //auto pp = mA.getComponents();
      //printf("abc\n");
      //printf("mA comp: %d\n", pp.Size());
      //printf("??\n");
      //printf("mB comp: %d\n", mB.getComponents().Size());
      //printf("mA2a comp: %d\n", mA2a.getComponents().Size());
      //printf("mB2a comp: %d\n", mB2a.getComponents().Size());

      assertEquals(mA, mA2a, 1.0e-10);
      assertEquals(mB, mB2a, 1.0e-10);

      assertEquals(UncertainValue2::multiply(2.0, mA), UncertainValue2::multiply(2.0, mA2a), 1.0e-10);
      assertEquals(UncertainValue2::multiply(4.0, mB), UncertainValue2::multiply(4.0, mB2a), 1.0e-10);
      assertEquals(UncertainValue2::multiply(UncertainValue2::UncertainValue2(4.0), mB), UncertainValue2::multiply(4.0, mB2a), 1.0e-10);
      assertEquals(UncertainValue2::multiply(4.0, mB), UncertainValue2::multiply(4.0, mB2a), 1.0e-10);
      assertEquals(UncertainValue2::multiply(4.0, mB), UncertainValue2::multiply(UncertainValue2::UncertainValue2(4.0), mB2a), 1.0e-10);
      assertEquals(UncertainValue2::multiply(2.0, UncertainValue2::multiply(2.0, mB)), UncertainValue2::multiply(4.0, mB2a), 1.0e-10);
      assertEquals(UncertainValue2::multiply(4.0, mB), UncertainValue2::multiply(UncertainValue2::UncertainValue2(4.0), mB2a), 1.0e-10);

      assertEquals(UncertainValue2::multiply(mA, mB), UncertainValue2::UncertainValue2(10.9368, 3.035697613), 1.0e-8);

      assertEquals(UncertainValue2::divide(mA, mB), UncertainValue2::UncertainValue2(0.1405895692, 0.03902306155), 1.0e-8);

      assertEquals(UncertainValue2::add(mA, mB), UncertainValue2::UncertainValue2(10.06, 1.236931688), 1.0e-8);

      assertEquals(UncertainValue2::subtract(mA, mB), UncertainValue2::UncertainValue2(-7.58, 1.236931688), 1.0e-8);

      assertEquals(UncertainValue2::divide(1.0, mB), UncertainValue2::invert(mB2a), 1.0e-6);
      assertEquals(UncertainValue2::divide(1.0, mB), UncertainValue2::divide(UncertainValue2::UncertainValue2(1.0), mB2a), 1.0e-6);

      assertEquals(UncertainValue2::exp(mA), UncertainValue2::UncertainValue2(3.455613465, 1.036684039), 1.0e-6);
      assertEquals(UncertainValue2::exp(mA), UncertainValue2::exp(mA2a), 1.0e-6);

      assertEquals(UncertainValue2::log(mA), UncertainValue2::UncertainValue2(0.2151113796, 0.2419354839), 1.0e-6);
      assertEquals(UncertainValue2::log(mA), UncertainValue2::log(mA2a), 1.0e-6);

      assertEquals(UncertainValue2::pow(mA, 2.5), UncertainValue2::UncertainValue2(1.712198897, 1.035604171), 1.0e-6);
      assertEquals(UncertainValue2::pow(mA, 2.5), UncertainValue2::pow(mA2a, 2.5), 1.0e-6);

      assertEquals(mA.sqrt(), UncertainValue2::UncertainValue2(1.113552873, 0.1347039765), 1.0e-6);
      assertEquals(mA.sqrt(), mA2a.sqrt(), 1.0e-6);

      assertEquals(UncertainValue2::atan(mA), UncertainValue2::UncertainValue2(0.892133836, 0.118221942), 1.0e-6);
      assertEquals(UncertainValue2::atan(mA), UncertainValue2::atan(mA2a), 1.0e-6);

      assertEquals(UncertainValue2::sqr(mA), UncertainValue2::UncertainValue2(1.5376, 0.744), 1.0e-6);
      assertEquals(UncertainValue2::sqr(mA), UncertainValue2::sqr(mA2a), 1.0e-6);

      printf("%s completed.\n", "UncertainValue2Test::testAB()");
   }

   void UncertainValue2Test::testAdd1()
   {
      UncertainValue2::UncertainValue2 a(1.0, "A", 0.1);
      UncertainValue2::UncertainValue2 b(2.0, "A", 0.2);

      assertEquals(UncertainValue2::add(a, b), UncertainValue2::UncertainValue2(3.0, 0.3), 0.0001);
      assertEquals(UncertainValue2::add(1.0, a, 1.0, b), UncertainValue2::UncertainValue2(3.0, 0.3), 0.0001);
      assertEquals(UncertainValue2::add(2.0, a, 1.0, b), UncertainValue2::UncertainValue2(4.0, 0.4), 0.0001);
      assertEquals(UncertainValue2::add(1.0, a, 2.0, b), UncertainValue2::UncertainValue2(5.0, 0.5), 0.0001);
      assertEquals(UncertainValue2::add(1.0, a, -1.0, b), UncertainValue2::UncertainValue2(-1.0, 0.1), 0.0001);
      assertEquals(UncertainValue2::add(-1.0, a, 1.0, b), UncertainValue2::UncertainValue2(1.0, 0.1), 0.0001);
      assertEquals(UncertainValue2::add(-2.0, a, 1.0, b), UncertainValue2::UncertainValue2(0.0, 0.0), 0.0001);

      // Example page 22 GUM
      UncertainValue2::UncertainValue2 rsc[] = {
         UncertainValue2::UncertainValue2(1000.0, "R", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R", 0.1),
      };
      size_t n = sizeof(rsc) / sizeof(rsc[0]);

      //LinkedList::Node<UncertainValue2::UncertainValue2>* rsc = NULL;
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsc, UncertainValue2::UncertainValue2(1000.0, "R", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsc, UncertainValue2::UncertainValue2(1000.0, "R", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsc, UncertainValue2::UncertainValue2(1000.0, "R", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsc, UncertainValue2::UncertainValue2(1000.0, "R", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsc, UncertainValue2::UncertainValue2(1000.0, "R", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsc, UncertainValue2::UncertainValue2(1000.0, "R", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsc, UncertainValue2::UncertainValue2(1000.0, "R", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsc, UncertainValue2::UncertainValue2(1000.0, "R", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsc, UncertainValue2::UncertainValue2(1000.0, "R", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsc, UncertainValue2::UncertainValue2(1000.0, "R", 0.1));

      assertEquals(UncertainValue2::add(rsc, n).uncertainty(), 1.0, 1.0e-10);
      //LinkedList::RemoveAll(&rsc);
      printf("%s completed.\n", "UncertainValue2Test::testAdd1()");
   }

   void UncertainValue2Test::testAdd2()
   {
      UncertainValue2::UncertainValue2 a(1.0, "A", 0.1);
      UncertainValue2::UncertainValue2 b(2.0, "B", 0.2);
      assertEquals(UncertainValue2::add(a, b), UncertainValue2::UncertainValue2(3.0, sqrt(0.05)), 0.0001);
      assertEquals(UncertainValue2::add(1.0, a, 1.0, b), UncertainValue2::UncertainValue2(3.0, sqrt(0.05)), 0.0001);
      assertEquals(UncertainValue2::add(2.0, a, 1.0, b), UncertainValue2::UncertainValue2(4.0, sqrt(0.08)), 0.0001);
      assertEquals(UncertainValue2::add(1.0, a, 2.0, b), UncertainValue2::UncertainValue2(5.0, sqrt(0.17)), 0.0001);
      assertEquals(UncertainValue2::add(1.0, a, -1.0, b), UncertainValue2::UncertainValue2(-1.0, sqrt(0.05)), 0.0001);
      assertEquals(UncertainValue2::add(-1.0, a, 1.0, b), UncertainValue2::UncertainValue2(1.0, sqrt(0.05)), 0.0001);
      assertEquals(UncertainValue2::add(-2.0, a, 1.0, b), UncertainValue2::UncertainValue2(0.0, sqrt(0.08)), 0.0001);
      // Example page 22 GUM
      UncertainValue2::UncertainValue2 rsu[] = {
         UncertainValue2::UncertainValue2(1000.0, "R0", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R1", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R2", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R3", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R4", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R5", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R6", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R7", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R8", 0.1),
            UncertainValue2::UncertainValue2(1000.0, "R9", 0.1),
      };
      size_t n = sizeof(rsu) / sizeof(rsu[0]);

      //LinkedList::Node<UncertainValue2::UncertainValue2>* rsu = NULL;
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsu, UncertainValue2::UncertainValue2(1000.0, "R0", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsu, UncertainValue2::UncertainValue2(1000.0, "R1", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsu, UncertainValue2::UncertainValue2(1000.0, "R2", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsu, UncertainValue2::UncertainValue2(1000.0, "R3", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsu, UncertainValue2::UncertainValue2(1000.0, "R4", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsu, UncertainValue2::UncertainValue2(1000.0, "R5", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsu, UncertainValue2::UncertainValue2(1000.0, "R6", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsu, UncertainValue2::UncertainValue2(1000.0, "R7", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsu, UncertainValue2::UncertainValue2(1000.0, "R8", 0.1));
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&rsu, UncertainValue2::UncertainValue2(1000.0, "R9", 0.1));

      assertEquals(UncertainValue2::add(rsu, n).uncertainty(), 0.32, 0.005);
      //LinkedList::RemoveAll(&rsu);
      printf("%s completed.\n", "UncertainValue2Test::testAdd2()");
   }

   void UncertainValue2Test::testAdd3()
   {
      UncertainValue2::UncertainValue2 a(1.0, "A", 0.1);
      UncertainValue2::UncertainValue2 b(2.0, "B", 0.2);
      UncertainValue2::UncertainValue2 c(3.0, "A", 0.15);

      UncertainValue2::UncertainValue2 uvs[] = {
         a,
         b,
         c
      };
      size_t n = sizeof(uvs) / sizeof(uvs[0]);

      //LinkedList::Node<UncertainValue2::UncertainValue2>* uvs = NULL;
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&uvs, a);
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&uvs, b);
      //LinkedList::InsertHead<UncertainValue2::UncertainValue2>(&uvs, c);

      assertEquals(UncertainValue2::add(uvs, n), UncertainValue2::UncertainValue2(6.0, sqrt(0.0625 + 0.04)), 1e-6);
      //LinkedList::RemoveAll(&uvs);
      printf("%s completed.\n", "UncertainValue2Test::testAdd3()");
   }

   void UncertainValue2Test::testMultiply()
   {
      UncertainValue2::UncertainValue2 a(1.1, "A", 0.1);
      UncertainValue2::UncertainValue2 b(2.3, "B", 0.2);
      UncertainValue2::UncertainValue2 c(3.6, "A", 0.15);

      assertEquals(UncertainValue2::multiply(a, b), UncertainValue2::UncertainValue2(2.53, 0.3182766093), 1.0e-6);
      assertEquals(UncertainValue2::multiply(a, c), UncertainValue2::UncertainValue2(3.96, 0.525), 1.0e-6);
      assertEquals(UncertainValue2::multiply(b, UncertainValue2::multiply(a, c)), UncertainValue2::UncertainValue2(9.108, 1.444063797), 1.0e-6);
      printf("%s completed.\n", "UncertainValue2Test::testMultiply()");
   }

   void UncertainValue2Test::testDivide()
   {
      UncertainValue2::UncertainValue2 a(1.1, "A", 0.1);
      UncertainValue2::UncertainValue2 b(2.3, "B", 0.2);
      UncertainValue2::UncertainValue2 c(3.6, "A", 0.15);

      assertEquals(UncertainValue2::divide(b, a), UncertainValue2::UncertainValue2(2.090909091, 0.26303852), 1.0e-6);
      assertEquals(UncertainValue2::divide(b, c), UncertainValue2::UncertainValue2(0.6388888889, 0.06160408973), 1.0e-6);
      assertEquals(UncertainValue2::divide(a, c), UncertainValue2::UncertainValue2(0.3055555556, (0.1 / 1.1 + 0.15 / 3.6) * 0.3055555556), 1.0e-6);
      assertEquals(UncertainValue2::divide(2.0, c), UncertainValue2::divide(UncertainValue2::UncertainValue2(2.0), c), 1.0e-6);
      assertEquals(UncertainValue2::divide(-2.0, c), UncertainValue2::divide(UncertainValue2::UncertainValue2(-2.0), c), 1.0e-6);
      assertEquals(UncertainValue2::divide(c, 2.0), UncertainValue2::divide(c, UncertainValue2::UncertainValue2(2.0)), 1.0e-6);
      assertEquals(UncertainValue2::divide(c, -2.0), UncertainValue2::divide(c, UncertainValue2::UncertainValue2(-2.0)), 1.0e-6);
      assertEquals(UncertainValue2::divide(c, -2.0), UncertainValue2::multiply(c, UncertainValue2::UncertainValue2(-0.5)), 1.0e-6);
      assertEquals(UncertainValue2::multiply(b, UncertainValue2::invert(a)), UncertainValue2::divide(b, a), 1.0e-7);
      assertEquals(UncertainValue2::multiply(c, UncertainValue2::invert(a)), UncertainValue2::divide(c, a), 1.0e-7);
      printf("%s completed.\n", "UncertainValue2Test::testDivide()");
   }

   void UncertainValue2Test::testFunctions()
   {
      UncertainValue2::UncertainValue2 a(1.1, "A", 0.1);
      UncertainValue2::UncertainValue2 b(2.3, "B", 0.2);
      UncertainValue2::UncertainValue2 c(-3.6, "A", 0.15);

      assertEquals(UncertainValue2::exp(a), UncertainValue2::UncertainValue2(::exp(1.1), ::exp(1.1) * 0.1), 1.0e-8);
      assertEquals(UncertainValue2::exp(b), UncertainValue2::UncertainValue2(::exp(2.3), ::exp(2.3) * 0.2), 1.0e-8);
      assertEquals(UncertainValue2::pow(a, 3), UncertainValue2::multiply(UncertainValue2::multiply(a, a), a), 1.0e-8);
      assertEquals(UncertainValue2::sqrt(a), UncertainValue2::UncertainValue2(1.048808848, 0.04767312946), 1.0e-8);
      assertEquals(UncertainValue2::pow(a, 0.5), UncertainValue2::UncertainValue2(1.048808848, 0.04767312946), 1.0e-8);
      assertEquals(UncertainValue2::pow(a, -0.5), UncertainValue2::UncertainValue2(0.9534625894, 0.04333920861), 1.0e-8);
      assertEquals(UncertainValue2::pow(c, 2.0), UncertainValue2::multiply(c, c), 1.0e-8);
      assertEquals(UncertainValue2::atan(a), UncertainValue2::UncertainValue2(0.8329812667, 0.04524886878), 1.0e-8);
      assertEquals(UncertainValue2::atan2(b, a), UncertainValue2::atan(UncertainValue2::divide(b, a)), 1.0e-8);
      printf("%s completed.\n", "UncertainValue2Test::testFunctions()");
   }
}
