#include "MaterialTest.cuh"

#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Composition.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"

namespace MaterialTest
{
   MaterialTest::MaterialTest()
   {
   }

   void assertEquals(double a, double b, double delta)
   {
      bool boo = false;
      boo = fabs(a - b) < delta;
      if (!boo) {
         printf("values are different: %lf, %lf\n", a, b);
      }
   }

   void assertTrue(bool expr)
   {
      if (!expr) {
         printf("false\n");
      }
   }

   void MaterialTest::testOne()
   {
      ::Composition::data_type eps = 1.0e-8f;
      {
         Composition::Composition mat;
         // Titanium dioxide - TiO2 (Z(Ti)=22,
         ::Composition::data_type wTi = Element::getAtomicWeight(Element::elmTi) / (Element::getAtomicWeight(Element::elmTi) + 2.0 * Element::getAtomicWeight(Element::elmO));
         ::Composition::data_type wO = 2.0 * Element::getAtomicWeight(Element::elmO) / (Element::getAtomicWeight(Element::elmTi) + 2.0 * Element::getAtomicWeight(Element::elmO));
         assertEquals(wTi + wO, 1.0, eps);
         ::Composition::data_type fac = 1.1f;
         mat.addElement(Element::elmTi, fac * wTi);
         mat.addElement(Element::elmO, fac * wO);
         {
            int elemLen = 2;
            const Element::Element* elms[] = {
               &Element::Ti,
               &Element::O
            };
            ::Composition::data_type wgts[] = {
               fac * wTi,
               fac * wO
            };

            Composition::Composition dup(elms, elemLen, wgts, elemLen);

            //printf("mat \n");
            //for (auto c : mat.GetConstituents()) {
            //   printf("%d: ", c.first.getAtomicNumber());
            //   auto uv = c.second;
            //   printf("%lf - ", uv.doubleValue());
            //   for (auto s : uv.getComponents()) {
            //      printf("(%s, %lf)", s.first.c_str(), s.second);
            //   }
            //   printf("\n");
            //}
            //printf("dup \n");
            //for (auto c : dup.GetConstituents()) {
            //   printf("%d: ", c.first.getAtomicNumber());
            //   auto uv = c.second;
            //   printf("%lf - ", uv.doubleValue());
            //   for (auto s : uv.getComponents()) {
            //      printf("(%s, %lf)", s.first.c_str(), s.second);
            //   }
            //   printf("\n");
            //}

            assertTrue(mat.equals(dup));
         }
         assertEquals(mat.weightFraction(Element::Ti, true), wTi, eps);
         assertEquals(mat.weightFraction(Element::O, true), wO, eps);
         assertEquals(mat.weightFraction(Element::Ti, false), fac * wTi, eps);
         assertEquals(mat.weightFraction(Element::O, false), fac * wO, eps);
         assertTrue(mat.equals(mat));
         {
            Material::Material mm(1.0);
            mm.addElement(Element::elmW, 2.0);
            mm.addElement(Element::elmO, 1.0);
            assertTrue(!mat.equals(mm));
         }
      }
      { // K3189 - NIST SRM glass
         int elemLen = 7;
         Element::Element const * elms[] = {
            &Element::Si,
            &Element::O,
            &Element::Al,
            &Element::Ca,
            &Element::Mg,
            &Element::Ti,
            &Element::Fe
         };
         ::Composition::data_type mF[] = {
            0.151971599,
            0.608811816,
            0.062690151,
            0.05699065,
            0.056638401,
            0.005716531,
            0.057180851
         };
         ::Composition::data_type wF[] = {
            0.186973968,
            0.426700667,
            0.074093109,
            0.100056707,
            0.06030359,
            0.011986858,
            0.139885101
         };
         Material::Material mat(1.0);
         mat.defineByMoleFraction(elms, elemLen, mF, elemLen);
         for (int i = 0; i < elemLen; ++i) {
            assertEquals(mat.weightFraction(*elms[i], true), wF[i], 1.0e-5);
            //printf("%d, %f\n", elms[i].getAtomicNumber(), mF[i]);
            assertEquals(mat.atomicPercent(*elms[i]), mF[i], 1.0e-5);
         }
         { // Try an alternative method to define a K3189
            Material::Material mat0(1.0);
            const Element::Element* elms0[] = { &Element::Si, &Element::O };
            ::Composition::data_type massFracs0[] = { 1.0, 2.0 };
            mat0.defineByMoleFraction(elms0, 2, massFracs0, 2);

            Material::Material mat1(1.0);
            const Element::Element* elms1[] = { &Element::Al, &Element::O };
            ::Composition::data_type massFracs1[] { 2.0, 3.0 };
            mat1.defineByMoleFraction(elms1, 2, massFracs1, 2);

            Material::Material mat2(1.0);
            const Element::Element* elms2[] { &Element::Ca, &Element::O };
            ::Composition::data_type massFracs2[] { 1.0, 1.0 };
            mat2.defineByMoleFraction(elms2, 2, massFracs2, 2);

            Material::Material mat3(1.0);
            const Element::Element* elm3[] = { &Element::Mg, &Element::O };
            ::Composition::data_type massFracs3[] = { 1.0, 1.0 };
            mat3.defineByMoleFraction(elm3, 2, massFracs3, 2);

            Material::Material mat4(1.0);
            const Element::Element* elm4[] = { &Element::Ti, &Element::O };
            ::Composition::data_type massFracs4[] = { 1.0, 2.0 };
            mat4.defineByMoleFraction(elm4, 2, massFracs4, 2);

            Material::Material mat5(1.0);
            const Element::Element* elm5[] = { &Element::Fe, &Element::O };
            ::Composition::data_type massFracs5[] = { 2.0, 3.0 };
            mat5.defineByMoleFraction(elm5, 2, massFracs5, 2);

            const Material::Material* mats[6] = { &mat0, &mat1, &mat2, &mat3, &mat4, &mat5 };
            //printf("%d\n", mat0.GetConstituents().size());
            //printf("%d\n", mat1.GetConstituents().size());
            //printf("%d\n", mat2.GetConstituents().size());
            //printf("%d\n", mat3.GetConstituents().size());
            //printf("%d\n", mat4.GetConstituents().size());
            //printf("%d\n", mat5.GetConstituents().size());

            Material::Material mat(1.0);
            ::Composition::data_type massFracs[] = {
               0.40,
               0.14,
               0.14,
               0.10,
               0.02,
               0.20
            };
            mat.defineByMaterialFraction(mats, 6, massFracs, 6);
            for (int i = 0; i < elemLen; ++i) {
               assertEquals(mat.weightFraction(*elms[i], true), wF[i], 1.0e-5);
            }
            Material::Material dup = mat.clone();
            assertTrue(mat == dup);
            dup.addElement(Element::elmAr, 0.12);
            assertTrue(!mat.equals(dup));
         }
      }

      //{
      //   String[] mats = {
      //      MaterialFactory.K2450,
      //      MaterialFactory.Ice,
      //      MaterialFactory.Al2O3
      //   };
      //   for (int i = 0; i < mats.length; ++i) {
      //      Composition::Composition comp = MaterialFactory.createMaterial(mats[i]);
      //      String ms = EPQXStream.getInstance().toXML(comp);
      //      if (false) {
      //         System.out.print(comp);
      //         System.out.print("\t");
      //         System.out.print(ms);
      //         System.out.println();
      //      }
      //      assertTrue(EPQXStream.getInstance().fromXML(ms).equals(comp));
      //   }
      //}
      printf("MaterialTest::testOne() completed\n");
   }
}
