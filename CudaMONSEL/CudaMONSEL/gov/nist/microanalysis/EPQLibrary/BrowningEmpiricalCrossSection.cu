#include "gov\nist\microanalysis\EPQLibrary\BrowningEmpiricalCrossSection.cuh"

#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\FromSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

#include "gov\nist\microanalysis\Utility\Math2.cuh"

namespace BrowningEmpiricalCrossSection
{
   BrowningEmpiricalCrossSection::BrowningEmpiricalCrossSection(const ElementT& elm) :
      mElement(elm),
      mZp17(::pow(elm.getAtomicNumber(), 1.7)),
      mZp2(::pow(elm.getAtomicNumber(), 2.0)),
      mZp3(::pow(elm.getAtomicNumber(), 3.0))
   {
   }

   const ElementT& BrowningEmpiricalCrossSection::getElement() const
   {
      return mElement;
   }

   double BrowningEmpiricalCrossSection::totalCrossSection(double energy) const
   {
      const double e = FromSI::keV(energy);
      const double re = ::sqrt(e);
      return 3.0e-22 * mZp17 / (e + 0.005 * mZp17 * re + 0.0007 * mZp2 / re);
   }

   double BrowningEmpiricalCrossSection::randomScatteringAngle(double energy) const
   {
      const double r1 = Math2::random(), r2 = Math2::random();
      const double z = mElement.getAtomicNumber();
      const double e = FromSI::keV(energy);
      const double r = (300.0 * e / z) + (mZp3 / (3.0e5 * e));
      if (r1 <= r / (r + 1)) {
         // Screened Rutherford scattering
         const double alpha = 7.0e-3 / e;
         return ::acos(1.0 - ((2.0 * alpha * r2) / (alpha - r2 + 1)));
      }
      else
         // Isotropic scattering
         return ::acos(1 - 2.0 * r2);
   }

   const BrowningEmpiricalCrossSection BECS1(1);
   const BrowningEmpiricalCrossSection BECS2(2);
   const BrowningEmpiricalCrossSection BECS3(3);
   const BrowningEmpiricalCrossSection BECS4(4);
   const BrowningEmpiricalCrossSection BECS5(5);
   const BrowningEmpiricalCrossSection BECS6(6);
   const BrowningEmpiricalCrossSection BECS7(7);
   const BrowningEmpiricalCrossSection BECS8(8);
   const BrowningEmpiricalCrossSection BECS9(9);
   const BrowningEmpiricalCrossSection BECS10(10);
   const BrowningEmpiricalCrossSection BECS11(11);
   const BrowningEmpiricalCrossSection BECS12(12);
   const BrowningEmpiricalCrossSection BECS13(13);
   const BrowningEmpiricalCrossSection BECS14(14);
   const BrowningEmpiricalCrossSection BECS15(15);
   const BrowningEmpiricalCrossSection BECS16(16);
   const BrowningEmpiricalCrossSection BECS17(17);
   const BrowningEmpiricalCrossSection BECS18(18);
   const BrowningEmpiricalCrossSection BECS19(19);
   const BrowningEmpiricalCrossSection BECS20(20);
   const BrowningEmpiricalCrossSection BECS21(21);
   const BrowningEmpiricalCrossSection BECS22(22);
   const BrowningEmpiricalCrossSection BECS23(23);
   const BrowningEmpiricalCrossSection BECS24(24);
   const BrowningEmpiricalCrossSection BECS25(25);
   const BrowningEmpiricalCrossSection BECS26(26);
   const BrowningEmpiricalCrossSection BECS27(27);
   const BrowningEmpiricalCrossSection BECS28(28);
   const BrowningEmpiricalCrossSection BECS29(29);
   const BrowningEmpiricalCrossSection BECS30(30);
   const BrowningEmpiricalCrossSection BECS31(31);
   const BrowningEmpiricalCrossSection BECS32(32);
   const BrowningEmpiricalCrossSection BECS33(33);
   const BrowningEmpiricalCrossSection BECS34(34);
   const BrowningEmpiricalCrossSection BECS35(35);
   const BrowningEmpiricalCrossSection BECS36(36);
   const BrowningEmpiricalCrossSection BECS37(37);
   const BrowningEmpiricalCrossSection BECS38(38);
   const BrowningEmpiricalCrossSection BECS39(39);
   const BrowningEmpiricalCrossSection BECS40(40);
   const BrowningEmpiricalCrossSection BECS41(41);
   const BrowningEmpiricalCrossSection BECS42(42);
   const BrowningEmpiricalCrossSection BECS43(43);
   const BrowningEmpiricalCrossSection BECS44(44);
   const BrowningEmpiricalCrossSection BECS45(45);
   const BrowningEmpiricalCrossSection BECS46(46);
   const BrowningEmpiricalCrossSection BECS47(47);
   const BrowningEmpiricalCrossSection BECS48(48);
   const BrowningEmpiricalCrossSection BECS49(49);
   const BrowningEmpiricalCrossSection BECS50(50);
   const BrowningEmpiricalCrossSection BECS51(51);
   const BrowningEmpiricalCrossSection BECS52(52);
   const BrowningEmpiricalCrossSection BECS53(53);
   const BrowningEmpiricalCrossSection BECS54(54);
   const BrowningEmpiricalCrossSection BECS55(55);
   const BrowningEmpiricalCrossSection BECS56(56);
   const BrowningEmpiricalCrossSection BECS57(57);
   const BrowningEmpiricalCrossSection BECS58(58);
   const BrowningEmpiricalCrossSection BECS59(59);
   const BrowningEmpiricalCrossSection BECS60(60);
   const BrowningEmpiricalCrossSection BECS61(61);
   const BrowningEmpiricalCrossSection BECS62(62);
   const BrowningEmpiricalCrossSection BECS63(63);
   const BrowningEmpiricalCrossSection BECS64(64);
   const BrowningEmpiricalCrossSection BECS65(65);
   const BrowningEmpiricalCrossSection BECS66(66);
   const BrowningEmpiricalCrossSection BECS67(67);
   const BrowningEmpiricalCrossSection BECS68(68);
   const BrowningEmpiricalCrossSection BECS69(69);
   const BrowningEmpiricalCrossSection BECS70(70);
   const BrowningEmpiricalCrossSection BECS71(71);
   const BrowningEmpiricalCrossSection BECS72(72);
   const BrowningEmpiricalCrossSection BECS73(73);
   const BrowningEmpiricalCrossSection BECS74(74);
   const BrowningEmpiricalCrossSection BECS75(75);
   const BrowningEmpiricalCrossSection BECS76(76);
   const BrowningEmpiricalCrossSection BECS77(77);
   const BrowningEmpiricalCrossSection BECS78(78);
   const BrowningEmpiricalCrossSection BECS79(79);
   const BrowningEmpiricalCrossSection BECS80(80);
   const BrowningEmpiricalCrossSection BECS81(81);
   const BrowningEmpiricalCrossSection BECS82(82);
   const BrowningEmpiricalCrossSection BECS83(83);
   const BrowningEmpiricalCrossSection BECS84(84);
   const BrowningEmpiricalCrossSection BECS85(85);
   const BrowningEmpiricalCrossSection BECS86(86);
   const BrowningEmpiricalCrossSection BECS87(87);
   const BrowningEmpiricalCrossSection BECS88(88);
   const BrowningEmpiricalCrossSection BECS89(89);
   const BrowningEmpiricalCrossSection BECS90(90);
   const BrowningEmpiricalCrossSection BECS91(91);
   const BrowningEmpiricalCrossSection BECS92(92);
   const BrowningEmpiricalCrossSection BECS93(93);
   const BrowningEmpiricalCrossSection BECS94(94);
   const BrowningEmpiricalCrossSection BECS95(95);
   const BrowningEmpiricalCrossSection BECS96(96);

   BrowningEmpiricalCrossSection const * mScatter[113] = {
      nullptr,
      &BECS1,
      &BECS2,
      &BECS3,
      &BECS4,
      &BECS5,
      &BECS6,
      &BECS7,
      &BECS8,
      &BECS9,
      &BECS10,
      &BECS11,
      &BECS12,
      &BECS13,
      &BECS14,
      &BECS15,
      &BECS16,
      &BECS17,
      &BECS18,
      &BECS19,
      &BECS20,
      &BECS21,
      &BECS22,
      &BECS23,
      &BECS24,
      &BECS25,
      &BECS26,
      &BECS27,
      &BECS28,
      &BECS29,
      &BECS30,
      &BECS31,
      &BECS32,
      &BECS33,
      &BECS34,
      &BECS35,
      &BECS36,
      &BECS37,
      &BECS38,
      &BECS39,
      &BECS40,
      &BECS41,
      &BECS42,
      &BECS43,
      &BECS44,
      &BECS45,
      &BECS46,
      &BECS47,
      &BECS48,
      &BECS49,
      &BECS50,
      &BECS51,
      &BECS52,
      &BECS53,
      &BECS54,
      &BECS55,
      &BECS56,
      &BECS57,
      &BECS58,
      &BECS59,
      &BECS60,
      &BECS61,
      &BECS62,
      &BECS63,
      &BECS64,
      &BECS65,
      &BECS66,
      &BECS67,
      &BECS68,
      &BECS69,
      &BECS70,
      &BECS71,
      &BECS72,
      &BECS73,
      &BECS74,
      &BECS75,
      &BECS76,
      &BECS77,
      &BECS78,
      &BECS79,
      &BECS80,
      &BECS81,
      &BECS82,
      &BECS83,
      &BECS84,
      &BECS85,
      &BECS86,
      &BECS87,
      &BECS88,
      &BECS89,
      &BECS90,
      &BECS91,
      &BECS92,
      &BECS93,
      &BECS94,
      &BECS95,
      &BECS96
   };

   const BrowningEmpiricalCrossSection& getBECS(int an)
   {
      return *mScatter[an];
   }
}
