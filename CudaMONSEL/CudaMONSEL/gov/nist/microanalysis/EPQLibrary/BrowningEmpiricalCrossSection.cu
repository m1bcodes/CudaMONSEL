#include "gov\nist\microanalysis\EPQLibrary\BrowningEmpiricalCrossSection.cuh"

#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\FromSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

#include "gov\nist\microanalysis\Utility\Math2.cuh"

namespace BrowningEmpiricalCrossSection
{
   __host__ __device__ BrowningEmpiricalCrossSection::BrowningEmpiricalCrossSection(const ElementT& elm) :
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

   __device__ const BrowningEmpiricalCrossSection *dBECS1;
   __device__ const BrowningEmpiricalCrossSection *dBECS2;
   __device__ const BrowningEmpiricalCrossSection *dBECS3;
   __device__ const BrowningEmpiricalCrossSection *dBECS4;
   __device__ const BrowningEmpiricalCrossSection *dBECS5;
   __device__ const BrowningEmpiricalCrossSection *dBECS6;
   __device__ const BrowningEmpiricalCrossSection *dBECS7;
   __device__ const BrowningEmpiricalCrossSection *dBECS8;
   __device__ const BrowningEmpiricalCrossSection *dBECS9;
   __device__ const BrowningEmpiricalCrossSection *dBECS10;
   __device__ const BrowningEmpiricalCrossSection *dBECS11;
   __device__ const BrowningEmpiricalCrossSection *dBECS12;
   __device__ const BrowningEmpiricalCrossSection *dBECS13;
   __device__ const BrowningEmpiricalCrossSection *dBECS14;
   __device__ const BrowningEmpiricalCrossSection *dBECS15;
   __device__ const BrowningEmpiricalCrossSection *dBECS16;
   __device__ const BrowningEmpiricalCrossSection *dBECS17;
   __device__ const BrowningEmpiricalCrossSection *dBECS18;
   __device__ const BrowningEmpiricalCrossSection *dBECS19;
   __device__ const BrowningEmpiricalCrossSection *dBECS20;
   __device__ const BrowningEmpiricalCrossSection *dBECS21;
   __device__ const BrowningEmpiricalCrossSection *dBECS22;
   __device__ const BrowningEmpiricalCrossSection *dBECS23;
   __device__ const BrowningEmpiricalCrossSection *dBECS24;
   __device__ const BrowningEmpiricalCrossSection *dBECS25;
   __device__ const BrowningEmpiricalCrossSection *dBECS26;
   __device__ const BrowningEmpiricalCrossSection *dBECS27;
   __device__ const BrowningEmpiricalCrossSection *dBECS28;
   __device__ const BrowningEmpiricalCrossSection *dBECS29;
   __device__ const BrowningEmpiricalCrossSection *dBECS30;
   __device__ const BrowningEmpiricalCrossSection *dBECS31;
   __device__ const BrowningEmpiricalCrossSection *dBECS32;
   __device__ const BrowningEmpiricalCrossSection *dBECS33;
   __device__ const BrowningEmpiricalCrossSection *dBECS34;
   __device__ const BrowningEmpiricalCrossSection *dBECS35;
   __device__ const BrowningEmpiricalCrossSection *dBECS36;
   __device__ const BrowningEmpiricalCrossSection *dBECS37;
   __device__ const BrowningEmpiricalCrossSection *dBECS38;
   __device__ const BrowningEmpiricalCrossSection *dBECS39;
   __device__ const BrowningEmpiricalCrossSection *dBECS40;
   __device__ const BrowningEmpiricalCrossSection *dBECS41;
   __device__ const BrowningEmpiricalCrossSection *dBECS42;
   __device__ const BrowningEmpiricalCrossSection *dBECS43;
   __device__ const BrowningEmpiricalCrossSection *dBECS44;
   __device__ const BrowningEmpiricalCrossSection *dBECS45;
   __device__ const BrowningEmpiricalCrossSection *dBECS46;
   __device__ const BrowningEmpiricalCrossSection *dBECS47;
   __device__ const BrowningEmpiricalCrossSection *dBECS48;
   __device__ const BrowningEmpiricalCrossSection *dBECS49;
   __device__ const BrowningEmpiricalCrossSection *dBECS50;
   __device__ const BrowningEmpiricalCrossSection *dBECS51;
   __device__ const BrowningEmpiricalCrossSection *dBECS52;
   __device__ const BrowningEmpiricalCrossSection *dBECS53;
   __device__ const BrowningEmpiricalCrossSection *dBECS54;
   __device__ const BrowningEmpiricalCrossSection *dBECS55;
   __device__ const BrowningEmpiricalCrossSection *dBECS56;
   __device__ const BrowningEmpiricalCrossSection *dBECS57;
   __device__ const BrowningEmpiricalCrossSection *dBECS58;
   __device__ const BrowningEmpiricalCrossSection *dBECS59;
   __device__ const BrowningEmpiricalCrossSection *dBECS60;
   __device__ const BrowningEmpiricalCrossSection *dBECS61;
   __device__ const BrowningEmpiricalCrossSection *dBECS62;
   __device__ const BrowningEmpiricalCrossSection *dBECS63;
   __device__ const BrowningEmpiricalCrossSection *dBECS64;
   __device__ const BrowningEmpiricalCrossSection *dBECS65;
   __device__ const BrowningEmpiricalCrossSection *dBECS66;
   __device__ const BrowningEmpiricalCrossSection *dBECS67;
   __device__ const BrowningEmpiricalCrossSection *dBECS68;
   __device__ const BrowningEmpiricalCrossSection *dBECS69;
   __device__ const BrowningEmpiricalCrossSection *dBECS70;
   __device__ const BrowningEmpiricalCrossSection *dBECS71;
   __device__ const BrowningEmpiricalCrossSection *dBECS72;
   __device__ const BrowningEmpiricalCrossSection *dBECS73;
   __device__ const BrowningEmpiricalCrossSection *dBECS74;
   __device__ const BrowningEmpiricalCrossSection *dBECS75;
   __device__ const BrowningEmpiricalCrossSection *dBECS76;
   __device__ const BrowningEmpiricalCrossSection *dBECS77;
   __device__ const BrowningEmpiricalCrossSection *dBECS78;
   __device__ const BrowningEmpiricalCrossSection *dBECS79;
   __device__ const BrowningEmpiricalCrossSection *dBECS80;
   __device__ const BrowningEmpiricalCrossSection *dBECS81;
   __device__ const BrowningEmpiricalCrossSection *dBECS82;
   __device__ const BrowningEmpiricalCrossSection *dBECS83;
   __device__ const BrowningEmpiricalCrossSection *dBECS84;
   __device__ const BrowningEmpiricalCrossSection *dBECS85;
   __device__ const BrowningEmpiricalCrossSection *dBECS86;
   __device__ const BrowningEmpiricalCrossSection *dBECS87;
   __device__ const BrowningEmpiricalCrossSection *dBECS88;
   __device__ const BrowningEmpiricalCrossSection *dBECS89;
   __device__ const BrowningEmpiricalCrossSection *dBECS90;
   __device__ const BrowningEmpiricalCrossSection *dBECS91;
   __device__ const BrowningEmpiricalCrossSection *dBECS92;
   __device__ const BrowningEmpiricalCrossSection *dBECS93;
   __device__ const BrowningEmpiricalCrossSection *dBECS94;
   __device__ const BrowningEmpiricalCrossSection *dBECS95;
   __device__ const BrowningEmpiricalCrossSection *dBECS96;

   __device__ BrowningEmpiricalCrossSection const * dScatter[113];

   __global__ void initCuda()
   {
      dScatter[1] = new BrowningEmpiricalCrossSection(1);
      dScatter[2] = new BrowningEmpiricalCrossSection(2);
      dScatter[3] = new BrowningEmpiricalCrossSection(3);
      dScatter[4] = new BrowningEmpiricalCrossSection(4);
      dScatter[5] = new BrowningEmpiricalCrossSection(5);
      dScatter[6] = new BrowningEmpiricalCrossSection(6);
      dScatter[7] = new BrowningEmpiricalCrossSection(7);
      dScatter[8] = new BrowningEmpiricalCrossSection(8);
      dScatter[9] = new BrowningEmpiricalCrossSection(9);
      dScatter[10] = new BrowningEmpiricalCrossSection(10);
      dScatter[11] = new BrowningEmpiricalCrossSection(11);
      dScatter[12] = new BrowningEmpiricalCrossSection(12);
      dScatter[13] = new BrowningEmpiricalCrossSection(13);
      dScatter[14] = new BrowningEmpiricalCrossSection(14);
      dScatter[15] = new BrowningEmpiricalCrossSection(15);
      dScatter[16] = new BrowningEmpiricalCrossSection(16);
      dScatter[17] = new BrowningEmpiricalCrossSection(17);
      dScatter[18] = new BrowningEmpiricalCrossSection(18);
      dScatter[19] = new BrowningEmpiricalCrossSection(19);
      dScatter[20] = new BrowningEmpiricalCrossSection(20);
      dScatter[21] = new BrowningEmpiricalCrossSection(21);
      dScatter[22] = new BrowningEmpiricalCrossSection(22);
      dScatter[23] = new BrowningEmpiricalCrossSection(23);
      dScatter[24] = new BrowningEmpiricalCrossSection(24);
      dScatter[25] = new BrowningEmpiricalCrossSection(25);
      dScatter[26] = new BrowningEmpiricalCrossSection(26);
      dScatter[27] = new BrowningEmpiricalCrossSection(27);
      dScatter[28] = new BrowningEmpiricalCrossSection(28);
      dScatter[29] = new BrowningEmpiricalCrossSection(29);
      dScatter[30] = new BrowningEmpiricalCrossSection(30);
      dScatter[31] = new BrowningEmpiricalCrossSection(31);
      dScatter[32] = new BrowningEmpiricalCrossSection(32);
      dScatter[33] = new BrowningEmpiricalCrossSection(33);
      dScatter[34] = new BrowningEmpiricalCrossSection(34);
      dScatter[35] = new BrowningEmpiricalCrossSection(35);
      dScatter[36] = new BrowningEmpiricalCrossSection(36);
      dScatter[37] = new BrowningEmpiricalCrossSection(37);
      dScatter[38] = new BrowningEmpiricalCrossSection(38);
      dScatter[39] = new BrowningEmpiricalCrossSection(39);
      dScatter[40] = new BrowningEmpiricalCrossSection(40);
      dScatter[41] = new BrowningEmpiricalCrossSection(41);
      dScatter[42] = new BrowningEmpiricalCrossSection(42);
      dScatter[43] = new BrowningEmpiricalCrossSection(43);
      dScatter[44] = new BrowningEmpiricalCrossSection(44);
      dScatter[45] = new BrowningEmpiricalCrossSection(45);
      dScatter[46] = new BrowningEmpiricalCrossSection(46);
      dScatter[47] = new BrowningEmpiricalCrossSection(47);
      dScatter[48] = new BrowningEmpiricalCrossSection(48);
      dScatter[49] = new BrowningEmpiricalCrossSection(49);
      dScatter[50] = new BrowningEmpiricalCrossSection(50);
      dScatter[51] = new BrowningEmpiricalCrossSection(51);
      dScatter[52] = new BrowningEmpiricalCrossSection(52);
      dScatter[53] = new BrowningEmpiricalCrossSection(53);
      dScatter[54] = new BrowningEmpiricalCrossSection(54);
      dScatter[55] = new BrowningEmpiricalCrossSection(55);
      dScatter[56] = new BrowningEmpiricalCrossSection(56);
      dScatter[57] = new BrowningEmpiricalCrossSection(57);
      dScatter[58] = new BrowningEmpiricalCrossSection(58);
      dScatter[59] = new BrowningEmpiricalCrossSection(59);
      dScatter[60] = new BrowningEmpiricalCrossSection(60);
      dScatter[61] = new BrowningEmpiricalCrossSection(61);
      dScatter[62] = new BrowningEmpiricalCrossSection(62);
      dScatter[63] = new BrowningEmpiricalCrossSection(63);
      dScatter[64] = new BrowningEmpiricalCrossSection(64);
      dScatter[65] = new BrowningEmpiricalCrossSection(65);
      dScatter[66] = new BrowningEmpiricalCrossSection(66);
      dScatter[67] = new BrowningEmpiricalCrossSection(67);
      dScatter[68] = new BrowningEmpiricalCrossSection(68);
      dScatter[69] = new BrowningEmpiricalCrossSection(69);
      dScatter[70] = new BrowningEmpiricalCrossSection(70);
      dScatter[71] = new BrowningEmpiricalCrossSection(71);
      dScatter[72] = new BrowningEmpiricalCrossSection(72);
      dScatter[73] = new BrowningEmpiricalCrossSection(73);
      dScatter[74] = new BrowningEmpiricalCrossSection(74);
      dScatter[75] = new BrowningEmpiricalCrossSection(75);
      dScatter[76] = new BrowningEmpiricalCrossSection(76);
      dScatter[77] = new BrowningEmpiricalCrossSection(77);
      dScatter[78] = new BrowningEmpiricalCrossSection(78);
      dScatter[79] = new BrowningEmpiricalCrossSection(79);
      dScatter[80] = new BrowningEmpiricalCrossSection(80);
      dScatter[81] = new BrowningEmpiricalCrossSection(81);
      dScatter[82] = new BrowningEmpiricalCrossSection(82);
      dScatter[83] = new BrowningEmpiricalCrossSection(83);
      dScatter[84] = new BrowningEmpiricalCrossSection(84);
      dScatter[85] = new BrowningEmpiricalCrossSection(85);
      dScatter[86] = new BrowningEmpiricalCrossSection(86);
      dScatter[87] = new BrowningEmpiricalCrossSection(87);
      dScatter[88] = new BrowningEmpiricalCrossSection(88);
      dScatter[89] = new BrowningEmpiricalCrossSection(89);
      dScatter[90] = new BrowningEmpiricalCrossSection(90);
      dScatter[91] = new BrowningEmpiricalCrossSection(91);
      dScatter[92] = new BrowningEmpiricalCrossSection(92);
      dScatter[93] = new BrowningEmpiricalCrossSection(93);
      dScatter[94] = new BrowningEmpiricalCrossSection(94);
      dScatter[95] = new BrowningEmpiricalCrossSection(95);
      dScatter[96] = new BrowningEmpiricalCrossSection(96);
   }

   const BrowningEmpiricalCrossSection& getBECS(int an)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      return *dScatter[an];
#else
      return *mScatter[an];
#endif
   }
}
