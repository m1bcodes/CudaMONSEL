#include "gov\nist\microanalysis\EPQLibrary\BrowningEmpiricalCrossSection.cuh"

#include "Amphibian\random.cuh"

#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\FromSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

#include "gov\nist\microanalysis\Utility\Math2.cuh"

namespace BrowningEmpiricalCrossSection
{
   __host__ __device__ BrowningEmpiricalCrossSection::BrowningEmpiricalCrossSection(const ElementT& elm) :
      mElement(elm),
      mZp17(::powf(elm.getAtomicNumber(), 1.7)),
      mZp2(::powf(elm.getAtomicNumber(), 2.0)),
      mZp3(::powf(elm.getAtomicNumber(), 3.0))
   {
   }

   __host__ __device__ const ElementT& BrowningEmpiricalCrossSection::getElement() const
   {
      return mElement;
   }

   __host__ __device__ double BrowningEmpiricalCrossSection::totalCrossSection(const double energy) const
   {
      const double e = FromSI::keV(energy);
      const double re = ::sqrt(e);
      return 3.0e-22 * mZp17 / (e + 0.005 * mZp17 * re + 0.0007 * mZp2 / re);
   }

   __host__ __device__ double BrowningEmpiricalCrossSection::randomScatteringAngle(const double energy) const
   {
      const double r1 = Random::random(), r2 = Random::random();
      const double z = mElement.getAtomicNumber();
      const double e = FromSI::keV(energy);
      const double r = (300.0 * e / z) + (mZp3 / (3.0e5 * e));
      if (r1 <= r / (r + 1)) {
         // Screened Rutherford scattering
         const double alpha = 7.0e-3 / e;
         return ::acos(1.0 - ((2.0 * alpha * r2) / (alpha - r2 + 1)));
      }
      else {
         // Isotropic scattering
         return ::acos(1 - 2.0 * r2);
      }
   }

   //const BrowningEmpiricalCrossSection BECS1(1);
   //const BrowningEmpiricalCrossSection BECS2(2);
   //const BrowningEmpiricalCrossSection BECS3(3);
   //const BrowningEmpiricalCrossSection BECS4(4);
   //const BrowningEmpiricalCrossSection BECS5(5);
   //const BrowningEmpiricalCrossSection BECS6(6);
   //const BrowningEmpiricalCrossSection BECS7(7);
   //const BrowningEmpiricalCrossSection BECS8(8);
   //const BrowningEmpiricalCrossSection BECS9(9);
   //const BrowningEmpiricalCrossSection BECS10(10);
   //const BrowningEmpiricalCrossSection BECS11(11);
   //const BrowningEmpiricalCrossSection BECS12(12);
   //const BrowningEmpiricalCrossSection BECS13(13);
   //const BrowningEmpiricalCrossSection BECS14(14);
   //const BrowningEmpiricalCrossSection BECS15(15);
   //const BrowningEmpiricalCrossSection BECS16(16);
   //const BrowningEmpiricalCrossSection BECS17(17);
   //const BrowningEmpiricalCrossSection BECS18(18);
   //const BrowningEmpiricalCrossSection BECS19(19);
   //const BrowningEmpiricalCrossSection BECS20(20);
   //const BrowningEmpiricalCrossSection BECS21(21);
   //const BrowningEmpiricalCrossSection BECS22(22);
   //const BrowningEmpiricalCrossSection BECS23(23);
   //const BrowningEmpiricalCrossSection BECS24(24);
   //const BrowningEmpiricalCrossSection BECS25(25);
   //const BrowningEmpiricalCrossSection BECS26(26);
   //const BrowningEmpiricalCrossSection BECS27(27);
   //const BrowningEmpiricalCrossSection BECS28(28);
   //const BrowningEmpiricalCrossSection BECS29(29);
   //const BrowningEmpiricalCrossSection BECS30(30);
   //const BrowningEmpiricalCrossSection BECS31(31);
   //const BrowningEmpiricalCrossSection BECS32(32);
   //const BrowningEmpiricalCrossSection BECS33(33);
   //const BrowningEmpiricalCrossSection BECS34(34);
   //const BrowningEmpiricalCrossSection BECS35(35);
   //const BrowningEmpiricalCrossSection BECS36(36);
   //const BrowningEmpiricalCrossSection BECS37(37);
   //const BrowningEmpiricalCrossSection BECS38(38);
   //const BrowningEmpiricalCrossSection BECS39(39);
   //const BrowningEmpiricalCrossSection BECS40(40);
   //const BrowningEmpiricalCrossSection BECS41(41);
   //const BrowningEmpiricalCrossSection BECS42(42);
   //const BrowningEmpiricalCrossSection BECS43(43);
   //const BrowningEmpiricalCrossSection BECS44(44);
   //const BrowningEmpiricalCrossSection BECS45(45);
   //const BrowningEmpiricalCrossSection BECS46(46);
   //const BrowningEmpiricalCrossSection BECS47(47);
   //const BrowningEmpiricalCrossSection BECS48(48);
   //const BrowningEmpiricalCrossSection BECS49(49);
   //const BrowningEmpiricalCrossSection BECS50(50);
   //const BrowningEmpiricalCrossSection BECS51(51);
   //const BrowningEmpiricalCrossSection BECS52(52);
   //const BrowningEmpiricalCrossSection BECS53(53);
   //const BrowningEmpiricalCrossSection BECS54(54);
   //const BrowningEmpiricalCrossSection BECS55(55);
   //const BrowningEmpiricalCrossSection BECS56(56);
   //const BrowningEmpiricalCrossSection BECS57(57);
   //const BrowningEmpiricalCrossSection BECS58(58);
   //const BrowningEmpiricalCrossSection BECS59(59);
   //const BrowningEmpiricalCrossSection BECS60(60);
   //const BrowningEmpiricalCrossSection BECS61(61);
   //const BrowningEmpiricalCrossSection BECS62(62);
   //const BrowningEmpiricalCrossSection BECS63(63);
   //const BrowningEmpiricalCrossSection BECS64(64);
   //const BrowningEmpiricalCrossSection BECS65(65);
   //const BrowningEmpiricalCrossSection BECS66(66);
   //const BrowningEmpiricalCrossSection BECS67(67);
   //const BrowningEmpiricalCrossSection BECS68(68);
   //const BrowningEmpiricalCrossSection BECS69(69);
   //const BrowningEmpiricalCrossSection BECS70(70);
   //const BrowningEmpiricalCrossSection BECS71(71);
   //const BrowningEmpiricalCrossSection BECS72(72);
   //const BrowningEmpiricalCrossSection BECS73(73);
   //const BrowningEmpiricalCrossSection BECS74(74);
   //const BrowningEmpiricalCrossSection BECS75(75);
   //const BrowningEmpiricalCrossSection BECS76(76);
   //const BrowningEmpiricalCrossSection BECS77(77);
   //const BrowningEmpiricalCrossSection BECS78(78);
   //const BrowningEmpiricalCrossSection BECS79(79);
   //const BrowningEmpiricalCrossSection BECS80(80);
   //const BrowningEmpiricalCrossSection BECS81(81);
   //const BrowningEmpiricalCrossSection BECS82(82);
   //const BrowningEmpiricalCrossSection BECS83(83);
   //const BrowningEmpiricalCrossSection BECS84(84);
   //const BrowningEmpiricalCrossSection BECS85(85);
   //const BrowningEmpiricalCrossSection BECS86(86);
   //const BrowningEmpiricalCrossSection BECS87(87);
   //const BrowningEmpiricalCrossSection BECS88(88);
   //const BrowningEmpiricalCrossSection BECS89(89);
   //const BrowningEmpiricalCrossSection BECS90(90);
   //const BrowningEmpiricalCrossSection BECS91(91);
   //const BrowningEmpiricalCrossSection BECS92(92);
   //const BrowningEmpiricalCrossSection BECS93(93);
   //const BrowningEmpiricalCrossSection BECS94(94);
   //const BrowningEmpiricalCrossSection BECS95(95);
   //const BrowningEmpiricalCrossSection BECS96(96);

   //BrowningEmpiricalCrossSection const * mScatter[113] = {
   //   nullptr,
   //   &BECS1,
   //   &BECS2,
   //   &BECS3,
   //   &BECS4,
   //   &BECS5,
   //   &BECS6,
   //   &BECS7,
   //   &BECS8,
   //   &BECS9,
   //   &BECS10,
   //   &BECS11,
   //   &BECS12,
   //   &BECS13,
   //   &BECS14,
   //   &BECS15,
   //   &BECS16,
   //   &BECS17,
   //   &BECS18,
   //   &BECS19,
   //   &BECS20,
   //   &BECS21,
   //   &BECS22,
   //   &BECS23,
   //   &BECS24,
   //   &BECS25,
   //   &BECS26,
   //   &BECS27,
   //   &BECS28,
   //   &BECS29,
   //   &BECS30,
   //   &BECS31,
   //   &BECS32,
   //   &BECS33,
   //   &BECS34,
   //   &BECS35,
   //   &BECS36,
   //   &BECS37,
   //   &BECS38,
   //   &BECS39,
   //   &BECS40,
   //   &BECS41,
   //   &BECS42,
   //   &BECS43,
   //   &BECS44,
   //   &BECS45,
   //   &BECS46,
   //   &BECS47,
   //   &BECS48,
   //   &BECS49,
   //   &BECS50,
   //   &BECS51,
   //   &BECS52,
   //   &BECS53,
   //   &BECS54,
   //   &BECS55,
   //   &BECS56,
   //   &BECS57,
   //   &BECS58,
   //   &BECS59,
   //   &BECS60,
   //   &BECS61,
   //   &BECS62,
   //   &BECS63,
   //   &BECS64,
   //   &BECS65,
   //   &BECS66,
   //   &BECS67,
   //   &BECS68,
   //   &BECS69,
   //   &BECS70,
   //   &BECS71,
   //   &BECS72,
   //   &BECS73,
   //   &BECS74,
   //   &BECS75,
   //   &BECS76,
   //   &BECS77,
   //   &BECS78,
   //   &BECS79,
   //   &BECS80,
   //   &BECS81,
   //   &BECS82,
   //   &BECS83,
   //   &BECS84,
   //   &BECS85,
   //   &BECS86,
   //   &BECS87,
   //   &BECS88,
   //   &BECS89,
   //   &BECS90,
   //   &BECS91,
   //   &BECS92,
   //   &BECS93,
   //   &BECS94,
   //   &BECS95,
   //   &BECS96
   //};

   BrowningEmpiricalCrossSection const * mScatter[113];

   void init()
   {
      mScatter[1] = new BrowningEmpiricalCrossSection(Element::H);
      mScatter[2] = new BrowningEmpiricalCrossSection(Element::He);
      mScatter[3] = new BrowningEmpiricalCrossSection(Element::Li);
      mScatter[4] = new BrowningEmpiricalCrossSection(Element::Be);
      mScatter[5] = new BrowningEmpiricalCrossSection(Element::B);
      mScatter[6] = new BrowningEmpiricalCrossSection(Element::C);
      mScatter[7] = new BrowningEmpiricalCrossSection(Element::N);
      mScatter[8] = new BrowningEmpiricalCrossSection(Element::O);
      mScatter[9] = new BrowningEmpiricalCrossSection(Element::F);
      mScatter[10] = new BrowningEmpiricalCrossSection(Element::Ne);
      mScatter[11] = new BrowningEmpiricalCrossSection(Element::Na);
      mScatter[12] = new BrowningEmpiricalCrossSection(Element::Mg);
      mScatter[13] = new BrowningEmpiricalCrossSection(Element::Al);
      mScatter[14] = new BrowningEmpiricalCrossSection(Element::Si);
      mScatter[15] = new BrowningEmpiricalCrossSection(Element::P);
      mScatter[16] = new BrowningEmpiricalCrossSection(Element::S);
      mScatter[17] = new BrowningEmpiricalCrossSection(Element::Cl);
      mScatter[18] = new BrowningEmpiricalCrossSection(Element::Ar);
      mScatter[19] = new BrowningEmpiricalCrossSection(Element::K);
      mScatter[20] = new BrowningEmpiricalCrossSection(Element::Ca);
      mScatter[21] = new BrowningEmpiricalCrossSection(Element::Sc);
      mScatter[22] = new BrowningEmpiricalCrossSection(Element::Ti);
      mScatter[23] = new BrowningEmpiricalCrossSection(Element::V);
      mScatter[24] = new BrowningEmpiricalCrossSection(Element::Cr);
      mScatter[25] = new BrowningEmpiricalCrossSection(Element::Mn);
      mScatter[26] = new BrowningEmpiricalCrossSection(Element::Fe);
      mScatter[27] = new BrowningEmpiricalCrossSection(Element::Co);
      mScatter[28] = new BrowningEmpiricalCrossSection(Element::Ni);
      mScatter[29] = new BrowningEmpiricalCrossSection(Element::Cu);
      mScatter[30] = new BrowningEmpiricalCrossSection(Element::Zn);
      mScatter[31] = new BrowningEmpiricalCrossSection(Element::Ga);
      mScatter[32] = new BrowningEmpiricalCrossSection(Element::Ge);
      mScatter[33] = new BrowningEmpiricalCrossSection(Element::As);
      mScatter[34] = new BrowningEmpiricalCrossSection(Element::Se);
      mScatter[35] = new BrowningEmpiricalCrossSection(Element::Br);
      mScatter[36] = new BrowningEmpiricalCrossSection(Element::Kr);
      mScatter[37] = new BrowningEmpiricalCrossSection(Element::Rb);
      mScatter[38] = new BrowningEmpiricalCrossSection(Element::Sr);
      mScatter[39] = new BrowningEmpiricalCrossSection(Element::Y);
      mScatter[40] = new BrowningEmpiricalCrossSection(Element::Zr);
      mScatter[41] = new BrowningEmpiricalCrossSection(Element::Nb);
      mScatter[42] = new BrowningEmpiricalCrossSection(Element::Mo);
      mScatter[43] = new BrowningEmpiricalCrossSection(Element::Tc);
      mScatter[44] = new BrowningEmpiricalCrossSection(Element::Ru);
      mScatter[45] = new BrowningEmpiricalCrossSection(Element::Rh);
      mScatter[46] = new BrowningEmpiricalCrossSection(Element::Pd);
      mScatter[47] = new BrowningEmpiricalCrossSection(Element::Ag);
      mScatter[48] = new BrowningEmpiricalCrossSection(Element::Cd);
      mScatter[49] = new BrowningEmpiricalCrossSection(Element::In);
      mScatter[50] = new BrowningEmpiricalCrossSection(Element::Sn);
      mScatter[51] = new BrowningEmpiricalCrossSection(Element::Sb);
      mScatter[52] = new BrowningEmpiricalCrossSection(Element::Te);
      mScatter[53] = new BrowningEmpiricalCrossSection(Element::I);
      mScatter[54] = new BrowningEmpiricalCrossSection(Element::Xe);
      mScatter[55] = new BrowningEmpiricalCrossSection(Element::Cs);
      mScatter[56] = new BrowningEmpiricalCrossSection(Element::Ba);
      mScatter[57] = new BrowningEmpiricalCrossSection(Element::La);
      mScatter[58] = new BrowningEmpiricalCrossSection(Element::Ce);
      mScatter[59] = new BrowningEmpiricalCrossSection(Element::Pr);
      mScatter[60] = new BrowningEmpiricalCrossSection(Element::Nd);
      mScatter[61] = new BrowningEmpiricalCrossSection(Element::Pm);
      mScatter[62] = new BrowningEmpiricalCrossSection(Element::Sm);
      mScatter[63] = new BrowningEmpiricalCrossSection(Element::Eu);
      mScatter[64] = new BrowningEmpiricalCrossSection(Element::Gd);
      mScatter[65] = new BrowningEmpiricalCrossSection(Element::Tb);
      mScatter[66] = new BrowningEmpiricalCrossSection(Element::Dy);
      mScatter[67] = new BrowningEmpiricalCrossSection(Element::Ho);
      mScatter[68] = new BrowningEmpiricalCrossSection(Element::Er);
      mScatter[69] = new BrowningEmpiricalCrossSection(Element::Tm);
      mScatter[70] = new BrowningEmpiricalCrossSection(Element::Yb);
      mScatter[71] = new BrowningEmpiricalCrossSection(Element::Lu);
      mScatter[72] = new BrowningEmpiricalCrossSection(Element::Hf);
      mScatter[73] = new BrowningEmpiricalCrossSection(Element::Ta);
      mScatter[74] = new BrowningEmpiricalCrossSection(Element::W);
      mScatter[75] = new BrowningEmpiricalCrossSection(Element::Re);
      mScatter[76] = new BrowningEmpiricalCrossSection(Element::Os);
      mScatter[77] = new BrowningEmpiricalCrossSection(Element::Ir);
      mScatter[78] = new BrowningEmpiricalCrossSection(Element::Pt);
      mScatter[79] = new BrowningEmpiricalCrossSection(Element::Au);
      mScatter[80] = new BrowningEmpiricalCrossSection(Element::Hg);
      mScatter[81] = new BrowningEmpiricalCrossSection(Element::Tl);
      mScatter[82] = new BrowningEmpiricalCrossSection(Element::Pb);
      mScatter[83] = new BrowningEmpiricalCrossSection(Element::Bi);
      mScatter[84] = new BrowningEmpiricalCrossSection(Element::Po);
      mScatter[85] = new BrowningEmpiricalCrossSection(Element::At);
      mScatter[86] = new BrowningEmpiricalCrossSection(Element::Rn);
      mScatter[87] = new BrowningEmpiricalCrossSection(Element::Fr);
      mScatter[88] = new BrowningEmpiricalCrossSection(Element::Ra);
      mScatter[89] = new BrowningEmpiricalCrossSection(Element::Ac);
      mScatter[90] = new BrowningEmpiricalCrossSection(Element::Th);
      mScatter[91] = new BrowningEmpiricalCrossSection(Element::Pa);
      mScatter[92] = new BrowningEmpiricalCrossSection(Element::U);
      mScatter[93] = new BrowningEmpiricalCrossSection(Element::Np);
      mScatter[94] = new BrowningEmpiricalCrossSection(Element::Pu);
      mScatter[95] = new BrowningEmpiricalCrossSection(Element::Am);
      mScatter[96] = new BrowningEmpiricalCrossSection(Element::Cm);
   }

   __device__ BrowningEmpiricalCrossSection const * d_mScatter[113];

   __global__ void initCuda()
   {
      d_mScatter[1] = new BrowningEmpiricalCrossSection(*Element::dH);
      d_mScatter[2] = new BrowningEmpiricalCrossSection(*Element::dHe);
      d_mScatter[3] = new BrowningEmpiricalCrossSection(*Element::dLi);
      d_mScatter[4] = new BrowningEmpiricalCrossSection(*Element::dBe);
      d_mScatter[5] = new BrowningEmpiricalCrossSection(*Element::dB);
      d_mScatter[6] = new BrowningEmpiricalCrossSection(*Element::dC);
      d_mScatter[7] = new BrowningEmpiricalCrossSection(*Element::dN);
      d_mScatter[8] = new BrowningEmpiricalCrossSection(*Element::dO);
      d_mScatter[9] = new BrowningEmpiricalCrossSection(*Element::dF);
      d_mScatter[10] = new BrowningEmpiricalCrossSection(*Element::dNe);
      d_mScatter[11] = new BrowningEmpiricalCrossSection(*Element::dNa);
      d_mScatter[12] = new BrowningEmpiricalCrossSection(*Element::dMg);
      d_mScatter[13] = new BrowningEmpiricalCrossSection(*Element::dAl);
      d_mScatter[14] = new BrowningEmpiricalCrossSection(*Element::dSi);
      d_mScatter[15] = new BrowningEmpiricalCrossSection(*Element::dP);
      d_mScatter[16] = new BrowningEmpiricalCrossSection(*Element::dS);
      d_mScatter[17] = new BrowningEmpiricalCrossSection(*Element::dCl);
      d_mScatter[18] = new BrowningEmpiricalCrossSection(*Element::dAr);
      d_mScatter[19] = new BrowningEmpiricalCrossSection(*Element::dK);
      d_mScatter[20] = new BrowningEmpiricalCrossSection(*Element::dCa);
      d_mScatter[21] = new BrowningEmpiricalCrossSection(*Element::dSc);
      d_mScatter[22] = new BrowningEmpiricalCrossSection(*Element::dTi);
      d_mScatter[23] = new BrowningEmpiricalCrossSection(*Element::dV);
      d_mScatter[24] = new BrowningEmpiricalCrossSection(*Element::dCr);
      d_mScatter[25] = new BrowningEmpiricalCrossSection(*Element::dMn);
      d_mScatter[26] = new BrowningEmpiricalCrossSection(*Element::dFe);
      d_mScatter[27] = new BrowningEmpiricalCrossSection(*Element::dCo);
      d_mScatter[28] = new BrowningEmpiricalCrossSection(*Element::dNi);
      d_mScatter[29] = new BrowningEmpiricalCrossSection(*Element::dCu);
      d_mScatter[30] = new BrowningEmpiricalCrossSection(*Element::dZn);
      d_mScatter[31] = new BrowningEmpiricalCrossSection(*Element::dGa);
      d_mScatter[32] = new BrowningEmpiricalCrossSection(*Element::dGe);
      d_mScatter[33] = new BrowningEmpiricalCrossSection(*Element::dAs);
      d_mScatter[34] = new BrowningEmpiricalCrossSection(*Element::dSe);
      d_mScatter[35] = new BrowningEmpiricalCrossSection(*Element::dBr);
      d_mScatter[36] = new BrowningEmpiricalCrossSection(*Element::dKr);
      d_mScatter[37] = new BrowningEmpiricalCrossSection(*Element::dRb);
      d_mScatter[38] = new BrowningEmpiricalCrossSection(*Element::dSr);
      d_mScatter[39] = new BrowningEmpiricalCrossSection(*Element::dY);
      d_mScatter[40] = new BrowningEmpiricalCrossSection(*Element::dZr);
      d_mScatter[41] = new BrowningEmpiricalCrossSection(*Element::dNb);
      d_mScatter[42] = new BrowningEmpiricalCrossSection(*Element::dMo);
      d_mScatter[43] = new BrowningEmpiricalCrossSection(*Element::dTc);
      d_mScatter[44] = new BrowningEmpiricalCrossSection(*Element::dRu);
      d_mScatter[45] = new BrowningEmpiricalCrossSection(*Element::dRh);
      d_mScatter[46] = new BrowningEmpiricalCrossSection(*Element::dPd);
      d_mScatter[47] = new BrowningEmpiricalCrossSection(*Element::dAg);
      d_mScatter[48] = new BrowningEmpiricalCrossSection(*Element::dCd);
      d_mScatter[49] = new BrowningEmpiricalCrossSection(*Element::dIn);
      d_mScatter[50] = new BrowningEmpiricalCrossSection(*Element::dSn);
      d_mScatter[51] = new BrowningEmpiricalCrossSection(*Element::dSb);
      d_mScatter[52] = new BrowningEmpiricalCrossSection(*Element::dTe);
      d_mScatter[53] = new BrowningEmpiricalCrossSection(*Element::dI);
      d_mScatter[54] = new BrowningEmpiricalCrossSection(*Element::dXe);
      d_mScatter[55] = new BrowningEmpiricalCrossSection(*Element::dCs);
      d_mScatter[56] = new BrowningEmpiricalCrossSection(*Element::dBa);
      d_mScatter[57] = new BrowningEmpiricalCrossSection(*Element::dLa);
      d_mScatter[58] = new BrowningEmpiricalCrossSection(*Element::dCe);
      d_mScatter[59] = new BrowningEmpiricalCrossSection(*Element::dPr);
      d_mScatter[60] = new BrowningEmpiricalCrossSection(*Element::dNd);
      d_mScatter[61] = new BrowningEmpiricalCrossSection(*Element::dPm);
      d_mScatter[62] = new BrowningEmpiricalCrossSection(*Element::dSm);
      d_mScatter[63] = new BrowningEmpiricalCrossSection(*Element::dEu);
      d_mScatter[64] = new BrowningEmpiricalCrossSection(*Element::dGd);
      d_mScatter[65] = new BrowningEmpiricalCrossSection(*Element::dTb);
      d_mScatter[66] = new BrowningEmpiricalCrossSection(*Element::dDy);
      d_mScatter[67] = new BrowningEmpiricalCrossSection(*Element::dHo);
      d_mScatter[68] = new BrowningEmpiricalCrossSection(*Element::dEr);
      d_mScatter[69] = new BrowningEmpiricalCrossSection(*Element::dTm);
      d_mScatter[70] = new BrowningEmpiricalCrossSection(*Element::dYb);
      d_mScatter[71] = new BrowningEmpiricalCrossSection(*Element::dLu);
      d_mScatter[72] = new BrowningEmpiricalCrossSection(*Element::dHf);
      d_mScatter[73] = new BrowningEmpiricalCrossSection(*Element::dTa);
      d_mScatter[74] = new BrowningEmpiricalCrossSection(*Element::dW);
      d_mScatter[75] = new BrowningEmpiricalCrossSection(*Element::dRe);
      d_mScatter[76] = new BrowningEmpiricalCrossSection(*Element::dOs);
      d_mScatter[77] = new BrowningEmpiricalCrossSection(*Element::dIr);
      d_mScatter[78] = new BrowningEmpiricalCrossSection(*Element::dPt);
      d_mScatter[79] = new BrowningEmpiricalCrossSection(*Element::dAu);
      d_mScatter[80] = new BrowningEmpiricalCrossSection(*Element::dHg);
      d_mScatter[81] = new BrowningEmpiricalCrossSection(*Element::dTl);
      d_mScatter[82] = new BrowningEmpiricalCrossSection(*Element::dPb);
      d_mScatter[83] = new BrowningEmpiricalCrossSection(*Element::dBi);
      d_mScatter[84] = new BrowningEmpiricalCrossSection(*Element::dPo);
      d_mScatter[85] = new BrowningEmpiricalCrossSection(*Element::dAt);
      d_mScatter[86] = new BrowningEmpiricalCrossSection(*Element::dRn);
      d_mScatter[87] = new BrowningEmpiricalCrossSection(*Element::dFr);
      d_mScatter[88] = new BrowningEmpiricalCrossSection(*Element::dRa);
      d_mScatter[89] = new BrowningEmpiricalCrossSection(*Element::dAc);
      d_mScatter[90] = new BrowningEmpiricalCrossSection(*Element::dTh);
      d_mScatter[91] = new BrowningEmpiricalCrossSection(*Element::dPa);
      d_mScatter[92] = new BrowningEmpiricalCrossSection(*Element::dU);
      d_mScatter[93] = new BrowningEmpiricalCrossSection(*Element::dNp);
      d_mScatter[94] = new BrowningEmpiricalCrossSection(*Element::dPu);
      d_mScatter[95] = new BrowningEmpiricalCrossSection(*Element::dAm);
      d_mScatter[96] = new BrowningEmpiricalCrossSection(*Element::dCm);
   }

   __host__ __device__ const BrowningEmpiricalCrossSection& getBECS(int an)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      return *d_mScatter[an];
#else
      return *mScatter[an];
#endif
   }
}
