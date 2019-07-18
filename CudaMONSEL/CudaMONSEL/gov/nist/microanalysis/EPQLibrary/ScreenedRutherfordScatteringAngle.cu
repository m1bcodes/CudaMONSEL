#include "gov\nist\microanalysis\EPQLibrary\ScreenedRutherfordScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

#include "gov\nist\microanalysis\Utility\Math2.cuh"

#include "Amphibian\random.cuh"

namespace ScreenedRutherfordScatteringAngle
{
   static const Reference::CrudeReference REFERENCE("NBSMONTE.FOR");

   __host__ __device__ ScreenedRutherfordScatteringAngle::ScreenedRutherfordScatteringAngle(const ElementT& elm) :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      RandomizedScatterT("Screened Rutherford", *Reference::d_NullReference),
#else
      RandomizedScatterT("Screened Rutherford", REFERENCE),
#endif
      mElement(elm)
   {
   }

   StringT ScreenedRutherfordScatteringAngle::toString() const
   {
      return "CrossSection[Screened-Rutherford," + StringT(mElement.toAbbrev()) + "]";
   }

   __host__ __device__ const ElementT& ScreenedRutherfordScatteringAngle::getElement() const
   {
      return mElement;
   }

   __host__ __device__ double ScreenedRutherfordScatteringAngle::totalCrossSection(const double energy) const
   {
      // Ref: Heinrich 1981 p 459 convert to SI units
      const double z = mElement.getAtomicNumber();
      const double zp = ::powf(z, 1.0 / 3.0);
      return (7.670843088080456e-38 * zp * (1.0 + z)) / ((energy + ((5.44967975966321e-19 * zp * zp))));
   }

   __host__ __device__ double ScreenedRutherfordScatteringAngle::randomScatteringAngle(const double energy) const
   {
      // This method for calculating the scattering angle is taken from
      // NBSMONTE.FOR
      const double alpha = (5.44968e-19 * ::powf(mElement.getAtomicNumber(), 2.0 / 3.0)) / energy;
      const double r = Random::random();
      return ::acos(1 - 2.0 * alpha * r / (1 + alpha - r));
   }

   //const ScreenedRutherfordScatteringAngle SRSA1(Element::H);
   //const ScreenedRutherfordScatteringAngle SRSA2(Element::He);
   //const ScreenedRutherfordScatteringAngle SRSA3(Element::Li);
   //const ScreenedRutherfordScatteringAngle SRSA4(Element::Be);
   //const ScreenedRutherfordScatteringAngle SRSA5(Element::B);
   //const ScreenedRutherfordScatteringAngle SRSA6(Element::C);
   //const ScreenedRutherfordScatteringAngle SRSA7(Element::N);
   //const ScreenedRutherfordScatteringAngle SRSA8(Element::O);
   //const ScreenedRutherfordScatteringAngle SRSA9(Element::F);
   //const ScreenedRutherfordScatteringAngle SRSA10(Element::Ne);
   //const ScreenedRutherfordScatteringAngle SRSA11(Element::Na);
   //const ScreenedRutherfordScatteringAngle SRSA12(Element::Mg);
   //const ScreenedRutherfordScatteringAngle SRSA13(Element::Al);
   //const ScreenedRutherfordScatteringAngle SRSA14(Element::Si);
   //const ScreenedRutherfordScatteringAngle SRSA15(Element::P);
   //const ScreenedRutherfordScatteringAngle SRSA16(Element::S);
   //const ScreenedRutherfordScatteringAngle SRSA17(Element::Cl);
   //const ScreenedRutherfordScatteringAngle SRSA18(Element::Ar);
   //const ScreenedRutherfordScatteringAngle SRSA19(Element::K);
   //const ScreenedRutherfordScatteringAngle SRSA20(Element::Ca);
   //const ScreenedRutherfordScatteringAngle SRSA21(Element::Sc);
   //const ScreenedRutherfordScatteringAngle SRSA22(Element::Ti);
   //const ScreenedRutherfordScatteringAngle SRSA23(Element::V);
   //const ScreenedRutherfordScatteringAngle SRSA24(Element::Cr);
   //const ScreenedRutherfordScatteringAngle SRSA25(Element::Mn);
   //const ScreenedRutherfordScatteringAngle SRSA26(Element::Fe);
   //const ScreenedRutherfordScatteringAngle SRSA27(Element::Co);
   //const ScreenedRutherfordScatteringAngle SRSA28(Element::Ni);
   //const ScreenedRutherfordScatteringAngle SRSA29(Element::Cu);
   //const ScreenedRutherfordScatteringAngle SRSA30(Element::Zn);
   //const ScreenedRutherfordScatteringAngle SRSA31(Element::Ga);
   //const ScreenedRutherfordScatteringAngle SRSA32(Element::Ge);
   //const ScreenedRutherfordScatteringAngle SRSA33(Element::As);
   //const ScreenedRutherfordScatteringAngle SRSA34(Element::Se);
   //const ScreenedRutherfordScatteringAngle SRSA35(Element::Br);
   //const ScreenedRutherfordScatteringAngle SRSA36(Element::Kr);
   //const ScreenedRutherfordScatteringAngle SRSA37(Element::Rb);
   //const ScreenedRutherfordScatteringAngle SRSA38(Element::Sr);
   //const ScreenedRutherfordScatteringAngle SRSA39(Element::Y);
   //const ScreenedRutherfordScatteringAngle SRSA40(Element::Zr);
   //const ScreenedRutherfordScatteringAngle SRSA41(Element::Nb);
   //const ScreenedRutherfordScatteringAngle SRSA42(Element::Mo);
   //const ScreenedRutherfordScatteringAngle SRSA43(Element::Tc);
   //const ScreenedRutherfordScatteringAngle SRSA44(Element::Ru);
   //const ScreenedRutherfordScatteringAngle SRSA45(Element::Rh);
   //const ScreenedRutherfordScatteringAngle SRSA46(Element::Pd);
   //const ScreenedRutherfordScatteringAngle SRSA47(Element::Ag);
   //const ScreenedRutherfordScatteringAngle SRSA48(Element::Cd);
   //const ScreenedRutherfordScatteringAngle SRSA49(Element::In);
   //const ScreenedRutherfordScatteringAngle SRSA50(Element::Sn);
   //const ScreenedRutherfordScatteringAngle SRSA51(Element::Sb);
   //const ScreenedRutherfordScatteringAngle SRSA52(Element::Te);
   //const ScreenedRutherfordScatteringAngle SRSA53(Element::I);
   //const ScreenedRutherfordScatteringAngle SRSA54(Element::Xe);
   //const ScreenedRutherfordScatteringAngle SRSA55(Element::Cs);
   //const ScreenedRutherfordScatteringAngle SRSA56(Element::Ba);
   //const ScreenedRutherfordScatteringAngle SRSA57(Element::La);
   //const ScreenedRutherfordScatteringAngle SRSA58(Element::Ce);
   //const ScreenedRutherfordScatteringAngle SRSA59(Element::Pr);
   //const ScreenedRutherfordScatteringAngle SRSA60(Element::Nd);
   //const ScreenedRutherfordScatteringAngle SRSA61(Element::Pm);
   //const ScreenedRutherfordScatteringAngle SRSA62(Element::Sm);
   //const ScreenedRutherfordScatteringAngle SRSA63(Element::Eu);
   //const ScreenedRutherfordScatteringAngle SRSA64(Element::Gd);
   //const ScreenedRutherfordScatteringAngle SRSA65(Element::Tb);
   //const ScreenedRutherfordScatteringAngle SRSA66(Element::Dy);
   //const ScreenedRutherfordScatteringAngle SRSA67(Element::Ho);
   //const ScreenedRutherfordScatteringAngle SRSA68(Element::Er);
   //const ScreenedRutherfordScatteringAngle SRSA69(Element::Tm);
   //const ScreenedRutherfordScatteringAngle SRSA70(Element::Yb);
   //const ScreenedRutherfordScatteringAngle SRSA71(Element::Lu);
   //const ScreenedRutherfordScatteringAngle SRSA72(Element::Hf);
   //const ScreenedRutherfordScatteringAngle SRSA73(Element::Ta);
   //const ScreenedRutherfordScatteringAngle SRSA74(Element::W);
   //const ScreenedRutherfordScatteringAngle SRSA75(Element::Re);
   //const ScreenedRutherfordScatteringAngle SRSA76(Element::Os);
   //const ScreenedRutherfordScatteringAngle SRSA77(Element::Ir);
   //const ScreenedRutherfordScatteringAngle SRSA78(Element::Pt);
   //const ScreenedRutherfordScatteringAngle SRSA79(Element::Au);
   //const ScreenedRutherfordScatteringAngle SRSA80(Element::Hg);
   //const ScreenedRutherfordScatteringAngle SRSA81(Element::Tl);
   //const ScreenedRutherfordScatteringAngle SRSA82(Element::Pb);
   //const ScreenedRutherfordScatteringAngle SRSA83(Element::Bi);
   //const ScreenedRutherfordScatteringAngle SRSA84(Element::Po);
   //const ScreenedRutherfordScatteringAngle SRSA85(Element::At);
   //const ScreenedRutherfordScatteringAngle SRSA86(Element::Rn);
   //const ScreenedRutherfordScatteringAngle SRSA87(Element::Fr);
   //const ScreenedRutherfordScatteringAngle SRSA88(Element::Ra);
   //const ScreenedRutherfordScatteringAngle SRSA89(Element::Ac);
   //const ScreenedRutherfordScatteringAngle SRSA90(Element::Th);
   //const ScreenedRutherfordScatteringAngle SRSA91(Element::Pa);
   //const ScreenedRutherfordScatteringAngle SRSA92(Element::U);
   //const ScreenedRutherfordScatteringAngle SRSA93(Element::Np);
   //const ScreenedRutherfordScatteringAngle SRSA94(Element::Pu);
   //const ScreenedRutherfordScatteringAngle SRSA95(Element::Am);
   //const ScreenedRutherfordScatteringAngle SRSA96(Element::Cm);

   //ScreenedRutherfordScatteringAngle const * mScatter[113] = {
   //   nullptr,
   //   &SRSA1,
   //   &SRSA2,
   //   &SRSA3,
   //   &SRSA4,
   //   &SRSA5,
   //   &SRSA6,
   //   &SRSA7,
   //   &SRSA8,
   //   &SRSA9,
   //   &SRSA10,
   //   &SRSA11,
   //   &SRSA12,
   //   &SRSA13,
   //   &SRSA14,
   //   &SRSA15,
   //   &SRSA16,
   //   &SRSA17,
   //   &SRSA18,
   //   &SRSA19,
   //   &SRSA20,
   //   &SRSA21,
   //   &SRSA22,
   //   &SRSA23,
   //   &SRSA24,
   //   &SRSA25,
   //   &SRSA26,
   //   &SRSA27,
   //   &SRSA28,
   //   &SRSA29,
   //   &SRSA30,
   //   &SRSA31,
   //   &SRSA32,
   //   &SRSA33,
   //   &SRSA34,
   //   &SRSA35,
   //   &SRSA36,
   //   &SRSA37,
   //   &SRSA38,
   //   &SRSA39,
   //   &SRSA40,
   //   &SRSA41,
   //   &SRSA42,
   //   &SRSA43,
   //   &SRSA44,
   //   &SRSA45,
   //   &SRSA46,
   //   &SRSA47,
   //   &SRSA48,
   //   &SRSA49,
   //   &SRSA50,
   //   &SRSA51,
   //   &SRSA52,
   //   &SRSA53,
   //   &SRSA54,
   //   &SRSA55,
   //   &SRSA56,
   //   &SRSA57,
   //   &SRSA58,
   //   &SRSA59,
   //   &SRSA60,
   //   &SRSA61,
   //   &SRSA62,
   //   &SRSA63,
   //   &SRSA64,
   //   &SRSA65,
   //   &SRSA66,
   //   &SRSA67,
   //   &SRSA68,
   //   &SRSA69,
   //   &SRSA70,
   //   &SRSA71,
   //   &SRSA72,
   //   &SRSA73,
   //   &SRSA74,
   //   &SRSA75,
   //   &SRSA76,
   //   &SRSA77,
   //   &SRSA78,
   //   &SRSA79,
   //   &SRSA80,
   //   &SRSA81,
   //   &SRSA82,
   //   &SRSA83,
   //   &SRSA84,
   //   &SRSA85,
   //   &SRSA86,
   //   &SRSA87,
   //   &SRSA88,
   //   &SRSA89,
   //   &SRSA90,
   //   &SRSA91,
   //   &SRSA92,
   //   &SRSA93,
   //   &SRSA94,
   //   &SRSA95,
   //   &SRSA96
   //};

   ScreenedRutherfordScatteringAngle const * mScatter[113];

   void init()
   {
      mScatter[1] = new ScreenedRutherfordScatteringAngle(Element::H);
      mScatter[2] = new ScreenedRutherfordScatteringAngle(Element::He);
      mScatter[3] = new ScreenedRutherfordScatteringAngle(Element::Li);
      mScatter[4] = new ScreenedRutherfordScatteringAngle(Element::Be);
      mScatter[5] = new ScreenedRutherfordScatteringAngle(Element::B);
      mScatter[6] = new ScreenedRutherfordScatteringAngle(Element::C);
      mScatter[7] = new ScreenedRutherfordScatteringAngle(Element::N);
      mScatter[8] = new ScreenedRutherfordScatteringAngle(Element::O);
      mScatter[9] = new ScreenedRutherfordScatteringAngle(Element::F);
      mScatter[10] = new ScreenedRutherfordScatteringAngle(Element::Ne);
      mScatter[11] = new ScreenedRutherfordScatteringAngle(Element::Na);
      mScatter[12] = new ScreenedRutherfordScatteringAngle(Element::Mg);
      mScatter[13] = new ScreenedRutherfordScatteringAngle(Element::Al);
      mScatter[14] = new ScreenedRutherfordScatteringAngle(Element::Si);
      mScatter[15] = new ScreenedRutherfordScatteringAngle(Element::P);
      mScatter[16] = new ScreenedRutherfordScatteringAngle(Element::S);
      mScatter[17] = new ScreenedRutherfordScatteringAngle(Element::Cl);
      mScatter[18] = new ScreenedRutherfordScatteringAngle(Element::Ar);
      mScatter[19] = new ScreenedRutherfordScatteringAngle(Element::K);
      mScatter[20] = new ScreenedRutherfordScatteringAngle(Element::Ca);
      mScatter[21] = new ScreenedRutherfordScatteringAngle(Element::Sc);
      mScatter[22] = new ScreenedRutherfordScatteringAngle(Element::Ti);
      mScatter[23] = new ScreenedRutherfordScatteringAngle(Element::V);
      mScatter[24] = new ScreenedRutherfordScatteringAngle(Element::Cr);
      mScatter[25] = new ScreenedRutherfordScatteringAngle(Element::Mn);
      mScatter[26] = new ScreenedRutherfordScatteringAngle(Element::Fe);
      mScatter[27] = new ScreenedRutherfordScatteringAngle(Element::Co);
      mScatter[28] = new ScreenedRutherfordScatteringAngle(Element::Ni);
      mScatter[29] = new ScreenedRutherfordScatteringAngle(Element::Cu);
      mScatter[30] = new ScreenedRutherfordScatteringAngle(Element::Zn);
      mScatter[31] = new ScreenedRutherfordScatteringAngle(Element::Ga);
      mScatter[32] = new ScreenedRutherfordScatteringAngle(Element::Ge);
      mScatter[33] = new ScreenedRutherfordScatteringAngle(Element::As);
      mScatter[34] = new ScreenedRutherfordScatteringAngle(Element::Se);
      mScatter[35] = new ScreenedRutherfordScatteringAngle(Element::Br);
      mScatter[36] = new ScreenedRutherfordScatteringAngle(Element::Kr);
      mScatter[37] = new ScreenedRutherfordScatteringAngle(Element::Rb);
      mScatter[38] = new ScreenedRutherfordScatteringAngle(Element::Sr);
      mScatter[39] = new ScreenedRutherfordScatteringAngle(Element::Y);
      mScatter[40] = new ScreenedRutherfordScatteringAngle(Element::Zr);
      mScatter[41] = new ScreenedRutherfordScatteringAngle(Element::Nb);
      mScatter[42] = new ScreenedRutherfordScatteringAngle(Element::Mo);
      mScatter[43] = new ScreenedRutherfordScatteringAngle(Element::Tc);
      mScatter[44] = new ScreenedRutherfordScatteringAngle(Element::Ru);
      mScatter[45] = new ScreenedRutherfordScatteringAngle(Element::Rh);
      mScatter[46] = new ScreenedRutherfordScatteringAngle(Element::Pd);
      mScatter[47] = new ScreenedRutherfordScatteringAngle(Element::Ag);
      mScatter[48] = new ScreenedRutherfordScatteringAngle(Element::Cd);
      mScatter[49] = new ScreenedRutherfordScatteringAngle(Element::In);
      mScatter[50] = new ScreenedRutherfordScatteringAngle(Element::Sn);
      mScatter[51] = new ScreenedRutherfordScatteringAngle(Element::Sb);
      mScatter[52] = new ScreenedRutherfordScatteringAngle(Element::Te);
      mScatter[53] = new ScreenedRutherfordScatteringAngle(Element::I);
      mScatter[54] = new ScreenedRutherfordScatteringAngle(Element::Xe);
      mScatter[55] = new ScreenedRutherfordScatteringAngle(Element::Cs);
      mScatter[56] = new ScreenedRutherfordScatteringAngle(Element::Ba);
      mScatter[57] = new ScreenedRutherfordScatteringAngle(Element::La);
      mScatter[58] = new ScreenedRutherfordScatteringAngle(Element::Ce);
      mScatter[59] = new ScreenedRutherfordScatteringAngle(Element::Pr);
      mScatter[60] = new ScreenedRutherfordScatteringAngle(Element::Nd);
      mScatter[61] = new ScreenedRutherfordScatteringAngle(Element::Pm);
      mScatter[62] = new ScreenedRutherfordScatteringAngle(Element::Sm);
      mScatter[63] = new ScreenedRutherfordScatteringAngle(Element::Eu);
      mScatter[64] = new ScreenedRutherfordScatteringAngle(Element::Gd);
      mScatter[65] = new ScreenedRutherfordScatteringAngle(Element::Tb);
      mScatter[66] = new ScreenedRutherfordScatteringAngle(Element::Dy);
      mScatter[67] = new ScreenedRutherfordScatteringAngle(Element::Ho);
      mScatter[68] = new ScreenedRutherfordScatteringAngle(Element::Er);
      mScatter[69] = new ScreenedRutherfordScatteringAngle(Element::Tm);
      mScatter[70] = new ScreenedRutherfordScatteringAngle(Element::Yb);
      mScatter[71] = new ScreenedRutherfordScatteringAngle(Element::Lu);
      mScatter[72] = new ScreenedRutherfordScatteringAngle(Element::Hf);
      mScatter[73] = new ScreenedRutherfordScatteringAngle(Element::Ta);
      mScatter[74] = new ScreenedRutherfordScatteringAngle(Element::W);
      mScatter[75] = new ScreenedRutherfordScatteringAngle(Element::Re);
      mScatter[76] = new ScreenedRutherfordScatteringAngle(Element::Os);
      mScatter[77] = new ScreenedRutherfordScatteringAngle(Element::Ir);
      mScatter[78] = new ScreenedRutherfordScatteringAngle(Element::Pt);
      mScatter[79] = new ScreenedRutherfordScatteringAngle(Element::Au);
      mScatter[80] = new ScreenedRutherfordScatteringAngle(Element::Hg);
      mScatter[81] = new ScreenedRutherfordScatteringAngle(Element::Tl);
      mScatter[82] = new ScreenedRutherfordScatteringAngle(Element::Pb);
      mScatter[83] = new ScreenedRutherfordScatteringAngle(Element::Bi);
      mScatter[84] = new ScreenedRutherfordScatteringAngle(Element::Po);
      mScatter[85] = new ScreenedRutherfordScatteringAngle(Element::At);
      mScatter[86] = new ScreenedRutherfordScatteringAngle(Element::Rn);
      mScatter[87] = new ScreenedRutherfordScatteringAngle(Element::Fr);
      mScatter[88] = new ScreenedRutherfordScatteringAngle(Element::Ra);
      mScatter[89] = new ScreenedRutherfordScatteringAngle(Element::Ac);
      mScatter[90] = new ScreenedRutherfordScatteringAngle(Element::Th);
      mScatter[91] = new ScreenedRutherfordScatteringAngle(Element::Pa);
      mScatter[92] = new ScreenedRutherfordScatteringAngle(Element::U);
      mScatter[93] = new ScreenedRutherfordScatteringAngle(Element::Np);
      mScatter[94] = new ScreenedRutherfordScatteringAngle(Element::Pu);
      mScatter[95] = new ScreenedRutherfordScatteringAngle(Element::Am);
      mScatter[96] = new ScreenedRutherfordScatteringAngle(Element::Cm);
   }

   __device__ ScreenedRutherfordScatteringAngle const * dScatter[113];

   __global__ void initCuda()
   {
      dScatter[1] = new ScreenedRutherfordScatteringAngle(*Element::dH);
      dScatter[2] = new ScreenedRutherfordScatteringAngle(*Element::dHe);
      dScatter[3] = new ScreenedRutherfordScatteringAngle(*Element::dLi);
      dScatter[4] = new ScreenedRutherfordScatteringAngle(*Element::dBe);
      dScatter[5] = new ScreenedRutherfordScatteringAngle(*Element::dB);
      dScatter[6] = new ScreenedRutherfordScatteringAngle(*Element::dC);
      dScatter[7] = new ScreenedRutherfordScatteringAngle(*Element::dN);
      dScatter[8] = new ScreenedRutherfordScatteringAngle(*Element::dO);
      dScatter[9] = new ScreenedRutherfordScatteringAngle(*Element::dF);
      dScatter[10] = new ScreenedRutherfordScatteringAngle(*Element::dNe);
      dScatter[11] = new ScreenedRutherfordScatteringAngle(*Element::dNa);
      dScatter[12] = new ScreenedRutherfordScatteringAngle(*Element::dMg);
      dScatter[13] = new ScreenedRutherfordScatteringAngle(*Element::dAl);
      dScatter[14] = new ScreenedRutherfordScatteringAngle(*Element::dSi);
      dScatter[15] = new ScreenedRutherfordScatteringAngle(*Element::dP);
      dScatter[16] = new ScreenedRutherfordScatteringAngle(*Element::dS);
      dScatter[17] = new ScreenedRutherfordScatteringAngle(*Element::dCl);
      dScatter[18] = new ScreenedRutherfordScatteringAngle(*Element::dAr);
      dScatter[19] = new ScreenedRutherfordScatteringAngle(*Element::dK);
      dScatter[20] = new ScreenedRutherfordScatteringAngle(*Element::dCa);
      dScatter[21] = new ScreenedRutherfordScatteringAngle(*Element::dSc);
      dScatter[22] = new ScreenedRutherfordScatteringAngle(*Element::dTi);
      dScatter[23] = new ScreenedRutherfordScatteringAngle(*Element::dV);
      dScatter[24] = new ScreenedRutherfordScatteringAngle(*Element::dCr);
      dScatter[25] = new ScreenedRutherfordScatteringAngle(*Element::dMn);
      dScatter[26] = new ScreenedRutherfordScatteringAngle(*Element::dFe);
      dScatter[27] = new ScreenedRutherfordScatteringAngle(*Element::dCo);
      dScatter[28] = new ScreenedRutherfordScatteringAngle(*Element::dNi);
      dScatter[29] = new ScreenedRutherfordScatteringAngle(*Element::dCu);
      dScatter[30] = new ScreenedRutherfordScatteringAngle(*Element::dZn);
      dScatter[31] = new ScreenedRutherfordScatteringAngle(*Element::dGa);
      dScatter[32] = new ScreenedRutherfordScatteringAngle(*Element::dGe);
      dScatter[33] = new ScreenedRutherfordScatteringAngle(*Element::dAs);
      dScatter[34] = new ScreenedRutherfordScatteringAngle(*Element::dSe);
      dScatter[35] = new ScreenedRutherfordScatteringAngle(*Element::dBr);
      dScatter[36] = new ScreenedRutherfordScatteringAngle(*Element::dKr);
      dScatter[37] = new ScreenedRutherfordScatteringAngle(*Element::dRb);
      dScatter[38] = new ScreenedRutherfordScatteringAngle(*Element::dSr);
      dScatter[39] = new ScreenedRutherfordScatteringAngle(*Element::dY);
      dScatter[40] = new ScreenedRutherfordScatteringAngle(*Element::dZr);
      dScatter[41] = new ScreenedRutherfordScatteringAngle(*Element::dNb);
      dScatter[42] = new ScreenedRutherfordScatteringAngle(*Element::dMo);
      dScatter[43] = new ScreenedRutherfordScatteringAngle(*Element::dTc);
      dScatter[44] = new ScreenedRutherfordScatteringAngle(*Element::dRu);
      dScatter[45] = new ScreenedRutherfordScatteringAngle(*Element::dRh);
      dScatter[46] = new ScreenedRutherfordScatteringAngle(*Element::dPd);
      dScatter[47] = new ScreenedRutherfordScatteringAngle(*Element::dAg);
      dScatter[48] = new ScreenedRutherfordScatteringAngle(*Element::dCd);
      dScatter[49] = new ScreenedRutherfordScatteringAngle(*Element::dIn);
      dScatter[50] = new ScreenedRutherfordScatteringAngle(*Element::dSn);
      dScatter[51] = new ScreenedRutherfordScatteringAngle(*Element::dSb);
      dScatter[52] = new ScreenedRutherfordScatteringAngle(*Element::dTe);
      dScatter[53] = new ScreenedRutherfordScatteringAngle(*Element::dI);
      dScatter[54] = new ScreenedRutherfordScatteringAngle(*Element::dXe);
      dScatter[55] = new ScreenedRutherfordScatteringAngle(*Element::dCs);
      dScatter[56] = new ScreenedRutherfordScatteringAngle(*Element::dBa);
      dScatter[57] = new ScreenedRutherfordScatteringAngle(*Element::dLa);
      dScatter[58] = new ScreenedRutherfordScatteringAngle(*Element::dCe);
      dScatter[59] = new ScreenedRutherfordScatteringAngle(*Element::dPr);
      dScatter[60] = new ScreenedRutherfordScatteringAngle(*Element::dNd);
      dScatter[61] = new ScreenedRutherfordScatteringAngle(*Element::dPm);
      dScatter[62] = new ScreenedRutherfordScatteringAngle(*Element::dSm);
      dScatter[63] = new ScreenedRutherfordScatteringAngle(*Element::dEu);
      dScatter[64] = new ScreenedRutherfordScatteringAngle(*Element::dGd);
      dScatter[65] = new ScreenedRutherfordScatteringAngle(*Element::dTb);
      dScatter[66] = new ScreenedRutherfordScatteringAngle(*Element::dDy);
      dScatter[67] = new ScreenedRutherfordScatteringAngle(*Element::dHo);
      dScatter[68] = new ScreenedRutherfordScatteringAngle(*Element::dEr);
      dScatter[69] = new ScreenedRutherfordScatteringAngle(*Element::dTm);
      dScatter[70] = new ScreenedRutherfordScatteringAngle(*Element::dYb);
      dScatter[71] = new ScreenedRutherfordScatteringAngle(*Element::dLu);
      dScatter[72] = new ScreenedRutherfordScatteringAngle(*Element::dHf);
      dScatter[73] = new ScreenedRutherfordScatteringAngle(*Element::dTa);
      dScatter[74] = new ScreenedRutherfordScatteringAngle(*Element::dW);
      dScatter[75] = new ScreenedRutherfordScatteringAngle(*Element::dRe);
      dScatter[76] = new ScreenedRutherfordScatteringAngle(*Element::dOs);
      dScatter[77] = new ScreenedRutherfordScatteringAngle(*Element::dIr);
      dScatter[78] = new ScreenedRutherfordScatteringAngle(*Element::dPt);
      dScatter[79] = new ScreenedRutherfordScatteringAngle(*Element::dAu);
      dScatter[80] = new ScreenedRutherfordScatteringAngle(*Element::dHg);
      dScatter[81] = new ScreenedRutherfordScatteringAngle(*Element::dTl);
      dScatter[82] = new ScreenedRutherfordScatteringAngle(*Element::dPb);
      dScatter[83] = new ScreenedRutherfordScatteringAngle(*Element::dBi);
      dScatter[84] = new ScreenedRutherfordScatteringAngle(*Element::dPo);
      dScatter[85] = new ScreenedRutherfordScatteringAngle(*Element::dAt);
      dScatter[86] = new ScreenedRutherfordScatteringAngle(*Element::dRn);
      dScatter[87] = new ScreenedRutherfordScatteringAngle(*Element::dFr);
      dScatter[88] = new ScreenedRutherfordScatteringAngle(*Element::dRa);
      dScatter[89] = new ScreenedRutherfordScatteringAngle(*Element::dAc);
      dScatter[90] = new ScreenedRutherfordScatteringAngle(*Element::dTh);
      dScatter[91] = new ScreenedRutherfordScatteringAngle(*Element::dPa);
      dScatter[92] = new ScreenedRutherfordScatteringAngle(*Element::dU);
      dScatter[93] = new ScreenedRutherfordScatteringAngle(*Element::dNp);
      dScatter[94] = new ScreenedRutherfordScatteringAngle(*Element::dPu);
      dScatter[95] = new ScreenedRutherfordScatteringAngle(*Element::dAm);
      dScatter[96] = new ScreenedRutherfordScatteringAngle(*Element::dCm);
   }

   __host__ __device__ const ScreenedRutherfordScatteringAngle& getSRSA(int an)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      return *dScatter[an];
#else
      return *mScatter[an];
#endif
   }

   __host__ __device__ ScreenedRutherfordRandomizedScatterFactory::ScreenedRutherfordRandomizedScatterFactory() :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      RandomizedScatterFactoryT("Screened Rutherford elastic cross-section", *Reference::d_NullReference)
#else
      RandomizedScatterFactoryT("Screened Rutherford elastic cross-section", REFERENCE)
#endif
   {
   }

   __host__ __device__ const RandomizedScatterT& ScreenedRutherfordRandomizedScatterFactory::get(const ElementT& elm) const
   {
      return getSRSA(elm.getAtomicNumber());
   }

   const ScreenedRutherfordRandomizedScatterFactory FactoryScreenedRutherford;
   const RandomizedScatterFactoryT& Factory = FactoryScreenedRutherford;
}
