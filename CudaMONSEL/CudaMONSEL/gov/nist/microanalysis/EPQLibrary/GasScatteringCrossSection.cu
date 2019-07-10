#include "gov\nist\microanalysis\EPQLibrary\GasScatteringCrossSection.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Reference.cuh"
#include "gov\nist\microanalysis\EPQLibrary\PhysicalConstants.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ScreenedRutherfordScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

#include "gov\nist\microanalysis\Utility\Math2.cuh"

#include "Amphibian\random.cuh"

namespace GasScatteringCrossSection
{
   const Reference::Author* al[] = { &Reference::RFEdgerton };
   const Reference::Book REFERENCE("Electron Energy-Loss Spectroscopy in the Electron Microscope, Second Edition", "Plenum Press, NY & London", 1996, al, 1);

   static const double E0 = PhysicalConstants::ElectronMass * PhysicalConstants::SpeedOfLight * PhysicalConstants::SpeedOfLight;

   __host__ __device__ GasScatteringCrossSection::GasScatteringCrossSection(const ElementT& elm) :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      RandomizedScatterT("Edgerton gas cross-section", *Reference::dNullReference), mElement(elm), mElastic(ScreenedRutherfordScatteringAngle::getSRSA(elm.getAtomicNumber()))
#else
      RandomizedScatterT("Edgerton gas cross-section", REFERENCE), mElement(elm), mElastic(ScreenedRutherfordScatteringAngle::getSRSA(elm.getAtomicNumber()))
#endif
   {
   }
   
   const RandomizedScatterT& GasScatteringCrossSection::getElasticModel()
   {
      return mElastic;
   }

   __host__ __device__ double GasScatteringCrossSection::ratioInelasticOverElastic() const
   {
      return 20.0 / mElement.getAtomicNumber();
   }

   const ElementT& GasScatteringCrossSection::getElement() const
   {
      return mElement;
   }

   __host__ __device__ double GasScatteringCrossSection::totalCrossSection(double energy) const
   {
      return (1.0 + ratioInelasticOverElastic()) * mElastic.totalCrossSection(energy);
   }

   double GasScatteringCrossSection::randomScatteringAngle(double energy) const
   {
      if (Random::random() * (1.0 + ratioInelasticOverElastic()) < 1.0)
         return mElastic.randomScatteringAngle(energy);
      else {
         // Electron velocity from energy
         double v = PhysicalConstants::SpeedOfLight * ::sqrt(1.0 - 1.0 / Math2::sqr((energy / E0) + 1.0));
         // Compute gamma, wave vector magnitude and atomic radius
         double g = ::sqrt(1.0 - Math2::sqr(v / PhysicalConstants::SpeedOfLight));
         double k0 = g * PhysicalConstants::ElectronMass * v / PhysicalConstants::PlanckReduced;
         double rz = PhysicalConstants::BohrRadius * ::pow(mElement.getAtomicNumber(), -0.3333333);
         // Two characteristic angles
         double thE2 = Math2::sqr(mElement.getIonizationEnergy() / (g * PhysicalConstants::ElectronMass * v * v)); // unitless
         double th02 = Math2::sqr(1.0 / (k0 * rz));
         // Compute the maximum integrated cross section
         double siInt = ::log((PhysicalConstants::PI * PhysicalConstants::PI + thE2) * (th02 + thE2) / (thE2 * (PhysicalConstants::PI * PhysicalConstants::PI + th02 + thE2)));
         if (!(siInt > 0.0)) printf("GasScatteringCrossSection::randomScatteringAngle siInt %lf\n", siInt);
         // Select a random integrated cross section
         double exp_si = ::exp(Random::random() * siInt);
         if (exp_si < 1.0) printf("GasScatteringCrossSection::randomScatteringAngle exp_si %lf\n", exp_si);
         // Solve for the angle that give us this (via Egerton 3.16)
         double beta = ::sqrt(((1 - exp_si) * thE2 * (th02 + thE2)) / ((exp_si - 1) * thE2 - th02));
         if (!((beta >= 0) && (beta <= PhysicalConstants::PI))) printf("GasScatteringCrossSection::randomScatteringAngle beta %lf\n", beta);
         return beta;
      }
   }

   //const GasScatteringCrossSection GSCS1(Element::H);
   //const GasScatteringCrossSection GSCS2(Element::He);
   //const GasScatteringCrossSection GSCS3(Element::Li);
   //const GasScatteringCrossSection GSCS4(Element::Be);
   //const GasScatteringCrossSection GSCS5(Element::B);
   //const GasScatteringCrossSection GSCS6(Element::C);
   //const GasScatteringCrossSection GSCS7(Element::N);
   //const GasScatteringCrossSection GSCS8(Element::O);
   //const GasScatteringCrossSection GSCS9(Element::F);
   //const GasScatteringCrossSection GSCS10(Element::Ne);
   //const GasScatteringCrossSection GSCS11(Element::Na);
   //const GasScatteringCrossSection GSCS12(Element::Mg);
   //const GasScatteringCrossSection GSCS13(Element::Al);
   //const GasScatteringCrossSection GSCS14(Element::Si);
   //const GasScatteringCrossSection GSCS15(Element::P);
   //const GasScatteringCrossSection GSCS16(Element::S);
   //const GasScatteringCrossSection GSCS17(Element::Cl);
   //const GasScatteringCrossSection GSCS18(Element::Ar);
   //const GasScatteringCrossSection GSCS19(Element::K);
   //const GasScatteringCrossSection GSCS20(Element::Ca);
   //const GasScatteringCrossSection GSCS21(Element::Sc);
   //const GasScatteringCrossSection GSCS22(Element::Ti);
   //const GasScatteringCrossSection GSCS23(Element::V);
   //const GasScatteringCrossSection GSCS24(Element::Cr);
   //const GasScatteringCrossSection GSCS25(Element::Mn);
   //const GasScatteringCrossSection GSCS26(Element::Fe);
   //const GasScatteringCrossSection GSCS27(Element::Co);
   //const GasScatteringCrossSection GSCS28(Element::Ni);
   //const GasScatteringCrossSection GSCS29(Element::Cu);
   //const GasScatteringCrossSection GSCS30(Element::Zn);
   //const GasScatteringCrossSection GSCS31(Element::Ga);
   //const GasScatteringCrossSection GSCS32(Element::Ge);
   //const GasScatteringCrossSection GSCS33(Element::As);
   //const GasScatteringCrossSection GSCS34(Element::Se);
   //const GasScatteringCrossSection GSCS35(Element::Br);
   //const GasScatteringCrossSection GSCS36(Element::Kr);
   //const GasScatteringCrossSection GSCS37(Element::Rb);
   //const GasScatteringCrossSection GSCS38(Element::Sr);
   //const GasScatteringCrossSection GSCS39(Element::Y);
   //const GasScatteringCrossSection GSCS40(Element::Zr);
   //const GasScatteringCrossSection GSCS41(Element::Nb);
   //const GasScatteringCrossSection GSCS42(Element::Mo);
   //const GasScatteringCrossSection GSCS43(Element::Tc);
   //const GasScatteringCrossSection GSCS44(Element::Ru);
   //const GasScatteringCrossSection GSCS45(Element::Rh);
   //const GasScatteringCrossSection GSCS46(Element::Pd);
   //const GasScatteringCrossSection GSCS47(Element::Ag);
   //const GasScatteringCrossSection GSCS48(Element::Cd);
   //const GasScatteringCrossSection GSCS49(Element::In);
   //const GasScatteringCrossSection GSCS50(Element::Sn);
   //const GasScatteringCrossSection GSCS51(Element::Sb);
   //const GasScatteringCrossSection GSCS52(Element::Te);
   //const GasScatteringCrossSection GSCS53(Element::I);
   //const GasScatteringCrossSection GSCS54(Element::Xe);
   //const GasScatteringCrossSection GSCS55(Element::Cs);
   //const GasScatteringCrossSection GSCS56(Element::Ba);
   //const GasScatteringCrossSection GSCS57(Element::La);
   //const GasScatteringCrossSection GSCS58(Element::Ce);
   //const GasScatteringCrossSection GSCS59(Element::Pr);
   //const GasScatteringCrossSection GSCS60(Element::Nd);
   //const GasScatteringCrossSection GSCS61(Element::Pm);
   //const GasScatteringCrossSection GSCS62(Element::Sm);
   //const GasScatteringCrossSection GSCS63(Element::Eu);
   //const GasScatteringCrossSection GSCS64(Element::Gd);
   //const GasScatteringCrossSection GSCS65(Element::Tb);
   //const GasScatteringCrossSection GSCS66(Element::Dy);
   //const GasScatteringCrossSection GSCS67(Element::Ho);
   //const GasScatteringCrossSection GSCS68(Element::Er);
   //const GasScatteringCrossSection GSCS69(Element::Tm);
   //const GasScatteringCrossSection GSCS70(Element::Yb);
   //const GasScatteringCrossSection GSCS71(Element::Lu);
   //const GasScatteringCrossSection GSCS72(Element::Hf);
   //const GasScatteringCrossSection GSCS73(Element::Ta);
   //const GasScatteringCrossSection GSCS74(Element::W);
   //const GasScatteringCrossSection GSCS75(Element::Re);
   //const GasScatteringCrossSection GSCS76(Element::Os);
   //const GasScatteringCrossSection GSCS77(Element::Ir);
   //const GasScatteringCrossSection GSCS78(Element::Pt);
   //const GasScatteringCrossSection GSCS79(Element::Au);
   //const GasScatteringCrossSection GSCS80(Element::Hg);
   //const GasScatteringCrossSection GSCS81(Element::Tl);
   //const GasScatteringCrossSection GSCS82(Element::Pb);
   //const GasScatteringCrossSection GSCS83(Element::Bi);
   //const GasScatteringCrossSection GSCS84(Element::Po);
   //const GasScatteringCrossSection GSCS85(Element::At);
   //const GasScatteringCrossSection GSCS86(Element::Rn);
   //const GasScatteringCrossSection GSCS87(Element::Fr);
   //const GasScatteringCrossSection GSCS88(Element::Ra);
   //const GasScatteringCrossSection GSCS89(Element::Ac);
   //const GasScatteringCrossSection GSCS90(Element::Th);
   //const GasScatteringCrossSection GSCS91(Element::Pa);
   //const GasScatteringCrossSection GSCS92(Element::U);
   //const GasScatteringCrossSection GSCS93(Element::Np);
   //const GasScatteringCrossSection GSCS94(Element::Pu);
   //const GasScatteringCrossSection GSCS95(Element::Am);
   //const GasScatteringCrossSection GSCS96(Element::Cm);

   GasScatteringCrossSection const * mScatter[113];

   void init()
   {
      mScatter[1] = new GasScatteringCrossSection(Element::H);
      mScatter[2] = new GasScatteringCrossSection(Element::He);
      mScatter[3] = new GasScatteringCrossSection(Element::Li);
      mScatter[4] = new GasScatteringCrossSection(Element::Be);
      mScatter[5] = new GasScatteringCrossSection(Element::B);
      mScatter[6] = new GasScatteringCrossSection(Element::C);
      mScatter[7] = new GasScatteringCrossSection(Element::N);
      mScatter[8] = new GasScatteringCrossSection(Element::O);
      mScatter[9] = new GasScatteringCrossSection(Element::F);
      mScatter[10] = new GasScatteringCrossSection(Element::Ne);
      mScatter[11] = new GasScatteringCrossSection(Element::Na);
      mScatter[12] = new GasScatteringCrossSection(Element::Mg);
      mScatter[13] = new GasScatteringCrossSection(Element::Al);
      mScatter[14] = new GasScatteringCrossSection(Element::Si);
      mScatter[15] = new GasScatteringCrossSection(Element::P);
      mScatter[16] = new GasScatteringCrossSection(Element::S);
      mScatter[17] = new GasScatteringCrossSection(Element::Cl);
      mScatter[18] = new GasScatteringCrossSection(Element::Ar);
      mScatter[19] = new GasScatteringCrossSection(Element::K);
      mScatter[20] = new GasScatteringCrossSection(Element::Ca);
      mScatter[21] = new GasScatteringCrossSection(Element::Sc);
      mScatter[22] = new GasScatteringCrossSection(Element::Ti);
      mScatter[23] = new GasScatteringCrossSection(Element::V);
      mScatter[24] = new GasScatteringCrossSection(Element::Cr);
      mScatter[25] = new GasScatteringCrossSection(Element::Mn);
      mScatter[26] = new GasScatteringCrossSection(Element::Fe);
      mScatter[27] = new GasScatteringCrossSection(Element::Co);
      mScatter[28] = new GasScatteringCrossSection(Element::Ni);
      mScatter[29] = new GasScatteringCrossSection(Element::Cu);
      mScatter[30] = new GasScatteringCrossSection(Element::Zn);
      mScatter[31] = new GasScatteringCrossSection(Element::Ga);
      mScatter[32] = new GasScatteringCrossSection(Element::Ge);
      mScatter[33] = new GasScatteringCrossSection(Element::As);
      mScatter[34] = new GasScatteringCrossSection(Element::Se);
      mScatter[35] = new GasScatteringCrossSection(Element::Br);
      mScatter[36] = new GasScatteringCrossSection(Element::Kr);
      mScatter[37] = new GasScatteringCrossSection(Element::Rb);
      mScatter[38] = new GasScatteringCrossSection(Element::Sr);
      mScatter[39] = new GasScatteringCrossSection(Element::Y);
      mScatter[40] = new GasScatteringCrossSection(Element::Zr);
      mScatter[41] = new GasScatteringCrossSection(Element::Nb);
      mScatter[42] = new GasScatteringCrossSection(Element::Mo);
      mScatter[43] = new GasScatteringCrossSection(Element::Tc);
      mScatter[44] = new GasScatteringCrossSection(Element::Ru);
      mScatter[45] = new GasScatteringCrossSection(Element::Rh);
      mScatter[46] = new GasScatteringCrossSection(Element::Pd);
      mScatter[47] = new GasScatteringCrossSection(Element::Ag);
      mScatter[48] = new GasScatteringCrossSection(Element::Cd);
      mScatter[49] = new GasScatteringCrossSection(Element::In);
      mScatter[50] = new GasScatteringCrossSection(Element::Sn);
      mScatter[51] = new GasScatteringCrossSection(Element::Sb);
      mScatter[52] = new GasScatteringCrossSection(Element::Te);
      mScatter[53] = new GasScatteringCrossSection(Element::I);
      mScatter[54] = new GasScatteringCrossSection(Element::Xe);
      mScatter[55] = new GasScatteringCrossSection(Element::Cs);
      mScatter[56] = new GasScatteringCrossSection(Element::Ba);
      mScatter[57] = new GasScatteringCrossSection(Element::La);
      mScatter[58] = new GasScatteringCrossSection(Element::Ce);
      mScatter[59] = new GasScatteringCrossSection(Element::Pr);
      mScatter[60] = new GasScatteringCrossSection(Element::Nd);
      mScatter[61] = new GasScatteringCrossSection(Element::Pm);
      mScatter[62] = new GasScatteringCrossSection(Element::Sm);
      mScatter[63] = new GasScatteringCrossSection(Element::Eu);
      mScatter[64] = new GasScatteringCrossSection(Element::Gd);
      mScatter[65] = new GasScatteringCrossSection(Element::Tb);
      mScatter[66] = new GasScatteringCrossSection(Element::Dy);
      mScatter[67] = new GasScatteringCrossSection(Element::Ho);
      mScatter[68] = new GasScatteringCrossSection(Element::Er);
      mScatter[69] = new GasScatteringCrossSection(Element::Tm);
      mScatter[70] = new GasScatteringCrossSection(Element::Yb);
      mScatter[71] = new GasScatteringCrossSection(Element::Lu);
      mScatter[72] = new GasScatteringCrossSection(Element::Hf);
      mScatter[73] = new GasScatteringCrossSection(Element::Ta);
      mScatter[74] = new GasScatteringCrossSection(Element::W);
      mScatter[75] = new GasScatteringCrossSection(Element::Re);
      mScatter[76] = new GasScatteringCrossSection(Element::Os);
      mScatter[77] = new GasScatteringCrossSection(Element::Ir);
      mScatter[78] = new GasScatteringCrossSection(Element::Pt);
      mScatter[79] = new GasScatteringCrossSection(Element::Au);
      mScatter[80] = new GasScatteringCrossSection(Element::Hg);
      mScatter[81] = new GasScatteringCrossSection(Element::Tl);
      mScatter[82] = new GasScatteringCrossSection(Element::Pb);
      mScatter[83] = new GasScatteringCrossSection(Element::Bi);
      mScatter[84] = new GasScatteringCrossSection(Element::Po);
      mScatter[85] = new GasScatteringCrossSection(Element::At);
      mScatter[86] = new GasScatteringCrossSection(Element::Rn);
      mScatter[87] = new GasScatteringCrossSection(Element::Fr);
      mScatter[88] = new GasScatteringCrossSection(Element::Ra);
      mScatter[89] = new GasScatteringCrossSection(Element::Ac);
      mScatter[90] = new GasScatteringCrossSection(Element::Th);
      mScatter[91] = new GasScatteringCrossSection(Element::Pa);
      mScatter[92] = new GasScatteringCrossSection(Element::U);
      mScatter[93] = new GasScatteringCrossSection(Element::Np);
      mScatter[94] = new GasScatteringCrossSection(Element::Pu);
      mScatter[95] = new GasScatteringCrossSection(Element::Am);
      mScatter[96] = new GasScatteringCrossSection(Element::Cm);
   }

   __device__ GasScatteringCrossSection const * dScatter[113];

   __global__ void initCuda()
   {
      dScatter[1] = new GasScatteringCrossSection(*Element::dH);
      dScatter[2] = new GasScatteringCrossSection(*Element::dHe);
      dScatter[3] = new GasScatteringCrossSection(*Element::dLi);
      dScatter[4] = new GasScatteringCrossSection(*Element::dBe);
      dScatter[5] = new GasScatteringCrossSection(*Element::dB);
      dScatter[6] = new GasScatteringCrossSection(*Element::dC);
      dScatter[7] = new GasScatteringCrossSection(*Element::dN);
      dScatter[8] = new GasScatteringCrossSection(*Element::dO);
      dScatter[9] = new GasScatteringCrossSection(*Element::dF);
      dScatter[10] = new GasScatteringCrossSection(*Element::dNe);
      dScatter[11] = new GasScatteringCrossSection(*Element::dNa);
      dScatter[12] = new GasScatteringCrossSection(*Element::dMg);
      dScatter[13] = new GasScatteringCrossSection(*Element::dAl);
      dScatter[14] = new GasScatteringCrossSection(*Element::dSi);
      dScatter[15] = new GasScatteringCrossSection(*Element::dP);
      dScatter[16] = new GasScatteringCrossSection(*Element::dS);
      dScatter[17] = new GasScatteringCrossSection(*Element::dCl);
      dScatter[18] = new GasScatteringCrossSection(*Element::dAr);
      dScatter[19] = new GasScatteringCrossSection(*Element::dK);
      dScatter[20] = new GasScatteringCrossSection(*Element::dCa);
      dScatter[21] = new GasScatteringCrossSection(*Element::dSc);
      dScatter[22] = new GasScatteringCrossSection(*Element::dTi);
      dScatter[23] = new GasScatteringCrossSection(*Element::dV);
      dScatter[24] = new GasScatteringCrossSection(*Element::dCr);
      dScatter[25] = new GasScatteringCrossSection(*Element::dMn);
      dScatter[26] = new GasScatteringCrossSection(*Element::dFe);
      dScatter[27] = new GasScatteringCrossSection(*Element::dCo);
      dScatter[28] = new GasScatteringCrossSection(*Element::dNi);
      dScatter[29] = new GasScatteringCrossSection(*Element::dCu);
      dScatter[30] = new GasScatteringCrossSection(*Element::dZn);
      dScatter[31] = new GasScatteringCrossSection(*Element::dGa);
      dScatter[32] = new GasScatteringCrossSection(*Element::dGe);
      dScatter[33] = new GasScatteringCrossSection(*Element::dAs);
      dScatter[34] = new GasScatteringCrossSection(*Element::dSe);
      dScatter[35] = new GasScatteringCrossSection(*Element::dBr);
      dScatter[36] = new GasScatteringCrossSection(*Element::dKr);
      dScatter[37] = new GasScatteringCrossSection(*Element::dRb);
      dScatter[38] = new GasScatteringCrossSection(*Element::dSr);
      dScatter[39] = new GasScatteringCrossSection(*Element::dY);
      dScatter[40] = new GasScatteringCrossSection(*Element::dZr);
      dScatter[41] = new GasScatteringCrossSection(*Element::dNb);
      dScatter[42] = new GasScatteringCrossSection(*Element::dMo);
      dScatter[43] = new GasScatteringCrossSection(*Element::dTc);
      dScatter[44] = new GasScatteringCrossSection(*Element::dRu);
      dScatter[45] = new GasScatteringCrossSection(*Element::dRh);
      dScatter[46] = new GasScatteringCrossSection(*Element::dPd);
      dScatter[47] = new GasScatteringCrossSection(*Element::dAg);
      dScatter[48] = new GasScatteringCrossSection(*Element::dCd);
      dScatter[49] = new GasScatteringCrossSection(*Element::dIn);
      dScatter[50] = new GasScatteringCrossSection(*Element::dSn);
      dScatter[51] = new GasScatteringCrossSection(*Element::dSb);
      dScatter[52] = new GasScatteringCrossSection(*Element::dTe);
      dScatter[53] = new GasScatteringCrossSection(*Element::dI);
      dScatter[54] = new GasScatteringCrossSection(*Element::dXe);
      dScatter[55] = new GasScatteringCrossSection(*Element::dCs);
      dScatter[56] = new GasScatteringCrossSection(*Element::dBa);
      dScatter[57] = new GasScatteringCrossSection(*Element::dLa);
      dScatter[58] = new GasScatteringCrossSection(*Element::dCe);
      dScatter[59] = new GasScatteringCrossSection(*Element::dPr);
      dScatter[60] = new GasScatteringCrossSection(*Element::dNd);
      dScatter[61] = new GasScatteringCrossSection(*Element::dPm);
      dScatter[62] = new GasScatteringCrossSection(*Element::dSm);
      dScatter[63] = new GasScatteringCrossSection(*Element::dEu);
      dScatter[64] = new GasScatteringCrossSection(*Element::dGd);
      dScatter[65] = new GasScatteringCrossSection(*Element::dTb);
      dScatter[66] = new GasScatteringCrossSection(*Element::dDy);
      dScatter[67] = new GasScatteringCrossSection(*Element::dHo);
      dScatter[68] = new GasScatteringCrossSection(*Element::dEr);
      dScatter[69] = new GasScatteringCrossSection(*Element::dTm);
      dScatter[70] = new GasScatteringCrossSection(*Element::dYb);
      dScatter[71] = new GasScatteringCrossSection(*Element::dLu);
      dScatter[72] = new GasScatteringCrossSection(*Element::dHf);
      dScatter[73] = new GasScatteringCrossSection(*Element::dTa);
      dScatter[74] = new GasScatteringCrossSection(*Element::dW);
      dScatter[75] = new GasScatteringCrossSection(*Element::dRe);
      dScatter[76] = new GasScatteringCrossSection(*Element::dOs);
      dScatter[77] = new GasScatteringCrossSection(*Element::dIr);
      dScatter[78] = new GasScatteringCrossSection(*Element::dPt);
      dScatter[79] = new GasScatteringCrossSection(*Element::dAu);
      dScatter[80] = new GasScatteringCrossSection(*Element::dHg);
      dScatter[81] = new GasScatteringCrossSection(*Element::dTl);
      dScatter[82] = new GasScatteringCrossSection(*Element::dPb);
      dScatter[83] = new GasScatteringCrossSection(*Element::dBi);
      dScatter[84] = new GasScatteringCrossSection(*Element::dPo);
      dScatter[85] = new GasScatteringCrossSection(*Element::dAt);
      dScatter[86] = new GasScatteringCrossSection(*Element::dRn);
      dScatter[87] = new GasScatteringCrossSection(*Element::dFr);
      dScatter[88] = new GasScatteringCrossSection(*Element::dRa);
      dScatter[89] = new GasScatteringCrossSection(*Element::dAc);
      dScatter[90] = new GasScatteringCrossSection(*Element::dTh);
      dScatter[91] = new GasScatteringCrossSection(*Element::dPa);
      dScatter[92] = new GasScatteringCrossSection(*Element::dU);
      dScatter[93] = new GasScatteringCrossSection(*Element::dNp);
      dScatter[94] = new GasScatteringCrossSection(*Element::dPu);
      dScatter[95] = new GasScatteringCrossSection(*Element::dAm);
      dScatter[96] = new GasScatteringCrossSection(*Element::dCm);
   }

   __host__ __device__ static const GasScatteringCrossSection& getGSCS(int an)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      return *dScatter[an];
#else
      return *mScatter[an];
#endif
   }

   __host__ __device__ GasScatteringRandomizedScatterFactory::GasScatteringRandomizedScatterFactory() :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      RandomizedScatterFactoryT("Gas scattering algorithm", *Reference::dNullReference)
#else
      RandomizedScatterFactoryT("Gas scattering algorithm", REFERENCE)
#endif
   {
   }

   __host__ __device__ const RandomizedScatterT& GasScatteringRandomizedScatterFactory::get(const ElementT& elm) const
   {
      return getGSCS(elm.getAtomicNumber());
   }

   void GasScatteringRandomizedScatterFactory::initializeDefaultStrategy()
   {
   }

   const GasScatteringRandomizedScatterFactory FactoryGasScattering;
   const RandomizedScatterFactoryT& Factory = FactoryGasScattering;
}