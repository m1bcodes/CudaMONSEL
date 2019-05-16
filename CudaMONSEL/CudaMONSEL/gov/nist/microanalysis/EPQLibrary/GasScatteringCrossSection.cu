#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\EPQLibrary\GasScatteringCrossSection.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Reference.cuh"
#include "gov\nist\microanalysis\EPQLibrary\PhysicalConstants.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ScreenedRutherfordScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

namespace GasScatteringCrossSection
{
   const Reference::Author* al[] = { &Reference::RFEdgerton };
   const Reference::Book REFERENCE("Electron Energy-Loss Spectroscopy in the Electron Microscope, Second Edition", "Plenum Press, NY & London", 1996, al, 1);

   static const double E0 = PhysicalConstants::ElectronMass * PhysicalConstants::SpeedOfLight * PhysicalConstants::SpeedOfLight;

   GasScatteringCrossSection::GasScatteringCrossSection(const ElementT& elm) : RandomizedScatterT("Edgerton gas cross-section", REFERENCE), mElement(elm), mElastic(ScreenedRutherfordScatteringAngle::getSRSA(elm.getAtomicNumber()))
   {
   }

   GasScatteringCrossSection::GasScatteringCrossSection(const GasScatteringCrossSection& other) : RandomizedScatterT("Edgerton gas cross-section", REFERENCE), mElement(other.mElement), mElastic(other.mElastic)
   {
   }
   
   const RandomizedScatterT& GasScatteringCrossSection::getElasticModel()
   {
      return mElastic;
   }

   double GasScatteringCrossSection::ratioInelasticOverElastic() const
   {
      return 20.0 / mElement.getAtomicNumber();
   }

   const ElementT& GasScatteringCrossSection::getElement() const
   {
      return mElement;
   }

   double GasScatteringCrossSection::totalCrossSection(double energy) const
   {
      return (1.0 + ratioInelasticOverElastic()) * mElastic.totalCrossSection(energy);
   }

   double GasScatteringCrossSection::randomScatteringAngle(double energy) const
   {
      if (Math2::random() * (1.0 + ratioInelasticOverElastic()) < 1.0)
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
         double exp_si = ::exp(Math2::random() * siInt);
         if (exp_si < 1.0) printf("GasScatteringCrossSection::randomScatteringAngle exp_si %lf\n", exp_si);
         // Solve for the angle that give us this (via Egerton 3.16)
         double beta = ::sqrt(((1 - exp_si) * thE2 * (th02 + thE2)) / ((exp_si - 1) * thE2 - th02));
         if (!((beta >= 0) && (beta <= PhysicalConstants::PI))) printf("GasScatteringCrossSection::randomScatteringAngle beta %lf\n", beta);
         return beta;
      }
   }

   const GasScatteringCrossSection GSCS1(Element::H);
   const GasScatteringCrossSection GSCS2(Element::He);
   const GasScatteringCrossSection GSCS3(Element::Li);
   const GasScatteringCrossSection GSCS4(Element::Be);
   const GasScatteringCrossSection GSCS5(Element::B);
   const GasScatteringCrossSection GSCS6(Element::C);
   const GasScatteringCrossSection GSCS7(Element::N);
   const GasScatteringCrossSection GSCS8(Element::O);
   const GasScatteringCrossSection GSCS9(Element::F);
   const GasScatteringCrossSection GSCS10(Element::Ne);
   const GasScatteringCrossSection GSCS11(Element::Na);
   const GasScatteringCrossSection GSCS12(Element::Mg);
   const GasScatteringCrossSection GSCS13(Element::Al);
   const GasScatteringCrossSection GSCS14(Element::Si);
   const GasScatteringCrossSection GSCS15(Element::P);
   const GasScatteringCrossSection GSCS16(Element::S);
   const GasScatteringCrossSection GSCS17(Element::Cl);
   const GasScatteringCrossSection GSCS18(Element::Ar);
   const GasScatteringCrossSection GSCS19(Element::K);
   const GasScatteringCrossSection GSCS20(Element::Ca);
   const GasScatteringCrossSection GSCS21(Element::Sc);
   const GasScatteringCrossSection GSCS22(Element::Ti);
   const GasScatteringCrossSection GSCS23(Element::V);
   const GasScatteringCrossSection GSCS24(Element::Cr);
   const GasScatteringCrossSection GSCS25(Element::Mn);
   const GasScatteringCrossSection GSCS26(Element::Fe);
   const GasScatteringCrossSection GSCS27(Element::Co);
   const GasScatteringCrossSection GSCS28(Element::Ni);
   const GasScatteringCrossSection GSCS29(Element::Cu);
   const GasScatteringCrossSection GSCS30(Element::Zn);
   const GasScatteringCrossSection GSCS31(Element::Ga);
   const GasScatteringCrossSection GSCS32(Element::Ge);
   const GasScatteringCrossSection GSCS33(Element::As);
   const GasScatteringCrossSection GSCS34(Element::Se);
   const GasScatteringCrossSection GSCS35(Element::Br);
   const GasScatteringCrossSection GSCS36(Element::Kr);
   const GasScatteringCrossSection GSCS37(Element::Rb);
   const GasScatteringCrossSection GSCS38(Element::Sr);
   const GasScatteringCrossSection GSCS39(Element::Y);
   const GasScatteringCrossSection GSCS40(Element::Zr);
   const GasScatteringCrossSection GSCS41(Element::Nb);
   const GasScatteringCrossSection GSCS42(Element::Mo);
   const GasScatteringCrossSection GSCS43(Element::Tc);
   const GasScatteringCrossSection GSCS44(Element::Ru);
   const GasScatteringCrossSection GSCS45(Element::Rh);
   const GasScatteringCrossSection GSCS46(Element::Pd);
   const GasScatteringCrossSection GSCS47(Element::Ag);
   const GasScatteringCrossSection GSCS48(Element::Cd);
   const GasScatteringCrossSection GSCS49(Element::In);
   const GasScatteringCrossSection GSCS50(Element::Sn);
   const GasScatteringCrossSection GSCS51(Element::Sb);
   const GasScatteringCrossSection GSCS52(Element::Te);
   const GasScatteringCrossSection GSCS53(Element::I);
   const GasScatteringCrossSection GSCS54(Element::Xe);
   const GasScatteringCrossSection GSCS55(Element::Cs);
   const GasScatteringCrossSection GSCS56(Element::Ba);
   const GasScatteringCrossSection GSCS57(Element::La);
   const GasScatteringCrossSection GSCS58(Element::Ce);
   const GasScatteringCrossSection GSCS59(Element::Pr);
   const GasScatteringCrossSection GSCS60(Element::Nd);
   const GasScatteringCrossSection GSCS61(Element::Pm);
   const GasScatteringCrossSection GSCS62(Element::Sm);
   const GasScatteringCrossSection GSCS63(Element::Eu);
   const GasScatteringCrossSection GSCS64(Element::Gd);
   const GasScatteringCrossSection GSCS65(Element::Tb);
   const GasScatteringCrossSection GSCS66(Element::Dy);
   const GasScatteringCrossSection GSCS67(Element::Ho);
   const GasScatteringCrossSection GSCS68(Element::Er);
   const GasScatteringCrossSection GSCS69(Element::Tm);
   const GasScatteringCrossSection GSCS70(Element::Yb);
   const GasScatteringCrossSection GSCS71(Element::Lu);
   const GasScatteringCrossSection GSCS72(Element::Hf);
   const GasScatteringCrossSection GSCS73(Element::Ta);
   const GasScatteringCrossSection GSCS74(Element::W);
   const GasScatteringCrossSection GSCS75(Element::Re);
   const GasScatteringCrossSection GSCS76(Element::Os);
   const GasScatteringCrossSection GSCS77(Element::Ir);
   const GasScatteringCrossSection GSCS78(Element::Pt);
   const GasScatteringCrossSection GSCS79(Element::Au);
   const GasScatteringCrossSection GSCS80(Element::Hg);
   const GasScatteringCrossSection GSCS81(Element::Tl);
   const GasScatteringCrossSection GSCS82(Element::Pb);
   const GasScatteringCrossSection GSCS83(Element::Bi);
   const GasScatteringCrossSection GSCS84(Element::Po);
   const GasScatteringCrossSection GSCS85(Element::At);
   const GasScatteringCrossSection GSCS86(Element::Rn);
   const GasScatteringCrossSection GSCS87(Element::Fr);
   const GasScatteringCrossSection GSCS88(Element::Ra);
   const GasScatteringCrossSection GSCS89(Element::Ac);
   const GasScatteringCrossSection GSCS90(Element::Th);
   const GasScatteringCrossSection GSCS91(Element::Pa);
   const GasScatteringCrossSection GSCS92(Element::U);
   const GasScatteringCrossSection GSCS93(Element::Np);
   const GasScatteringCrossSection GSCS94(Element::Pu);
   const GasScatteringCrossSection GSCS95(Element::Am);
   const GasScatteringCrossSection GSCS96(Element::Cm);

   GasScatteringCrossSection const * mScatter[113] = {
      nullptr,
      &GSCS1,
      &GSCS2,
      &GSCS3,
      &GSCS4,
      &GSCS5,
      &GSCS6,
      &GSCS7,
      &GSCS8,
      &GSCS9,
      &GSCS10,
      &GSCS11,
      &GSCS12,
      &GSCS13,
      &GSCS14,
      &GSCS15,
      &GSCS16,
      &GSCS17,
      &GSCS18,
      &GSCS19,
      &GSCS20,
      &GSCS21,
      &GSCS22,
      &GSCS23,
      &GSCS24,
      &GSCS25,
      &GSCS26,
      &GSCS27,
      &GSCS28,
      &GSCS29,
      &GSCS30,
      &GSCS31,
      &GSCS32,
      &GSCS33,
      &GSCS34,
      &GSCS35,
      &GSCS36,
      &GSCS37,
      &GSCS38,
      &GSCS39,
      &GSCS40,
      &GSCS41,
      &GSCS42,
      &GSCS43,
      &GSCS44,
      &GSCS45,
      &GSCS46,
      &GSCS47,
      &GSCS48,
      &GSCS49,
      &GSCS50,
      &GSCS51,
      &GSCS52,
      &GSCS53,
      &GSCS54,
      &GSCS55,
      &GSCS56,
      &GSCS57,
      &GSCS58,
      &GSCS59,
      &GSCS60,
      &GSCS61,
      &GSCS62,
      &GSCS63,
      &GSCS64,
      &GSCS65,
      &GSCS66,
      &GSCS67,
      &GSCS68,
      &GSCS69,
      &GSCS70,
      &GSCS71,
      &GSCS72,
      &GSCS73,
      &GSCS74,
      &GSCS75,
      &GSCS76,
      &GSCS77,
      &GSCS78,
      &GSCS79,
      &GSCS80,
      &GSCS81,
      &GSCS82,
      &GSCS83,
      &GSCS84,
      &GSCS85,
      &GSCS86,
      &GSCS87,
      &GSCS88,
      &GSCS89,
      &GSCS90,
      &GSCS91,
      &GSCS92,
      &GSCS93,
      &GSCS94,
      &GSCS95,
      &GSCS96
   };

   const GasScatteringCrossSection& getGSCS(int an)
   {
      return *mScatter[an];
   }

   GasScatteringRandomizedScatterFactory::GasScatteringRandomizedScatterFactory() : RandomizedScatterFactoryT("Gas scattering algorithm", REFERENCE)
   {
   }

   const RandomizedScatterT& GasScatteringRandomizedScatterFactory::get(const ElementT& elm) const
   {
      return getGSCS(elm.getAtomicNumber());
   }

   void GasScatteringRandomizedScatterFactory::initializeDefaultStrategy()
   {
   }

   const GasScatteringRandomizedScatterFactory FactoryGasScattering;
   const RandomizedScatterFactoryT& Factory = FactoryGasScattering;
}