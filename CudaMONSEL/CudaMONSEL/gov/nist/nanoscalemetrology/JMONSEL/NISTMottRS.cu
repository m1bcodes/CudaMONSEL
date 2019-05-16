#include "gov\nist\nanoscalemetrology\JMONSEL\NISTMottRS.cuh"

#include "gov\nist\microanalysis\EPQLibrary\Reference.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\PhysicalConstants.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ScreenedRutherfordScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\BrowningEmpiricalCrossSection.cuh"
#include "gov\nist\microanalysis\EPQLibrary\NISTMottScatteringAngle.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

#include "gov\nist\nanoscalemetrology\JMONSELutils\ULagrangeInterpolation.cuh"

namespace NISTMottRS
{
   const double MAX_NISTMOTT = ToSI::keV(20.0);
   const double MIN_NISTMOTT = ToSI::keV(0.050);

   static const int qINTERPOLATIONORDER = 3;
   static const int sigmaINTERPOLATIONORDER = 3;
   static const double scale = PhysicalConstants::BohrRadius * PhysicalConstants::BohrRadius;

   static const int SPWEM_LEN = 61;
   static const int X1_LEN = 201;
   static const double DL50 = ::log(MIN_NISTMOTT);
   static const double PARAM = (::log(MAX_NISTMOTT) - DL50) / 60.0;

   static const Reference::Author* al[] = {
      &Reference::FSalvat,
      &Reference::AJablonski,
      &Reference::CPowell
   };
   static const Reference::WebSite mReferenceWebsite("http://www.nist.gov/srd/nist64.htm", "NIST Electron Elastic-Scattering Cross-Section Database version 3.1", "AUGUST 24, 2007", al, 3);

   static double sciToDub(const std::string& str)
   {
      std::string tmp = str.substr(str.find_first_not_of(" "));
      std::stringstream ss(tmp);
      double d = 0;
      ss >> d;

      if (ss.fail()) {
         std::string s = "Unable to format ";
         s += tmp;
         s += " as a number!";
         throw s;
      }

      return d;
   }

   //void NISTMottRS::loadData(int an)
   //{
   //   const std::string name(an < 10 ? ".\\gov\\nist\\microanalysis\\EPQLibrary\\NistXSec/E0" + std::to_string(an) + ".D64" : ".\\gov\\nist\\microanalysis\\EPQLibrary\\NistXSec/E" + std::to_string(an) + ".D64");
   //   printf("Reading: %s\n", name.c_str());
   //   try {
   //      std::ifstream t(name);
   //      if (!t.good()) throw 0;
   //      std::string line;
   //      std::getline(t, line);
   //      for (int j = 0; j < SPWEM_LEN; ++j) {
   //         std::getline(t, line);
   //         mSpwem[j] = sciToDub(line);
   //         for (int i = 0; i < X1_LEN; ++i){
   //            std::getline(t, line);
   //            mX1[j][i] = sciToDub(line);
   //         }
   //      }
   //      t.close();
   //   }
   //   catch (std::exception& ex) {
   //      printf("Unable to construct NISTMottRS: %s\n", name.c_str());
   //   }
   //}

   NISTMottRS::NISTMottRS(const ElementT& elm, int method) :
      RandomizedScatterT("NIST Elastic cross-section", mReferenceWebsite),
      mElement(elm),
      method(method),
      mRutherford(ScreenedRutherfordScatteringAngle::getSRSA(elm.getAtomicNumber())),
      mBrowning(BrowningEmpiricalCrossSection::getBECS(elm.getAtomicNumber())),
      extrapolateBelowEnergy(method == 1 ? ToSI::eV(50.) : ToSI::eV(100.)),
      mSpwem(NISTMottScatteringAngle::getNISTMSA(elm.getAtomicNumber()).getSpwem()),
      mX1(NISTMottScatteringAngle::getNISTMSA(elm.getAtomicNumber()).getX1()),
      MottXSatMinEnergy(totalCrossSection(extrapolateBelowEnergy)),
      sfBrowning(MottXSatMinEnergy / mBrowning.totalCrossSection(extrapolateBelowEnergy))
   {
      //loadData(elm.getAtomicNumber());
      if (!(method >= 1 && method <= 3))
         printf("NISTMottRS::setMethod: Invalid NISTMottRS method: method must = 1, 2, or 3.");
   }

   StringT NISTMottRS::toString() const
   {
      return "CrossSection[NIST-Mott," + StringT(mElement.toAbbrev()) + "]";
   }

   const ElementT& NISTMottRS::getElement() const
   {
      return mElement;
   }

   double NISTMottRS::totalCrossSection(double energy) const
   {
      if (energy < extrapolateBelowEnergy) {
         if (method == 3) { // linear interpolation
            return MottXSatMinEnergy * energy / extrapolateBelowEnergy;
         }
         else { // Browning interpolation
            return sfBrowning * mBrowning.totalCrossSection(energy);
         }
      }
      else if (energy < MAX_NISTMOTT) {
         return scale * ULagrangeInterpolation::d1(mSpwem, DL50, PARAM, sigmaINTERPOLATIONORDER, ::log(energy))[0];
      }
      else {
         return mRutherford.totalCrossSection(energy);
      }
   }

   double NISTMottRS::randomScatteringAngle(double energy) const
   {
      /*
      * Even in method 3 (linear interpolation) we use Browning for the angular
      * distribution.
      */
      if (energy < extrapolateBelowEnergy) {
         //if (mBrowning == null) // never happens
         //   mBrowning = new BrowningEmpiricalCrossSection(mElement);
         // sfBrowning = this.totalCrossSection(MIN_NISTMOTT) / mBrowning.totalCrossSection(MIN_NISTMOTT);
         return mBrowning.randomScatteringAngle(energy);

      }
      else if (energy < MAX_NISTMOTT) {
         const double x0[] =  {
            DL50,
               0.
         };
         const double xinc[] = {
            PARAM,
            0.005
         };
         const double x[] = {
            ::log(energy),
            Math2::random()
         };
         const double q = ULagrangeInterpolation::d2(mX1, x0, 2, xinc, 2, qINTERPOLATIONORDER, x, 2)[0];
         const double com = 1.0 - (2.0 * q * q);
         return com > -1.0 ? (com < 1.0 ? ::acos(com) : 0.0) : Math2::PI;
      }
      else {
         return mRutherford.randomScatteringAngle(energy);
      }
   }

   int NISTMottRS::getMethod() const
   {
      return method;
   }

   //void NISTMottRS::setMethod(int method)
   //{
   //   if (method >= 1 && method <= 3)
   //      this->method = method;
   //   else
   //      printf("NISTMottRS::setMethod: Invalid NISTMottRS method: method must = 1, 2, or 3.");
   //   if (method == 2 || method == 3)
   //      extrapolateBelowEnergy = ToSI::eV(100.);
   //   MottXSatMinEnergy = totalCrossSection(extrapolateBelowEnergy);
   //}

   //const NISTMottRS NMRS1_1(Element::H, 1);
   //const NISTMottRS NMRS2_1(Element::He, 1);
   //const NISTMottRS NMRS3_1(Element::Li, 1);
   //const NISTMottRS NMRS4_1(Element::Be, 1);
   //const NISTMottRS NMRS5_1(Element::B, 1);
   //const NISTMottRS NMRS6_1(Element::C, 1);
   //const NISTMottRS NMRS7_1(Element::N, 1);
   //const NISTMottRS NMRS8_1(Element::O, 1);
   //const NISTMottRS NMRS9_1(Element::F, 1);
   //const NISTMottRS NMRS10_1(Element::Ne, 1);
   //const NISTMottRS NMRS11_1(Element::Na, 1);
   //const NISTMottRS NMRS12_1(Element::Mg, 1);
   //const NISTMottRS NMRS13_1(Element::Al, 1);
   //const NISTMottRS NMRS14_1(Element::Si, 1);
   //const NISTMottRS NMRS15_1(Element::P, 1);
   //const NISTMottRS NMRS16_1(Element::S, 1);
   //const NISTMottRS NMRS17_1(Element::Cl, 1);
   //const NISTMottRS NMRS18_1(Element::Ar, 1);
   //const NISTMottRS NMRS19_1(Element::K, 1);
   //const NISTMottRS NMRS20_1(Element::Ca, 1);
   //const NISTMottRS NMRS21_1(Element::Sc, 1);
   //const NISTMottRS NMRS22_1(Element::Ti, 1);
   //const NISTMottRS NMRS23_1(Element::V, 1);
   //const NISTMottRS NMRS24_1(Element::Cr, 1);
   //const NISTMottRS NMRS25_1(Element::Mn, 1);
   //const NISTMottRS NMRS26_1(Element::Fe, 1);
   //const NISTMottRS NMRS27_1(Element::Co, 1);
   //const NISTMottRS NMRS28_1(Element::Ni, 1);
   //const NISTMottRS NMRS29_1(Element::Cu, 1);
   //const NISTMottRS NMRS30_1(Element::Zn, 1);
   //const NISTMottRS NMRS31_1(Element::Ga, 1);
   //const NISTMottRS NMRS32_1(Element::Ge, 1);
   //const NISTMottRS NMRS33_1(Element::As, 1);
   //const NISTMottRS NMRS34_1(Element::Se, 1);
   //const NISTMottRS NMRS35_1(Element::Br, 1);
   //const NISTMottRS NMRS36_1(Element::Kr, 1);
   //const NISTMottRS NMRS37_1(Element::Rb, 1);
   //const NISTMottRS NMRS38_1(Element::Sr, 1);
   //const NISTMottRS NMRS39_1(Element::Y, 1);
   //const NISTMottRS NMRS40_1(Element::Zr, 1);
   //const NISTMottRS NMRS41_1(Element::Nb, 1);
   //const NISTMottRS NMRS42_1(Element::Mo, 1);
   //const NISTMottRS NMRS43_1(Element::Tc, 1);
   //const NISTMottRS NMRS44_1(Element::Ru, 1);
   //const NISTMottRS NMRS45_1(Element::Rh, 1);
   //const NISTMottRS NMRS46_1(Element::Pd, 1);
   //const NISTMottRS NMRS47_1(Element::Ag, 1);
   //const NISTMottRS NMRS48_1(Element::Cd, 1);
   //const NISTMottRS NMRS49_1(Element::In, 1);
   //const NISTMottRS NMRS50_1(Element::Sn, 1);
   //const NISTMottRS NMRS51_1(Element::Sb, 1);
   //const NISTMottRS NMRS52_1(Element::Te, 1);
   //const NISTMottRS NMRS53_1(Element::I, 1);
   //const NISTMottRS NMRS54_1(Element::Xe, 1);
   //const NISTMottRS NMRS55_1(Element::Cs, 1);
   //const NISTMottRS NMRS56_1(Element::Ba, 1);
   //const NISTMottRS NMRS57_1(Element::La, 1);
   //const NISTMottRS NMRS58_1(Element::Ce, 1);
   //const NISTMottRS NMRS59_1(Element::Pr, 1);
   //const NISTMottRS NMRS60_1(Element::Nd, 1);
   //const NISTMottRS NMRS61_1(Element::Pm, 1);
   //const NISTMottRS NMRS62_1(Element::Sm, 1);
   //const NISTMottRS NMRS63_1(Element::Eu, 1);
   //const NISTMottRS NMRS64_1(Element::Gd, 1);
   //const NISTMottRS NMRS65_1(Element::Tb, 1);
   //const NISTMottRS NMRS66_1(Element::Dy, 1);
   //const NISTMottRS NMRS67_1(Element::Ho, 1);
   //const NISTMottRS NMRS68_1(Element::Er, 1);
   //const NISTMottRS NMRS69_1(Element::Tm, 1);
   //const NISTMottRS NMRS70_1(Element::Yb, 1);
   //const NISTMottRS NMRS71_1(Element::Lu, 1);
   //const NISTMottRS NMRS72_1(Element::Hf, 1);
   //const NISTMottRS NMRS73_1(Element::Ta, 1);
   //const NISTMottRS NMRS74_1(Element::W, 1);
   //const NISTMottRS NMRS75_1(Element::Re, 1);
   //const NISTMottRS NMRS76_1(Element::Os, 1);
   //const NISTMottRS NMRS77_1(Element::Ir, 1);
   //const NISTMottRS NMRS78_1(Element::Pt, 1);
   //const NISTMottRS NMRS79_1(Element::Au, 1);
   //const NISTMottRS NMRS80_1(Element::Hg, 1);
   //const NISTMottRS NMRS81_1(Element::Tl, 1);
   //const NISTMottRS NMRS82_1(Element::Pb, 1);
   //const NISTMottRS NMRS83_1(Element::Bi, 1);
   //const NISTMottRS NMRS84_1(Element::Po, 1);
   //const NISTMottRS NMRS85_1(Element::At, 1);
   //const NISTMottRS NMRS86_1(Element::Rn, 1);
   //const NISTMottRS NMRS87_1(Element::Fr, 1);
   //const NISTMottRS NMRS88_1(Element::Ra, 1);
   //const NISTMottRS NMRS89_1(Element::Ac, 1);
   //const NISTMottRS NMRS90_1(Element::Th, 1);
   //const NISTMottRS NMRS91_1(Element::Pa, 1);
   //const NISTMottRS NMRS92_1(Element::U, 1);
   //const NISTMottRS NMRS93_1(Element::Np, 1);
   //const NISTMottRS NMRS94_1(Element::Pu, 1);
   //const NISTMottRS NMRS95_1(Element::Am, 1);
   //const NISTMottRS NMRS96_1(Element::Cm, 1);

   const NISTMottRS* NMRS1_1;
   const NISTMottRS* NMRS2_1;
   const NISTMottRS* NMRS3_1;
   const NISTMottRS* NMRS4_1;
   const NISTMottRS* NMRS5_1;
   const NISTMottRS* NMRS6_1;
   const NISTMottRS* NMRS7_1;
   const NISTMottRS* NMRS8_1;
   const NISTMottRS* NMRS9_1;
   const NISTMottRS* NMRS10_1;
   const NISTMottRS* NMRS11_1;
   const NISTMottRS* NMRS12_1;
   const NISTMottRS* NMRS13_1;
   const NISTMottRS* NMRS14_1;
   const NISTMottRS* NMRS15_1;
   const NISTMottRS* NMRS16_1;
   const NISTMottRS* NMRS17_1;
   const NISTMottRS* NMRS18_1;
   const NISTMottRS* NMRS19_1;
   const NISTMottRS* NMRS20_1;
   const NISTMottRS* NMRS21_1;
   const NISTMottRS* NMRS22_1;
   const NISTMottRS* NMRS23_1;
   const NISTMottRS* NMRS24_1;
   const NISTMottRS* NMRS25_1;
   const NISTMottRS* NMRS26_1;
   const NISTMottRS* NMRS27_1;
   const NISTMottRS* NMRS28_1;
   const NISTMottRS* NMRS29_1;
   const NISTMottRS* NMRS30_1;
   const NISTMottRS* NMRS31_1;
   const NISTMottRS* NMRS32_1;
   const NISTMottRS* NMRS33_1;
   const NISTMottRS* NMRS34_1;
   const NISTMottRS* NMRS35_1;
   const NISTMottRS* NMRS36_1;
   const NISTMottRS* NMRS37_1;
   const NISTMottRS* NMRS38_1;
   const NISTMottRS* NMRS39_1;
   const NISTMottRS* NMRS40_1;
   const NISTMottRS* NMRS41_1;
   const NISTMottRS* NMRS42_1;
   const NISTMottRS* NMRS43_1;
   const NISTMottRS* NMRS44_1;
   const NISTMottRS* NMRS45_1;
   const NISTMottRS* NMRS46_1;
   const NISTMottRS* NMRS47_1;
   const NISTMottRS* NMRS48_1;
   const NISTMottRS* NMRS49_1;
   const NISTMottRS* NMRS50_1;
   const NISTMottRS* NMRS51_1;
   const NISTMottRS* NMRS52_1;
   const NISTMottRS* NMRS53_1;
   const NISTMottRS* NMRS54_1;
   const NISTMottRS* NMRS55_1;
   const NISTMottRS* NMRS56_1;
   const NISTMottRS* NMRS57_1;
   const NISTMottRS* NMRS58_1;
   const NISTMottRS* NMRS59_1;
   const NISTMottRS* NMRS60_1;
   const NISTMottRS* NMRS61_1;
   const NISTMottRS* NMRS62_1;
   const NISTMottRS* NMRS63_1;
   const NISTMottRS* NMRS64_1;
   const NISTMottRS* NMRS65_1;
   const NISTMottRS* NMRS66_1;
   const NISTMottRS* NMRS67_1;
   const NISTMottRS* NMRS68_1;
   const NISTMottRS* NMRS69_1;
   const NISTMottRS* NMRS70_1;
   const NISTMottRS* NMRS71_1;
   const NISTMottRS* NMRS72_1;
   const NISTMottRS* NMRS73_1;
   const NISTMottRS* NMRS74_1;
   const NISTMottRS* NMRS75_1;
   const NISTMottRS* NMRS76_1;
   const NISTMottRS* NMRS77_1;
   const NISTMottRS* NMRS78_1;
   const NISTMottRS* NMRS79_1;
   const NISTMottRS* NMRS80_1;
   const NISTMottRS* NMRS81_1;
   const NISTMottRS* NMRS82_1;
   const NISTMottRS* NMRS83_1;
   const NISTMottRS* NMRS84_1;
   const NISTMottRS* NMRS85_1;
   const NISTMottRS* NMRS86_1;
   const NISTMottRS* NMRS87_1;
   const NISTMottRS* NMRS88_1;
   const NISTMottRS* NMRS89_1;
   const NISTMottRS* NMRS90_1;
   const NISTMottRS* NMRS91_1;
   const NISTMottRS* NMRS92_1;
   const NISTMottRS* NMRS93_1;
   const NISTMottRS* NMRS94_1;
   const NISTMottRS* NMRS95_1;
   const NISTMottRS* NMRS96_1;

   const NISTMottRS* mScatter1[113] = {
      nullptr,
      NMRS1_1,
      NMRS2_1,
      NMRS3_1,
      NMRS4_1,
      NMRS5_1,
      NMRS6_1,
      NMRS7_1,
      NMRS8_1,
      NMRS9_1,
      NMRS10_1,
      NMRS11_1,
      NMRS12_1,
      NMRS13_1,
      NMRS14_1,
      NMRS15_1,
      NMRS16_1,
      NMRS17_1,
      NMRS18_1,
      NMRS19_1,
      NMRS20_1,
      NMRS21_1,
      NMRS22_1,
      NMRS23_1,
      NMRS24_1,
      NMRS25_1,
      NMRS26_1,
      NMRS27_1,
      NMRS28_1,
      NMRS29_1,
      NMRS30_1,
      NMRS31_1,
      NMRS32_1,
      NMRS33_1,
      NMRS34_1,
      NMRS35_1,
      NMRS36_1,
      NMRS37_1,
      NMRS38_1,
      NMRS39_1,
      NMRS40_1,
      NMRS41_1,
      NMRS42_1,
      NMRS43_1,
      NMRS44_1,
      NMRS45_1,
      NMRS46_1,
      NMRS47_1,
      NMRS48_1,
      NMRS49_1,
      NMRS50_1,
      NMRS51_1,
      NMRS52_1,
      NMRS53_1,
      NMRS54_1,
      NMRS55_1,
      NMRS56_1,
      NMRS57_1,
      NMRS58_1,
      NMRS59_1,
      NMRS60_1,
      NMRS61_1,
      NMRS62_1,
      NMRS63_1,
      NMRS64_1,
      NMRS65_1,
      NMRS66_1,
      NMRS67_1,
      NMRS68_1,
      NMRS69_1,
      NMRS70_1,
      NMRS71_1,
      NMRS72_1,
      NMRS73_1,
      NMRS74_1,
      NMRS75_1,
      NMRS76_1,
      NMRS77_1,
      NMRS78_1,
      NMRS79_1,
      NMRS80_1,
      NMRS81_1,
      NMRS82_1,
      NMRS83_1,
      NMRS84_1,
      NMRS85_1,
      NMRS86_1,
      NMRS87_1,
      NMRS88_1,
      NMRS89_1,
      NMRS90_1,
      NMRS91_1,
      NMRS92_1,
      NMRS93_1,
      NMRS94_1,
      NMRS95_1,
      NMRS96_1
   };

   //const NISTMottRS NMRS1_2(Element::H, 2);
   //const NISTMottRS NMRS2_2(Element::He, 2);
   //const NISTMottRS NMRS3_2(Element::Li, 2);
   //const NISTMottRS NMRS4_2(Element::Be, 2);
   //const NISTMottRS NMRS5_2(Element::B, 2);
   //const NISTMottRS NMRS6_2(Element::C, 2);
   //const NISTMottRS NMRS7_2(Element::N, 2);
   //const NISTMottRS NMRS8_2(Element::O, 2);
   //const NISTMottRS NMRS9_2(Element::F, 2);
   //const NISTMottRS NMRS10_2(Element::Ne, 2);
   //const NISTMottRS NMRS11_2(Element::Na, 2);
   //const NISTMottRS NMRS12_2(Element::Mg, 2);
   //const NISTMottRS NMRS13_2(Element::Al, 2);
   //const NISTMottRS NMRS14_2(Element::Si, 2);
   //const NISTMottRS NMRS15_2(Element::P, 2);
   //const NISTMottRS NMRS16_2(Element::S, 2);
   //const NISTMottRS NMRS17_2(Element::Cl, 2);
   //const NISTMottRS NMRS18_2(Element::Ar, 2);
   //const NISTMottRS NMRS19_2(Element::K, 2);
   //const NISTMottRS NMRS20_2(Element::Ca, 2);
   //const NISTMottRS NMRS21_2(Element::Sc, 2);
   //const NISTMottRS NMRS22_2(Element::Ti, 2);
   //const NISTMottRS NMRS23_2(Element::V, 2);
   //const NISTMottRS NMRS24_2(Element::Cr, 2);
   //const NISTMottRS NMRS25_2(Element::Mn, 2);
   //const NISTMottRS NMRS26_2(Element::Fe, 2);
   //const NISTMottRS NMRS27_2(Element::Co, 2);
   //const NISTMottRS NMRS28_2(Element::Ni, 2);
   //const NISTMottRS NMRS29_2(Element::Cu, 2);
   //const NISTMottRS NMRS30_2(Element::Zn, 2);
   //const NISTMottRS NMRS31_2(Element::Ga, 2);
   //const NISTMottRS NMRS32_2(Element::Ge, 2);
   //const NISTMottRS NMRS33_2(Element::As, 2);
   //const NISTMottRS NMRS34_2(Element::Se, 2);
   //const NISTMottRS NMRS35_2(Element::Br, 2);
   //const NISTMottRS NMRS36_2(Element::Kr, 2);
   //const NISTMottRS NMRS37_2(Element::Rb, 2);
   //const NISTMottRS NMRS38_2(Element::Sr, 2);
   //const NISTMottRS NMRS39_2(Element::Y, 2);
   //const NISTMottRS NMRS40_2(Element::Zr, 2);
   //const NISTMottRS NMRS41_2(Element::Nb, 2);
   //const NISTMottRS NMRS42_2(Element::Mo, 2);
   //const NISTMottRS NMRS43_2(Element::Tc, 2);
   //const NISTMottRS NMRS44_2(Element::Ru, 2);
   //const NISTMottRS NMRS45_2(Element::Rh, 2);
   //const NISTMottRS NMRS46_2(Element::Pd, 2);
   //const NISTMottRS NMRS47_2(Element::Ag, 2);
   //const NISTMottRS NMRS48_2(Element::Cd, 2);
   //const NISTMottRS NMRS49_2(Element::In, 2);
   //const NISTMottRS NMRS50_2(Element::Sn, 2);
   //const NISTMottRS NMRS51_2(Element::Sb, 2);
   //const NISTMottRS NMRS52_2(Element::Te, 2);
   //const NISTMottRS NMRS53_2(Element::I, 2);
   //const NISTMottRS NMRS54_2(Element::Xe, 2);
   //const NISTMottRS NMRS55_2(Element::Cs, 2);
   //const NISTMottRS NMRS56_2(Element::Ba, 2);
   //const NISTMottRS NMRS57_2(Element::La, 2);
   //const NISTMottRS NMRS58_2(Element::Ce, 2);
   //const NISTMottRS NMRS59_2(Element::Pr, 2);
   //const NISTMottRS NMRS60_2(Element::Nd, 2);
   //const NISTMottRS NMRS61_2(Element::Pm, 2);
   //const NISTMottRS NMRS62_2(Element::Sm, 2);
   //const NISTMottRS NMRS63_2(Element::Eu, 2);
   //const NISTMottRS NMRS64_2(Element::Gd, 2);
   //const NISTMottRS NMRS65_2(Element::Tb, 2);
   //const NISTMottRS NMRS66_2(Element::Dy, 2);
   //const NISTMottRS NMRS67_2(Element::Ho, 2);
   //const NISTMottRS NMRS68_2(Element::Er, 2);
   //const NISTMottRS NMRS69_2(Element::Tm, 2);
   //const NISTMottRS NMRS70_2(Element::Yb, 2);
   //const NISTMottRS NMRS71_2(Element::Lu, 2);
   //const NISTMottRS NMRS72_2(Element::Hf, 2);
   //const NISTMottRS NMRS73_2(Element::Ta, 2);
   //const NISTMottRS NMRS74_2(Element::W, 2);
   //const NISTMottRS NMRS75_2(Element::Re, 2);
   //const NISTMottRS NMRS76_2(Element::Os, 2);
   //const NISTMottRS NMRS77_2(Element::Ir, 2);
   //const NISTMottRS NMRS78_2(Element::Pt, 2);
   //const NISTMottRS NMRS79_2(Element::Au, 2);
   //const NISTMottRS NMRS80_2(Element::Hg, 2);
   //const NISTMottRS NMRS81_2(Element::Tl, 2);
   //const NISTMottRS NMRS82_2(Element::Pb, 2);
   //const NISTMottRS NMRS83_2(Element::Bi, 2);
   //const NISTMottRS NMRS84_2(Element::Po, 2);
   //const NISTMottRS NMRS85_2(Element::At, 2);
   //const NISTMottRS NMRS86_2(Element::Rn, 2);
   //const NISTMottRS NMRS87_2(Element::Fr, 2);
   //const NISTMottRS NMRS88_2(Element::Ra, 2);
   //const NISTMottRS NMRS89_2(Element::Ac, 2);
   //const NISTMottRS NMRS90_2(Element::Th, 2);
   //const NISTMottRS NMRS91_2(Element::Pa, 2);
   //const NISTMottRS NMRS92_2(Element::U, 2);
   //const NISTMottRS NMRS93_2(Element::Np, 2);
   //const NISTMottRS NMRS94_2(Element::Pu, 2);
   //const NISTMottRS NMRS95_2(Element::Am, 2);
   //const NISTMottRS NMRS96_2(Element::Cm, 2);

   const NISTMottRS* NMRS1_2;
   const NISTMottRS* NMRS2_2;
   const NISTMottRS* NMRS3_2;
   const NISTMottRS* NMRS4_2;
   const NISTMottRS* NMRS5_2;
   const NISTMottRS* NMRS6_2;
   const NISTMottRS* NMRS7_2;
   const NISTMottRS* NMRS8_2;
   const NISTMottRS* NMRS9_2;
   const NISTMottRS* NMRS10_2;
   const NISTMottRS* NMRS11_2;
   const NISTMottRS* NMRS12_2;
   const NISTMottRS* NMRS13_2;
   const NISTMottRS* NMRS14_2;
   const NISTMottRS* NMRS15_2;
   const NISTMottRS* NMRS16_2;
   const NISTMottRS* NMRS17_2;
   const NISTMottRS* NMRS18_2;
   const NISTMottRS* NMRS19_2;
   const NISTMottRS* NMRS20_2;
   const NISTMottRS* NMRS21_2;
   const NISTMottRS* NMRS22_2;
   const NISTMottRS* NMRS23_2;
   const NISTMottRS* NMRS24_2;
   const NISTMottRS* NMRS25_2;
   const NISTMottRS* NMRS26_2;
   const NISTMottRS* NMRS27_2;
   const NISTMottRS* NMRS28_2;
   const NISTMottRS* NMRS29_2;
   const NISTMottRS* NMRS30_2;
   const NISTMottRS* NMRS31_2;
   const NISTMottRS* NMRS32_2;
   const NISTMottRS* NMRS33_2;
   const NISTMottRS* NMRS34_2;
   const NISTMottRS* NMRS35_2;
   const NISTMottRS* NMRS36_2;
   const NISTMottRS* NMRS37_2;
   const NISTMottRS* NMRS38_2;
   const NISTMottRS* NMRS39_2;
   const NISTMottRS* NMRS40_2;
   const NISTMottRS* NMRS41_2;
   const NISTMottRS* NMRS42_2;
   const NISTMottRS* NMRS43_2;
   const NISTMottRS* NMRS44_2;
   const NISTMottRS* NMRS45_2;
   const NISTMottRS* NMRS46_2;
   const NISTMottRS* NMRS47_2;
   const NISTMottRS* NMRS48_2;
   const NISTMottRS* NMRS49_2;
   const NISTMottRS* NMRS50_2;
   const NISTMottRS* NMRS51_2;
   const NISTMottRS* NMRS52_2;
   const NISTMottRS* NMRS53_2;
   const NISTMottRS* NMRS54_2;
   const NISTMottRS* NMRS55_2;
   const NISTMottRS* NMRS56_2;
   const NISTMottRS* NMRS57_2;
   const NISTMottRS* NMRS58_2;
   const NISTMottRS* NMRS59_2;
   const NISTMottRS* NMRS60_2;
   const NISTMottRS* NMRS61_2;
   const NISTMottRS* NMRS62_2;
   const NISTMottRS* NMRS63_2;
   const NISTMottRS* NMRS64_2;
   const NISTMottRS* NMRS65_2;
   const NISTMottRS* NMRS66_2;
   const NISTMottRS* NMRS67_2;
   const NISTMottRS* NMRS68_2;
   const NISTMottRS* NMRS69_2;
   const NISTMottRS* NMRS70_2;
   const NISTMottRS* NMRS71_2;
   const NISTMottRS* NMRS72_2;
   const NISTMottRS* NMRS73_2;
   const NISTMottRS* NMRS74_2;
   const NISTMottRS* NMRS75_2;
   const NISTMottRS* NMRS76_2;
   const NISTMottRS* NMRS77_2;
   const NISTMottRS* NMRS78_2;
   const NISTMottRS* NMRS79_2;
   const NISTMottRS* NMRS80_2;
   const NISTMottRS* NMRS81_2;
   const NISTMottRS* NMRS82_2;
   const NISTMottRS* NMRS83_2;
   const NISTMottRS* NMRS84_2;
   const NISTMottRS* NMRS85_2;
   const NISTMottRS* NMRS86_2;
   const NISTMottRS* NMRS87_2;
   const NISTMottRS* NMRS88_2;
   const NISTMottRS* NMRS89_2;
   const NISTMottRS* NMRS90_2;
   const NISTMottRS* NMRS91_2;
   const NISTMottRS* NMRS92_2;
   const NISTMottRS* NMRS93_2;
   const NISTMottRS* NMRS94_2;
   const NISTMottRS* NMRS95_2;
   const NISTMottRS* NMRS96_2;

   const NISTMottRS* mScatter2[113] = {
      nullptr,
      NMRS1_2,
      NMRS2_2,
      NMRS3_2,
      NMRS4_2,
      NMRS5_2,
      NMRS6_2,
      NMRS7_2,
      NMRS8_2,
      NMRS9_2,
      NMRS10_2,
      NMRS11_2,
      NMRS12_2,
      NMRS13_2,
      NMRS14_2,
      NMRS15_2,
      NMRS16_2,
      NMRS17_2,
      NMRS18_2,
      NMRS19_2,
      NMRS20_2,
      NMRS21_2,
      NMRS22_2,
      NMRS23_2,
      NMRS24_2,
      NMRS25_2,
      NMRS26_2,
      NMRS27_2,
      NMRS28_2,
      NMRS29_2,
      NMRS30_2,
      NMRS31_2,
      NMRS32_2,
      NMRS33_2,
      NMRS34_2,
      NMRS35_2,
      NMRS36_2,
      NMRS37_2,
      NMRS38_2,
      NMRS39_2,
      NMRS40_2,
      NMRS41_2,
      NMRS42_2,
      NMRS43_2,
      NMRS44_2,
      NMRS45_2,
      NMRS46_2,
      NMRS47_2,
      NMRS48_2,
      NMRS49_2,
      NMRS50_2,
      NMRS51_2,
      NMRS52_2,
      NMRS53_2,
      NMRS54_2,
      NMRS55_2,
      NMRS56_2,
      NMRS57_2,
      NMRS58_2,
      NMRS59_2,
      NMRS60_2,
      NMRS61_2,
      NMRS62_2,
      NMRS63_2,
      NMRS64_2,
      NMRS65_2,
      NMRS66_2,
      NMRS67_2,
      NMRS68_2,
      NMRS69_2,
      NMRS70_2,
      NMRS71_2,
      NMRS72_2,
      NMRS73_2,
      NMRS74_2,
      NMRS75_2,
      NMRS76_2,
      NMRS77_2,
      NMRS78_2,
      NMRS79_2,
      NMRS80_2,
      NMRS81_2,
      NMRS82_2,
      NMRS83_2,
      NMRS84_2,
      NMRS85_2,
      NMRS86_2,
      NMRS87_2,
      NMRS88_2,
      NMRS89_2,
      NMRS90_2,
      NMRS91_2,
      NMRS92_2,
      NMRS93_2,
      NMRS94_2,
      NMRS95_2,
      NMRS96_2
   };

   //const NISTMottRS NMRS1_3(Element::H, 3);
   //const NISTMottRS NMRS2_3(Element::He, 3);
   //const NISTMottRS NMRS3_3(Element::Li, 3);
   //const NISTMottRS NMRS4_3(Element::Be, 3);
   //const NISTMottRS NMRS5_3(Element::B, 3);
   //const NISTMottRS NMRS6_3(Element::C, 3);
   //const NISTMottRS NMRS7_3(Element::N, 3);
   //const NISTMottRS NMRS8_3(Element::O, 3);
   //const NISTMottRS NMRS9_3(Element::F, 3);
   //const NISTMottRS NMRS10_3(Element::Ne, 3);
   //const NISTMottRS NMRS11_3(Element::Na, 3);
   //const NISTMottRS NMRS12_3(Element::Mg, 3);
   //const NISTMottRS NMRS13_3(Element::Al, 3);
   //const NISTMottRS NMRS14_3(Element::Si, 3);
   //const NISTMottRS NMRS15_3(Element::P, 3);
   //const NISTMottRS NMRS16_3(Element::S, 3);
   //const NISTMottRS NMRS17_3(Element::Cl, 3);
   //const NISTMottRS NMRS18_3(Element::Ar, 3);
   //const NISTMottRS NMRS19_3(Element::K, 3);
   //const NISTMottRS NMRS20_3(Element::Ca, 3);
   //const NISTMottRS NMRS21_3(Element::Sc, 3);
   //const NISTMottRS NMRS22_3(Element::Ti, 3);
   //const NISTMottRS NMRS23_3(Element::V, 3);
   //const NISTMottRS NMRS24_3(Element::Cr, 3);
   //const NISTMottRS NMRS25_3(Element::Mn, 3);
   //const NISTMottRS NMRS26_3(Element::Fe, 3);
   //const NISTMottRS NMRS27_3(Element::Co, 3);
   //const NISTMottRS NMRS28_3(Element::Ni, 3);
   //const NISTMottRS NMRS29_3(Element::Cu, 3);
   //const NISTMottRS NMRS30_3(Element::Zn, 3);
   //const NISTMottRS NMRS31_3(Element::Ga, 3);
   //const NISTMottRS NMRS32_3(Element::Ge, 3);
   //const NISTMottRS NMRS33_3(Element::As, 3);
   //const NISTMottRS NMRS34_3(Element::Se, 3);
   //const NISTMottRS NMRS35_3(Element::Br, 3);
   //const NISTMottRS NMRS36_3(Element::Kr, 3);
   //const NISTMottRS NMRS37_3(Element::Rb, 3);
   //const NISTMottRS NMRS38_3(Element::Sr, 3);
   //const NISTMottRS NMRS39_3(Element::Y, 3);
   //const NISTMottRS NMRS40_3(Element::Zr, 3);
   //const NISTMottRS NMRS41_3(Element::Nb, 3);
   //const NISTMottRS NMRS42_3(Element::Mo, 3);
   //const NISTMottRS NMRS43_3(Element::Tc, 3);
   //const NISTMottRS NMRS44_3(Element::Ru, 3);
   //const NISTMottRS NMRS45_3(Element::Rh, 3);
   //const NISTMottRS NMRS46_3(Element::Pd, 3);
   //const NISTMottRS NMRS47_3(Element::Ag, 3);
   //const NISTMottRS NMRS48_3(Element::Cd, 3);
   //const NISTMottRS NMRS49_3(Element::In, 3);
   //const NISTMottRS NMRS50_3(Element::Sn, 3);
   //const NISTMottRS NMRS51_3(Element::Sb, 3);
   //const NISTMottRS NMRS52_3(Element::Te, 3);
   //const NISTMottRS NMRS53_3(Element::I, 3);
   //const NISTMottRS NMRS54_3(Element::Xe, 3);
   //const NISTMottRS NMRS55_3(Element::Cs, 3);
   //const NISTMottRS NMRS56_3(Element::Ba, 3);
   //const NISTMottRS NMRS57_3(Element::La, 3);
   //const NISTMottRS NMRS58_3(Element::Ce, 3);
   //const NISTMottRS NMRS59_3(Element::Pr, 3);
   //const NISTMottRS NMRS60_3(Element::Nd, 3);
   //const NISTMottRS NMRS61_3(Element::Pm, 3);
   //const NISTMottRS NMRS62_3(Element::Sm, 3);
   //const NISTMottRS NMRS63_3(Element::Eu, 3);
   //const NISTMottRS NMRS64_3(Element::Gd, 3);
   //const NISTMottRS NMRS65_3(Element::Tb, 3);
   //const NISTMottRS NMRS66_3(Element::Dy, 3);
   //const NISTMottRS NMRS67_3(Element::Ho, 3);
   //const NISTMottRS NMRS68_3(Element::Er, 3);
   //const NISTMottRS NMRS69_3(Element::Tm, 3);
   //const NISTMottRS NMRS70_3(Element::Yb, 3);
   //const NISTMottRS NMRS71_3(Element::Lu, 3);
   //const NISTMottRS NMRS72_3(Element::Hf, 3);
   //const NISTMottRS NMRS73_3(Element::Ta, 3);
   //const NISTMottRS NMRS74_3(Element::W, 3);
   //const NISTMottRS NMRS75_3(Element::Re, 3);
   //const NISTMottRS NMRS76_3(Element::Os, 3);
   //const NISTMottRS NMRS77_3(Element::Ir, 3);
   //const NISTMottRS NMRS78_3(Element::Pt, 3);
   //const NISTMottRS NMRS79_3(Element::Au, 3);
   //const NISTMottRS NMRS80_3(Element::Hg, 3);
   //const NISTMottRS NMRS81_3(Element::Tl, 3);
   //const NISTMottRS NMRS82_3(Element::Pb, 3);
   //const NISTMottRS NMRS83_3(Element::Bi, 3);
   //const NISTMottRS NMRS84_3(Element::Po, 3);
   //const NISTMottRS NMRS85_3(Element::At, 3);
   //const NISTMottRS NMRS86_3(Element::Rn, 3);
   //const NISTMottRS NMRS87_3(Element::Fr, 3);
   //const NISTMottRS NMRS88_3(Element::Ra, 3);
   //const NISTMottRS NMRS89_3(Element::Ac, 3);
   //const NISTMottRS NMRS90_3(Element::Th, 3);
   //const NISTMottRS NMRS91_3(Element::Pa, 3);
   //const NISTMottRS NMRS92_3(Element::U, 3);
   //const NISTMottRS NMRS93_3(Element::Np, 3);
   //const NISTMottRS NMRS94_3(Element::Pu, 3);
   //const NISTMottRS NMRS95_3(Element::Am, 3);
   //const NISTMottRS NMRS96_3(Element::Cm, 3);

   const NISTMottRS* NMRS1_3;
   const NISTMottRS* NMRS2_3;
   const NISTMottRS* NMRS3_3;
   const NISTMottRS* NMRS4_3;
   const NISTMottRS* NMRS5_3;
   const NISTMottRS* NMRS6_3;
   const NISTMottRS* NMRS7_3;
   const NISTMottRS* NMRS8_3;
   const NISTMottRS* NMRS9_3;
   const NISTMottRS* NMRS10_3;
   const NISTMottRS* NMRS11_3;
   const NISTMottRS* NMRS12_3;
   const NISTMottRS* NMRS13_3;
   const NISTMottRS* NMRS14_3;
   const NISTMottRS* NMRS15_3;
   const NISTMottRS* NMRS16_3;
   const NISTMottRS* NMRS17_3;
   const NISTMottRS* NMRS18_3;
   const NISTMottRS* NMRS19_3;
   const NISTMottRS* NMRS20_3;
   const NISTMottRS* NMRS21_3;
   const NISTMottRS* NMRS22_3;
   const NISTMottRS* NMRS23_3;
   const NISTMottRS* NMRS24_3;
   const NISTMottRS* NMRS25_3;
   const NISTMottRS* NMRS26_3;
   const NISTMottRS* NMRS27_3;
   const NISTMottRS* NMRS28_3;
   const NISTMottRS* NMRS29_3;
   const NISTMottRS* NMRS30_3;
   const NISTMottRS* NMRS31_3;
   const NISTMottRS* NMRS32_3;
   const NISTMottRS* NMRS33_3;
   const NISTMottRS* NMRS34_3;
   const NISTMottRS* NMRS35_3;
   const NISTMottRS* NMRS36_3;
   const NISTMottRS* NMRS37_3;
   const NISTMottRS* NMRS38_3;
   const NISTMottRS* NMRS39_3;
   const NISTMottRS* NMRS40_3;
   const NISTMottRS* NMRS41_3;
   const NISTMottRS* NMRS42_3;
   const NISTMottRS* NMRS43_3;
   const NISTMottRS* NMRS44_3;
   const NISTMottRS* NMRS45_3;
   const NISTMottRS* NMRS46_3;
   const NISTMottRS* NMRS47_3;
   const NISTMottRS* NMRS48_3;
   const NISTMottRS* NMRS49_3;
   const NISTMottRS* NMRS50_3;
   const NISTMottRS* NMRS51_3;
   const NISTMottRS* NMRS52_3;
   const NISTMottRS* NMRS53_3;
   const NISTMottRS* NMRS54_3;
   const NISTMottRS* NMRS55_3;
   const NISTMottRS* NMRS56_3;
   const NISTMottRS* NMRS57_3;
   const NISTMottRS* NMRS58_3;
   const NISTMottRS* NMRS59_3;
   const NISTMottRS* NMRS60_3;
   const NISTMottRS* NMRS61_3;
   const NISTMottRS* NMRS62_3;
   const NISTMottRS* NMRS63_3;
   const NISTMottRS* NMRS64_3;
   const NISTMottRS* NMRS65_3;
   const NISTMottRS* NMRS66_3;
   const NISTMottRS* NMRS67_3;
   const NISTMottRS* NMRS68_3;
   const NISTMottRS* NMRS69_3;
   const NISTMottRS* NMRS70_3;
   const NISTMottRS* NMRS71_3;
   const NISTMottRS* NMRS72_3;
   const NISTMottRS* NMRS73_3;
   const NISTMottRS* NMRS74_3;
   const NISTMottRS* NMRS75_3;
   const NISTMottRS* NMRS76_3;
   const NISTMottRS* NMRS77_3;
   const NISTMottRS* NMRS78_3;
   const NISTMottRS* NMRS79_3;
   const NISTMottRS* NMRS80_3;
   const NISTMottRS* NMRS81_3;
   const NISTMottRS* NMRS82_3;
   const NISTMottRS* NMRS83_3;
   const NISTMottRS* NMRS84_3;
   const NISTMottRS* NMRS85_3;
   const NISTMottRS* NMRS86_3;
   const NISTMottRS* NMRS87_3;
   const NISTMottRS* NMRS88_3;
   const NISTMottRS* NMRS89_3;
   const NISTMottRS* NMRS90_3;
   const NISTMottRS* NMRS91_3;
   const NISTMottRS* NMRS92_3;
   const NISTMottRS* NMRS93_3;
   const NISTMottRS* NMRS94_3;
   const NISTMottRS* NMRS95_3;
   const NISTMottRS* NMRS96_3;

   const NISTMottRS* mScatter3[113] = {
      nullptr,
      NMRS1_3,
      NMRS2_3,
      NMRS3_3,
      NMRS4_3,
      NMRS5_3,
      NMRS6_3,
      NMRS7_3,
      NMRS8_3,
      NMRS9_3,
      NMRS10_3,
      NMRS11_3,
      NMRS12_3,
      NMRS13_3,
      NMRS14_3,
      NMRS15_3,
      NMRS16_3,
      NMRS17_3,
      NMRS18_3,
      NMRS19_3,
      NMRS20_3,
      NMRS21_3,
      NMRS22_3,
      NMRS23_3,
      NMRS24_3,
      NMRS25_3,
      NMRS26_3,
      NMRS27_3,
      NMRS28_3,
      NMRS29_3,
      NMRS30_3,
      NMRS31_3,
      NMRS32_3,
      NMRS33_3,
      NMRS34_3,
      NMRS35_3,
      NMRS36_3,
      NMRS37_3,
      NMRS38_3,
      NMRS39_3,
      NMRS40_3,
      NMRS41_3,
      NMRS42_3,
      NMRS43_3,
      NMRS44_3,
      NMRS45_3,
      NMRS46_3,
      NMRS47_3,
      NMRS48_3,
      NMRS49_3,
      NMRS50_3,
      NMRS51_3,
      NMRS52_3,
      NMRS53_3,
      NMRS54_3,
      NMRS55_3,
      NMRS56_3,
      NMRS57_3,
      NMRS58_3,
      NMRS59_3,
      NMRS60_3,
      NMRS61_3,
      NMRS62_3,
      NMRS63_3,
      NMRS64_3,
      NMRS65_3,
      NMRS66_3,
      NMRS67_3,
      NMRS68_3,
      NMRS69_3,
      NMRS70_3,
      NMRS71_3,
      NMRS72_3,
      NMRS73_3,
      NMRS74_3,
      NMRS75_3,
      NMRS76_3,
      NMRS77_3,
      NMRS78_3,
      NMRS79_3,
      NMRS80_3,
      NMRS81_3,
      NMRS82_3,
      NMRS83_3,
      NMRS84_3,
      NMRS85_3,
      NMRS86_3,
      NMRS87_3,
      NMRS88_3,
      NMRS89_3,
      NMRS90_3,
      NMRS91_3,
      NMRS92_3,
      NMRS93_3,
      NMRS94_3,
      NMRS95_3,
      NMRS96_3
   };

   const NISTMottRS& getNMRS1(int an)
   {
      return *mScatter1[an];
   }

   const NISTMottRS& getNMRS2(int an)
   {
      return *mScatter2[an];
   }

   const NISTMottRS& getNMRS3(int an)
   {
      return *mScatter3[an];
   }

   NISTMottRSFactory::NISTMottRSFactory(int method) : RandomizedScatterFactoryT("NIST Mott Inelastic Cross-Section", mReferenceWebsite), method(method >= 1 && method <= 3 ? method : 1)
   {
   }

   const RandomizedScatterT& NISTMottRSFactory::get(const ElementT& elm) const
   {
      switch (method) {
      case 1:
         return getNMRS1(elm.getAtomicNumber());
      case 2:
         return getNMRS2(elm.getAtomicNumber());
      case 3:
         return getNMRS3(elm.getAtomicNumber());
      default:
         return getNMRS1(elm.getAtomicNumber());
      }
   }

   void NISTMottRSFactory::initializeDefaultStrategy()
   {
   }

   const NISTMottRSFactory FactoryRef = NISTMottRSFactory(1);
   const NISTMottRSFactory Factory100Ref = NISTMottRSFactory(2);
   const NISTMottRSFactory Factory100LinRef = NISTMottRSFactory(3);

   const RandomizedScatterFactoryT& Factory = FactoryRef;
   const RandomizedScatterFactoryT& Factory100 = Factory100Ref;
   const RandomizedScatterFactoryT& Factory100Lin = Factory100LinRef;

   void init()
   {
      NMRS1_1 = new NISTMottRS(Element::H, 1);
      NMRS2_1 = new NISTMottRS(Element::He, 1);
      NMRS3_1 = new NISTMottRS(Element::Li, 1);
      NMRS4_1 = new NISTMottRS(Element::Be, 1);
      NMRS5_1 = new NISTMottRS(Element::B, 1);
      NMRS6_1 = new NISTMottRS(Element::C, 1);
      NMRS7_1 = new NISTMottRS(Element::N, 1);
      NMRS8_1 = new NISTMottRS(Element::O, 1);
      NMRS9_1 = new NISTMottRS(Element::F, 1);
      NMRS10_1 = new NISTMottRS(Element::Ne, 1);
      NMRS11_1 = new NISTMottRS(Element::Na, 1);
      NMRS12_1 = new NISTMottRS(Element::Mg, 1);
      NMRS13_1 = new NISTMottRS(Element::Al, 1);
      NMRS14_1 = new NISTMottRS(Element::Si, 1);
      NMRS15_1 = new NISTMottRS(Element::P, 1);
      NMRS16_1 = new NISTMottRS(Element::S, 1);
      NMRS17_1 = new NISTMottRS(Element::Cl, 1);
      NMRS18_1 = new NISTMottRS(Element::Ar, 1);
      NMRS19_1 = new NISTMottRS(Element::K, 1);
      NMRS20_1 = new NISTMottRS(Element::Ca, 1);
      NMRS21_1 = new NISTMottRS(Element::Sc, 1);
      NMRS22_1 = new NISTMottRS(Element::Ti, 1);
      NMRS23_1 = new NISTMottRS(Element::V, 1);
      NMRS24_1 = new NISTMottRS(Element::Cr, 1);
      NMRS25_1 = new NISTMottRS(Element::Mn, 1);
      NMRS26_1 = new NISTMottRS(Element::Fe, 1);
      NMRS27_1 = new NISTMottRS(Element::Co, 1);
      NMRS28_1 = new NISTMottRS(Element::Ni, 1);
      NMRS29_1 = new NISTMottRS(Element::Cu, 1);
      NMRS30_1 = new NISTMottRS(Element::Zn, 1);
      NMRS31_1 = new NISTMottRS(Element::Ga, 1);
      NMRS32_1 = new NISTMottRS(Element::Ge, 1);
      NMRS33_1 = new NISTMottRS(Element::As, 1);
      NMRS34_1 = new NISTMottRS(Element::Se, 1);
      NMRS35_1 = new NISTMottRS(Element::Br, 1);
      NMRS36_1 = new NISTMottRS(Element::Kr, 1);
      NMRS37_1 = new NISTMottRS(Element::Rb, 1);
      NMRS38_1 = new NISTMottRS(Element::Sr, 1);
      NMRS39_1 = new NISTMottRS(Element::Y, 1);
      NMRS40_1 = new NISTMottRS(Element::Zr, 1);
      NMRS41_1 = new NISTMottRS(Element::Nb, 1);
      NMRS42_1 = new NISTMottRS(Element::Mo, 1);
      NMRS43_1 = new NISTMottRS(Element::Tc, 1);
      NMRS44_1 = new NISTMottRS(Element::Ru, 1);
      NMRS45_1 = new NISTMottRS(Element::Rh, 1);
      NMRS46_1 = new NISTMottRS(Element::Pd, 1);
      NMRS47_1 = new NISTMottRS(Element::Ag, 1);
      NMRS48_1 = new NISTMottRS(Element::Cd, 1);
      NMRS49_1 = new NISTMottRS(Element::In, 1);
      NMRS50_1 = new NISTMottRS(Element::Sn, 1);
      NMRS51_1 = new NISTMottRS(Element::Sb, 1);
      NMRS52_1 = new NISTMottRS(Element::Te, 1);
      NMRS53_1 = new NISTMottRS(Element::I, 1);
      NMRS54_1 = new NISTMottRS(Element::Xe, 1);
      NMRS55_1 = new NISTMottRS(Element::Cs, 1);
      NMRS56_1 = new NISTMottRS(Element::Ba, 1);
      NMRS57_1 = new NISTMottRS(Element::La, 1);
      NMRS58_1 = new NISTMottRS(Element::Ce, 1);
      NMRS59_1 = new NISTMottRS(Element::Pr, 1);
      NMRS60_1 = new NISTMottRS(Element::Nd, 1);
      NMRS61_1 = new NISTMottRS(Element::Pm, 1);
      NMRS62_1 = new NISTMottRS(Element::Sm, 1);
      NMRS63_1 = new NISTMottRS(Element::Eu, 1);
      NMRS64_1 = new NISTMottRS(Element::Gd, 1);
      NMRS65_1 = new NISTMottRS(Element::Tb, 1);
      NMRS66_1 = new NISTMottRS(Element::Dy, 1);
      NMRS67_1 = new NISTMottRS(Element::Ho, 1);
      NMRS68_1 = new NISTMottRS(Element::Er, 1);
      NMRS69_1 = new NISTMottRS(Element::Tm, 1);
      NMRS70_1 = new NISTMottRS(Element::Yb, 1);
      NMRS71_1 = new NISTMottRS(Element::Lu, 1);
      NMRS72_1 = new NISTMottRS(Element::Hf, 1);
      NMRS73_1 = new NISTMottRS(Element::Ta, 1);
      NMRS74_1 = new NISTMottRS(Element::W, 1);
      NMRS75_1 = new NISTMottRS(Element::Re, 1);
      NMRS76_1 = new NISTMottRS(Element::Os, 1);
      NMRS77_1 = new NISTMottRS(Element::Ir, 1);
      NMRS78_1 = new NISTMottRS(Element::Pt, 1);
      NMRS79_1 = new NISTMottRS(Element::Au, 1);
      NMRS80_1 = new NISTMottRS(Element::Hg, 1);
      NMRS81_1 = new NISTMottRS(Element::Tl, 1);
      NMRS82_1 = new NISTMottRS(Element::Pb, 1);
      NMRS83_1 = new NISTMottRS(Element::Bi, 1);
      NMRS84_1 = new NISTMottRS(Element::Po, 1);
      NMRS85_1 = new NISTMottRS(Element::At, 1);
      NMRS86_1 = new NISTMottRS(Element::Rn, 1);
      NMRS87_1 = new NISTMottRS(Element::Fr, 1);
      NMRS88_1 = new NISTMottRS(Element::Ra, 1);
      NMRS89_1 = new NISTMottRS(Element::Ac, 1);
      NMRS90_1 = new NISTMottRS(Element::Th, 1);
      NMRS91_1 = new NISTMottRS(Element::Pa, 1);
      NMRS92_1 = new NISTMottRS(Element::U, 1);
      NMRS93_1 = new NISTMottRS(Element::Np, 1);
      NMRS94_1 = new NISTMottRS(Element::Pu, 1);
      NMRS95_1 = new NISTMottRS(Element::Am, 1);
      NMRS96_1 = new NISTMottRS(Element::Cm, 1);

      NMRS1_2 = new NISTMottRS(Element::H, 2);
      NMRS2_2 = new NISTMottRS(Element::He, 2);
      NMRS3_2 = new NISTMottRS(Element::Li, 2);
      NMRS4_2 = new NISTMottRS(Element::Be, 2);
      NMRS5_2 = new NISTMottRS(Element::B, 2);
      NMRS6_2 = new NISTMottRS(Element::C, 2);
      NMRS7_2 = new NISTMottRS(Element::N, 2);
      NMRS8_2 = new NISTMottRS(Element::O, 2);
      NMRS9_2 = new NISTMottRS(Element::F, 2);
      NMRS10_2 = new NISTMottRS(Element::Ne, 2);
      NMRS11_2 = new NISTMottRS(Element::Na, 2);
      NMRS12_2 = new NISTMottRS(Element::Mg, 2);
      NMRS13_2 = new NISTMottRS(Element::Al, 2);
      NMRS14_2 = new NISTMottRS(Element::Si, 2);
      NMRS15_2 = new NISTMottRS(Element::P, 2);
      NMRS16_2 = new NISTMottRS(Element::S, 2);
      NMRS17_2 = new NISTMottRS(Element::Cl, 2);
      NMRS18_2 = new NISTMottRS(Element::Ar, 2);
      NMRS19_2 = new NISTMottRS(Element::K, 2);
      NMRS20_2 = new NISTMottRS(Element::Ca, 2);
      NMRS21_2 = new NISTMottRS(Element::Sc, 2);
      NMRS22_2 = new NISTMottRS(Element::Ti, 2);
      NMRS23_2 = new NISTMottRS(Element::V, 2);
      NMRS24_2 = new NISTMottRS(Element::Cr, 2);
      NMRS25_2 = new NISTMottRS(Element::Mn, 2);
      NMRS26_2 = new NISTMottRS(Element::Fe, 2);
      NMRS27_2 = new NISTMottRS(Element::Co, 2);
      NMRS28_2 = new NISTMottRS(Element::Ni, 2);
      NMRS29_2 = new NISTMottRS(Element::Cu, 2);
      NMRS30_2 = new NISTMottRS(Element::Zn, 2);
      NMRS31_2 = new NISTMottRS(Element::Ga, 2);
      NMRS32_2 = new NISTMottRS(Element::Ge, 2);
      NMRS33_2 = new NISTMottRS(Element::As, 2);
      NMRS34_2 = new NISTMottRS(Element::Se, 2);
      NMRS35_2 = new NISTMottRS(Element::Br, 2);
      NMRS36_2 = new NISTMottRS(Element::Kr, 2);
      NMRS37_2 = new NISTMottRS(Element::Rb, 2);
      NMRS38_2 = new NISTMottRS(Element::Sr, 2);
      NMRS39_2 = new NISTMottRS(Element::Y, 2);
      NMRS40_2 = new NISTMottRS(Element::Zr, 2);
      NMRS41_2 = new NISTMottRS(Element::Nb, 2);
      NMRS42_2 = new NISTMottRS(Element::Mo, 2);
      NMRS43_2 = new NISTMottRS(Element::Tc, 2);
      NMRS44_2 = new NISTMottRS(Element::Ru, 2);
      NMRS45_2 = new NISTMottRS(Element::Rh, 2);
      NMRS46_2 = new NISTMottRS(Element::Pd, 2);
      NMRS47_2 = new NISTMottRS(Element::Ag, 2);
      NMRS48_2 = new NISTMottRS(Element::Cd, 2);
      NMRS49_2 = new NISTMottRS(Element::In, 2);
      NMRS50_2 = new NISTMottRS(Element::Sn, 2);
      NMRS51_2 = new NISTMottRS(Element::Sb, 2);
      NMRS52_2 = new NISTMottRS(Element::Te, 2);
      NMRS53_2 = new NISTMottRS(Element::I, 2);
      NMRS54_2 = new NISTMottRS(Element::Xe, 2);
      NMRS55_2 = new NISTMottRS(Element::Cs, 2);
      NMRS56_2 = new NISTMottRS(Element::Ba, 2);
      NMRS57_2 = new NISTMottRS(Element::La, 2);
      NMRS58_2 = new NISTMottRS(Element::Ce, 2);
      NMRS59_2 = new NISTMottRS(Element::Pr, 2);
      NMRS60_2 = new NISTMottRS(Element::Nd, 2);
      NMRS61_2 = new NISTMottRS(Element::Pm, 2);
      NMRS62_2 = new NISTMottRS(Element::Sm, 2);
      NMRS63_2 = new NISTMottRS(Element::Eu, 2);
      NMRS64_2 = new NISTMottRS(Element::Gd, 2);
      NMRS65_2 = new NISTMottRS(Element::Tb, 2);
      NMRS66_2 = new NISTMottRS(Element::Dy, 2);
      NMRS67_2 = new NISTMottRS(Element::Ho, 2);
      NMRS68_2 = new NISTMottRS(Element::Er, 2);
      NMRS69_2 = new NISTMottRS(Element::Tm, 2);
      NMRS70_2 = new NISTMottRS(Element::Yb, 2);
      NMRS71_2 = new NISTMottRS(Element::Lu, 2);
      NMRS72_2 = new NISTMottRS(Element::Hf, 2);
      NMRS73_2 = new NISTMottRS(Element::Ta, 2);
      NMRS74_2 = new NISTMottRS(Element::W, 2);
      NMRS75_2 = new NISTMottRS(Element::Re, 2);
      NMRS76_2 = new NISTMottRS(Element::Os, 2);
      NMRS77_2 = new NISTMottRS(Element::Ir, 2);
      NMRS78_2 = new NISTMottRS(Element::Pt, 2);
      NMRS79_2 = new NISTMottRS(Element::Au, 2);
      NMRS80_2 = new NISTMottRS(Element::Hg, 2);
      NMRS81_2 = new NISTMottRS(Element::Tl, 2);
      NMRS82_2 = new NISTMottRS(Element::Pb, 2);
      NMRS83_2 = new NISTMottRS(Element::Bi, 2);
      NMRS84_2 = new NISTMottRS(Element::Po, 2);
      NMRS85_2 = new NISTMottRS(Element::At, 2);
      NMRS86_2 = new NISTMottRS(Element::Rn, 2);
      NMRS87_2 = new NISTMottRS(Element::Fr, 2);
      NMRS88_2 = new NISTMottRS(Element::Ra, 2);
      NMRS89_2 = new NISTMottRS(Element::Ac, 2);
      NMRS90_2 = new NISTMottRS(Element::Th, 2);
      NMRS91_2 = new NISTMottRS(Element::Pa, 2);
      NMRS92_2 = new NISTMottRS(Element::U, 2);
      NMRS93_2 = new NISTMottRS(Element::Np, 2);
      NMRS94_2 = new NISTMottRS(Element::Pu, 2);
      NMRS95_2 = new NISTMottRS(Element::Am, 2);
      NMRS96_2 = new NISTMottRS(Element::Cm, 2);

      NMRS1_3 = new NISTMottRS(Element::H, 3);
      NMRS2_3 = new NISTMottRS(Element::He, 3);
      NMRS3_3 = new NISTMottRS(Element::Li, 3);
      NMRS4_3 = new NISTMottRS(Element::Be, 3);
      NMRS5_3 = new NISTMottRS(Element::B, 3);
      NMRS6_3 = new NISTMottRS(Element::C, 3);
      NMRS7_3 = new NISTMottRS(Element::N, 3);
      NMRS8_3 = new NISTMottRS(Element::O, 3);
      NMRS9_3 = new NISTMottRS(Element::F, 3);
      NMRS10_3 = new NISTMottRS(Element::Ne, 3);
      NMRS11_3 = new NISTMottRS(Element::Na, 3);
      NMRS12_3 = new NISTMottRS(Element::Mg, 3);
      NMRS13_3 = new NISTMottRS(Element::Al, 3);
      NMRS14_3 = new NISTMottRS(Element::Si, 3);
      NMRS15_3 = new NISTMottRS(Element::P, 3);
      NMRS16_3 = new NISTMottRS(Element::S, 3);
      NMRS17_3 = new NISTMottRS(Element::Cl, 3);
      NMRS18_3 = new NISTMottRS(Element::Ar, 3);
      NMRS19_3 = new NISTMottRS(Element::K, 3);
      NMRS20_3 = new NISTMottRS(Element::Ca, 3);
      NMRS21_3 = new NISTMottRS(Element::Sc, 3);
      NMRS22_3 = new NISTMottRS(Element::Ti, 3);
      NMRS23_3 = new NISTMottRS(Element::V, 3);
      NMRS24_3 = new NISTMottRS(Element::Cr, 3);
      NMRS25_3 = new NISTMottRS(Element::Mn, 3);
      NMRS26_3 = new NISTMottRS(Element::Fe, 3);
      NMRS27_3 = new NISTMottRS(Element::Co, 3);
      NMRS28_3 = new NISTMottRS(Element::Ni, 3);
      NMRS29_3 = new NISTMottRS(Element::Cu, 3);
      NMRS30_3 = new NISTMottRS(Element::Zn, 3);
      NMRS31_3 = new NISTMottRS(Element::Ga, 3);
      NMRS32_3 = new NISTMottRS(Element::Ge, 3);
      NMRS33_3 = new NISTMottRS(Element::As, 3);
      NMRS34_3 = new NISTMottRS(Element::Se, 3);
      NMRS35_3 = new NISTMottRS(Element::Br, 3);
      NMRS36_3 = new NISTMottRS(Element::Kr, 3);
      NMRS37_3 = new NISTMottRS(Element::Rb, 3);
      NMRS38_3 = new NISTMottRS(Element::Sr, 3);
      NMRS39_3 = new NISTMottRS(Element::Y, 3);
      NMRS40_3 = new NISTMottRS(Element::Zr, 3);
      NMRS41_3 = new NISTMottRS(Element::Nb, 3);
      NMRS42_3 = new NISTMottRS(Element::Mo, 3);
      NMRS43_3 = new NISTMottRS(Element::Tc, 3);
      NMRS44_3 = new NISTMottRS(Element::Ru, 3);
      NMRS45_3 = new NISTMottRS(Element::Rh, 3);
      NMRS46_3 = new NISTMottRS(Element::Pd, 3);
      NMRS47_3 = new NISTMottRS(Element::Ag, 3);
      NMRS48_3 = new NISTMottRS(Element::Cd, 3);
      NMRS49_3 = new NISTMottRS(Element::In, 3);
      NMRS50_3 = new NISTMottRS(Element::Sn, 3);
      NMRS51_3 = new NISTMottRS(Element::Sb, 3);
      NMRS52_3 = new NISTMottRS(Element::Te, 3);
      NMRS53_3 = new NISTMottRS(Element::I, 3);
      NMRS54_3 = new NISTMottRS(Element::Xe, 3);
      NMRS55_3 = new NISTMottRS(Element::Cs, 3);
      NMRS56_3 = new NISTMottRS(Element::Ba, 3);
      NMRS57_3 = new NISTMottRS(Element::La, 3);
      NMRS58_3 = new NISTMottRS(Element::Ce, 3);
      NMRS59_3 = new NISTMottRS(Element::Pr, 3);
      NMRS60_3 = new NISTMottRS(Element::Nd, 3);
      NMRS61_3 = new NISTMottRS(Element::Pm, 3);
      NMRS62_3 = new NISTMottRS(Element::Sm, 3);
      NMRS63_3 = new NISTMottRS(Element::Eu, 3);
      NMRS64_3 = new NISTMottRS(Element::Gd, 3);
      NMRS65_3 = new NISTMottRS(Element::Tb, 3);
      NMRS66_3 = new NISTMottRS(Element::Dy, 3);
      NMRS67_3 = new NISTMottRS(Element::Ho, 3);
      NMRS68_3 = new NISTMottRS(Element::Er, 3);
      NMRS69_3 = new NISTMottRS(Element::Tm, 3);
      NMRS70_3 = new NISTMottRS(Element::Yb, 3);
      NMRS71_3 = new NISTMottRS(Element::Lu, 3);
      NMRS72_3 = new NISTMottRS(Element::Hf, 3);
      NMRS73_3 = new NISTMottRS(Element::Ta, 3);
      NMRS74_3 = new NISTMottRS(Element::W, 3);
      NMRS75_3 = new NISTMottRS(Element::Re, 3);
      NMRS76_3 = new NISTMottRS(Element::Os, 3);
      NMRS77_3 = new NISTMottRS(Element::Ir, 3);
      NMRS78_3 = new NISTMottRS(Element::Pt, 3);
      NMRS79_3 = new NISTMottRS(Element::Au, 3);
      NMRS80_3 = new NISTMottRS(Element::Hg, 3);
      NMRS81_3 = new NISTMottRS(Element::Tl, 3);
      NMRS82_3 = new NISTMottRS(Element::Pb, 3);
      NMRS83_3 = new NISTMottRS(Element::Bi, 3);
      NMRS84_3 = new NISTMottRS(Element::Po, 3);
      NMRS85_3 = new NISTMottRS(Element::At, 3);
      NMRS86_3 = new NISTMottRS(Element::Rn, 3);
      NMRS87_3 = new NISTMottRS(Element::Fr, 3);
      NMRS88_3 = new NISTMottRS(Element::Ra, 3);
      NMRS89_3 = new NISTMottRS(Element::Ac, 3);
      NMRS90_3 = new NISTMottRS(Element::Th, 3);
      NMRS91_3 = new NISTMottRS(Element::Pa, 3);
      NMRS92_3 = new NISTMottRS(Element::U, 3);
      NMRS93_3 = new NISTMottRS(Element::Np, 3);
      NMRS94_3 = new NISTMottRS(Element::Pu, 3);
      NMRS95_3 = new NISTMottRS(Element::Am, 3);
      NMRS96_3 = new NISTMottRS(Element::Cm, 3);
   }
}