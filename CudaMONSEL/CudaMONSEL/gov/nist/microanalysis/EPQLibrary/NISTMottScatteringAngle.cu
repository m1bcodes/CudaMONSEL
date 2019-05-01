#include "gov\nist\microanalysis\EPQLibrary\NISTMottScatteringAngle.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Reference.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\PhysicalConstants.cuh"
#include "gov\nist\microanalysis\EPQLibrary\FromSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ScreenedRutherfordScatteringAngle.cuh"

namespace NISTMottScatteringAngle
{
   static const Reference::Author auRef[] = { Reference::CPowell, Reference::FSalvat, Reference::AJablonski };
   static const Reference::WebSite REFERENCE("http://www.nist.gov/srd/nist64.htm", "NIST Electron Elastic-Scattering Cross-Section Database version 3.1", "2007 AUGUST 24", auRef, 3);

   const int SPWEM_LEN = 61;
   const int X1_LEN = 201;
   const double DL50 = ::log(50.0);
   const double PARAM = (::log(2.0e4) - DL50) / 60.0;

   const double MAX_NISTMOTT = ToSI::keV(20.0);

   static double value(double a, double b, double c, double y0, double y1, double y2, double x)
   {
      return (x - b) * (x - c) * y0 / ((a - b) * (a - c)) + (x - a) * (x - c) * y1 / ((b - a) * (b - c)) + (x - a) * (x - b)
         * y2 / ((c - a) * (c - b));
   }


   // https://www.oreilly.com/library/view/c-cookbook/0596007612/ch03s06.html
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
         throw (s);
      }

      return (d);
   }

   NISTMottScatteringAngle::NISTMottScatteringAngle(const ElementT& elm) : RandomizedScatterT("NIST Elastic cross-section", REFERENCE), mElement(elm), mSpwem(SPWEM_LEN, 0), mX1(SPWEM_LEN, VectorXd(X1_LEN, 0)), mRutherford(*ScreenedRutherfordScatteringAngle::mScatter[elm.getAtomicNumber()])
   {
      try {
         std::string name = elm.getAtomicNumber() < 10 ? "NistXSec/E0" + std::to_string(elm.getAtomicNumber()) + ".D64" : "NistXSec/E" + std::to_string(elm.getAtomicNumber()) + ".D64";
         std::ifstream t(name);
         std::string line;
         std::getline(t, line);
         for (int j = 0; j < SPWEM_LEN; ++j) {
            std::getline(t, line);
            mSpwem[j] = sciToDub(line);
            for (int i = 0; i < X1_LEN; ++i) {
               std::getline(t, line);
               mX1[j][i] = sciToDub(line);
            }
         }
      }
      catch (std::exception ex) {
         printf("Unable to construct NISTMottScatteringAngle: ");
      }
   }

   StringT NISTMottScatteringAngle::toString() const
   {
      return "CrossSection[NIST-Mott," + StringT(mElement.toAbbrev()) + "]";
   }

   const ElementT& NISTMottScatteringAngle::getElement() const
   {
      return mElement;
   }

   
   double NISTMottScatteringAngle::totalCrossSection(double energy) const
   {
      if (energy < MAX_NISTMOTT) {
         double scale = PhysicalConstants::BohrRadius * PhysicalConstants::BohrRadius;
         double logE = ::log(FromSI::eV(energy));
         int j = 1 + (int)((logE - DL50) / PARAM);
         if (j == 1)
            return value(DL50, DL50 + PARAM, DL50 + 2.0 * PARAM, mSpwem[0], mSpwem[1], mSpwem[2], logE) * scale;
         else if (j == SPWEM_LEN)
            return value(DL50 + 58.0 * PARAM, DL50 + 59.0 * PARAM, DL50 + 60.0 * PARAM, mSpwem[SPWEM_LEN - 3], mSpwem[SPWEM_LEN - 2], mSpwem[SPWEM_LEN - 1], logE)
            * scale;
         else {
            double e0 = DL50 + (j - 2) * PARAM;
            return value(e0, e0 + PARAM, e0 + 2.0 * PARAM, mSpwem[j - 2], mSpwem[j - 1], mSpwem[j], logE) * scale;
         }
      }
      else {
         return mRutherford.totalCrossSection(energy);
      }
   }

   double NISTMottScatteringAngle::randomScatteringAngle(double energy) const
   {
      if (energy < MAX_NISTMOTT) {
         double logE = ::log(FromSI::eV(energy));
         int j = (int)((logE - DL50) / PARAM); // offset to zero-based
         double e2 = DL50 + (j + 1) * PARAM;
         double e1 = e2 - PARAM;
         int i = (logE - e1 < e2 - logE ? j : j + 1); // offset to
         // zero-based
         if (!((i >= 0) && (i < SPWEM_LEN))) printf("%s\n", std::to_string(i) + "\t" + std::to_string(FromSI::eV(energy)) + "\t" + std::to_string(e1) + "\t" + std::to_string(e2));
         // via j
         int k = (int)(200.0 * Math2::random()); // offset to
         // zero-based
         double x = (mX1[i][k + 1] - mX1[i][k]) * Math2::random();
         double q = mX1[i][k] + x;
         double com = 1.0 - 2.0 * q * q;
         return com > -1.0 ? (com < 1.0 ? ::acos(com) : 0.0) : PhysicalConstants::PI;
      }
      else {
         return mRutherford.randomScatteringAngle(energy);
      }
   }

   const NISTMottScatteringAngle NISTMSA1(Element::H);
   const NISTMottScatteringAngle NISTMSA2(Element::He);
   const NISTMottScatteringAngle NISTMSA3(Element::Li);
   const NISTMottScatteringAngle NISTMSA4(Element::Be);
   const NISTMottScatteringAngle NISTMSA5(Element::B);
   const NISTMottScatteringAngle NISTMSA6(Element::C);
   const NISTMottScatteringAngle NISTMSA7(Element::N);
   const NISTMottScatteringAngle NISTMSA8(Element::O);
   const NISTMottScatteringAngle NISTMSA9(Element::F);
   const NISTMottScatteringAngle NISTMSA10(Element::Ne);
   const NISTMottScatteringAngle NISTMSA11(Element::Na);
   const NISTMottScatteringAngle NISTMSA12(Element::Mg);
   const NISTMottScatteringAngle NISTMSA13(Element::Al);
   const NISTMottScatteringAngle NISTMSA14(Element::Si);
   const NISTMottScatteringAngle NISTMSA15(Element::P);
   const NISTMottScatteringAngle NISTMSA16(Element::S);
   const NISTMottScatteringAngle NISTMSA17(Element::Cl);
   const NISTMottScatteringAngle NISTMSA18(Element::Ar);
   const NISTMottScatteringAngle NISTMSA19(Element::K);
   const NISTMottScatteringAngle NISTMSA20(Element::Ca);
   const NISTMottScatteringAngle NISTMSA21(Element::Sc);
   const NISTMottScatteringAngle NISTMSA22(Element::Ti);
   const NISTMottScatteringAngle NISTMSA23(Element::V);
   const NISTMottScatteringAngle NISTMSA24(Element::Cr);
   const NISTMottScatteringAngle NISTMSA25(Element::Mn);
   const NISTMottScatteringAngle NISTMSA26(Element::Fe);
   const NISTMottScatteringAngle NISTMSA27(Element::Co);
   const NISTMottScatteringAngle NISTMSA28(Element::Ni);
   const NISTMottScatteringAngle NISTMSA29(Element::Cu);
   const NISTMottScatteringAngle NISTMSA30(Element::Zn);
   const NISTMottScatteringAngle NISTMSA31(Element::Ga);
   const NISTMottScatteringAngle NISTMSA32(Element::Ge);
   const NISTMottScatteringAngle NISTMSA33(Element::As);
   const NISTMottScatteringAngle NISTMSA34(Element::Se);
   const NISTMottScatteringAngle NISTMSA35(Element::Br);
   const NISTMottScatteringAngle NISTMSA36(Element::Kr);
   const NISTMottScatteringAngle NISTMSA37(Element::Rb);
   const NISTMottScatteringAngle NISTMSA38(Element::Sr);
   const NISTMottScatteringAngle NISTMSA39(Element::Y);
   const NISTMottScatteringAngle NISTMSA40(Element::Zr);
   const NISTMottScatteringAngle NISTMSA41(Element::Nb);
   const NISTMottScatteringAngle NISTMSA42(Element::Mo);
   const NISTMottScatteringAngle NISTMSA43(Element::Tc);
   const NISTMottScatteringAngle NISTMSA44(Element::Ru);
   const NISTMottScatteringAngle NISTMSA45(Element::Rh);
   const NISTMottScatteringAngle NISTMSA46(Element::Pd);
   const NISTMottScatteringAngle NISTMSA47(Element::Ag);
   const NISTMottScatteringAngle NISTMSA48(Element::Cd);
   const NISTMottScatteringAngle NISTMSA49(Element::In);
   const NISTMottScatteringAngle NISTMSA50(Element::Sn);
   const NISTMottScatteringAngle NISTMSA51(Element::Sb);
   const NISTMottScatteringAngle NISTMSA52(Element::Te);
   const NISTMottScatteringAngle NISTMSA53(Element::I);
   const NISTMottScatteringAngle NISTMSA54(Element::Xe);
   const NISTMottScatteringAngle NISTMSA55(Element::Cs);
   const NISTMottScatteringAngle NISTMSA56(Element::Ba);
   const NISTMottScatteringAngle NISTMSA57(Element::La);
   const NISTMottScatteringAngle NISTMSA58(Element::Ce);
   const NISTMottScatteringAngle NISTMSA59(Element::Pr);
   const NISTMottScatteringAngle NISTMSA60(Element::Nd);
   const NISTMottScatteringAngle NISTMSA61(Element::Pm);
   const NISTMottScatteringAngle NISTMSA62(Element::Sm);
   const NISTMottScatteringAngle NISTMSA63(Element::Eu);
   const NISTMottScatteringAngle NISTMSA64(Element::Gd);
   const NISTMottScatteringAngle NISTMSA65(Element::Tb);
   const NISTMottScatteringAngle NISTMSA66(Element::Dy);
   const NISTMottScatteringAngle NISTMSA67(Element::Ho);
   const NISTMottScatteringAngle NISTMSA68(Element::Er);
   const NISTMottScatteringAngle NISTMSA69(Element::Tm);
   const NISTMottScatteringAngle NISTMSA70(Element::Yb);
   const NISTMottScatteringAngle NISTMSA71(Element::Lu);
   const NISTMottScatteringAngle NISTMSA72(Element::Hf);
   const NISTMottScatteringAngle NISTMSA73(Element::Ta);
   const NISTMottScatteringAngle NISTMSA74(Element::W);
   const NISTMottScatteringAngle NISTMSA75(Element::Re);
   const NISTMottScatteringAngle NISTMSA76(Element::Os);
   const NISTMottScatteringAngle NISTMSA77(Element::Ir);
   const NISTMottScatteringAngle NISTMSA78(Element::Pt);
   const NISTMottScatteringAngle NISTMSA79(Element::Au);
   const NISTMottScatteringAngle NISTMSA80(Element::Hg);
   const NISTMottScatteringAngle NISTMSA81(Element::Tl);
   const NISTMottScatteringAngle NISTMSA82(Element::Pb);
   const NISTMottScatteringAngle NISTMSA83(Element::Bi);
   const NISTMottScatteringAngle NISTMSA84(Element::Po);
   const NISTMottScatteringAngle NISTMSA85(Element::At);
   const NISTMottScatteringAngle NISTMSA86(Element::Rn);
   const NISTMottScatteringAngle NISTMSA87(Element::Fr);
   const NISTMottScatteringAngle NISTMSA88(Element::Ra);
   const NISTMottScatteringAngle NISTMSA89(Element::Ac);
   const NISTMottScatteringAngle NISTMSA90(Element::Th);
   const NISTMottScatteringAngle NISTMSA91(Element::Pa);
   const NISTMottScatteringAngle NISTMSA92(Element::U);
   const NISTMottScatteringAngle NISTMSA93(Element::Np);
   const NISTMottScatteringAngle NISTMSA94(Element::Pu);

   const NISTMottScatteringAngle* mScatter[113] = {
      &NISTMSA1,
      &NISTMSA2,
      &NISTMSA3,
      &NISTMSA4,
      &NISTMSA5,
      &NISTMSA6,
      &NISTMSA7,
      &NISTMSA8,
      &NISTMSA9,
      &NISTMSA10,
      &NISTMSA11,
      &NISTMSA12,
      &NISTMSA13,
      &NISTMSA14,
      &NISTMSA15,
      &NISTMSA16,
      &NISTMSA17,
      &NISTMSA18,
      &NISTMSA19,
      &NISTMSA20,
      &NISTMSA21,
      &NISTMSA22,
      &NISTMSA23,
      &NISTMSA24,
      &NISTMSA25,
      &NISTMSA26,
      &NISTMSA27,
      &NISTMSA28,
      &NISTMSA29,
      &NISTMSA30,
      &NISTMSA31,
      &NISTMSA32,
      &NISTMSA33,
      &NISTMSA34,
      &NISTMSA35,
      &NISTMSA36,
      &NISTMSA37,
      &NISTMSA38,
      &NISTMSA39,
      &NISTMSA40,
      &NISTMSA41,
      &NISTMSA42,
      &NISTMSA43,
      &NISTMSA44,
      &NISTMSA45,
      &NISTMSA46,
      &NISTMSA47,
      &NISTMSA48,
      &NISTMSA49,
      &NISTMSA50,
      &NISTMSA51,
      &NISTMSA52,
      &NISTMSA53,
      &NISTMSA54,
      &NISTMSA55,
      &NISTMSA56,
      &NISTMSA57,
      &NISTMSA58,
      &NISTMSA59,
      &NISTMSA60,
      &NISTMSA61,
      &NISTMSA62,
      &NISTMSA63,
      &NISTMSA64,
      &NISTMSA65,
      &NISTMSA66,
      &NISTMSA67,
      &NISTMSA68,
      &NISTMSA69,
      &NISTMSA70,
      &NISTMSA71,
      &NISTMSA72,
      &NISTMSA73,
      &NISTMSA74,
      &NISTMSA75,
      &NISTMSA76,
      &NISTMSA77,
      &NISTMSA78,
      &NISTMSA79,
      &NISTMSA80,
      &NISTMSA81,
      &NISTMSA82,
      &NISTMSA83,
      &NISTMSA84,
      &NISTMSA85,
      &NISTMSA86,
      &NISTMSA87,
      &NISTMSA88,
      &NISTMSA89,
      &NISTMSA90,
      &NISTMSA91,
      &NISTMSA92,
      &NISTMSA93,
      &NISTMSA94
   };

   NISTMottRandomizedScatterFactory::NISTMottRandomizedScatterFactory() : RandomizedScatterFactoryT("NIST Mott Inelastic Cross-Section", REFERENCE)
   {
   }

   void NISTMottRandomizedScatterFactory::initializeDefaultStrategy()
   {
   }

   const RandomizedScatterT& NISTMottRandomizedScatterFactory::get(const ElementT& elm) const
   {
      return *mScatter[elm.getAtomicNumber()];
   }

   const NISTMottRandomizedScatterFactory NISTMottRandomizedFactory;
   const RandomizedScatterFactoryT& Factory = NISTMottRandomizedFactory;
}
