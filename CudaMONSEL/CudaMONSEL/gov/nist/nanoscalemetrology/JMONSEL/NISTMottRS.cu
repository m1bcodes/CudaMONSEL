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

#include "Amphibian\random.cuh"

namespace NISTMottRS
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ static const float MAX_NISTMOTT = 3.2043531e-15f;
   __constant__ static const float MIN_NISTMOTT = 8.0108827e-18f;

   __constant__ static const int qINTERPOLATIONORDER = 3;
   __constant__ static const int sigmaINTERPOLATIONORDER = 3;
   __constant__ static const float scale = 2.8002852e-21f;

   __constant__ static const int SPWEM_LEN = 61;
   __constant__ static const int X1_LEN = 201;
   __constant__ static const float DL50 = -39.3657307254f;
   __constant__ static const float PARAM = 0.09985774245f;
#else
   static const float MAX_NISTMOTT = ToSI::keV(20.0);
   static const float MIN_NISTMOTT = ToSI::keV(0.050);

   static const int qINTERPOLATIONORDER = 3;
   static const int sigmaINTERPOLATIONORDER = 3;
   static const float scale = PhysicalConstants::BohrRadius * PhysicalConstants::BohrRadius;

   static const int SPWEM_LEN = 61;
   static const int X1_LEN = 201;
   static const float DL50 = ::log(MIN_NISTMOTT);
   static const float PARAM = (::log(MAX_NISTMOTT) - DL50) / 60.0;
#endif

   static const Reference::Author* al[] = {
      &Reference::FSalvat,
      &Reference::AJablonski,
      &Reference::CPowell
   };
   static const Reference::WebSite mReferenceWebsite("http://www.nist.gov/srd/nist64.htm", "NIST Electron Elastic-Scattering Cross-Section Database version 3.1", "AUGUST 24, 2007", al, 3);

   //static double sciToDub(const std::string& str)
   //{
   //   std::string tmp = str.substr(str.find_first_not_of(" "));
   //   std::stringstream ss(tmp);
   //   double d = 0;
   //   ss >> d;

   //   if (ss.fail()) {
   //      std::string s = "Unable to format ";
   //      s += tmp;
   //      s += " as a number!";
   //      throw s;
   //   }

   //   return d;
   //}

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

   __host__ __device__ NISTMottRS::NISTMottRS(const ElementT& elm, int method) :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      RandomizedScatterT("NIST-Mott Elastic cross-section", *Reference::d_NullReference),
#else
      RandomizedScatterT("NIST-Mott Elastic cross-section", mReferenceWebsite),
#endif
      mElement(elm),
      method(method >= 1 && method <= 3 ? method : 1),
      mRutherford(ScreenedRutherfordScatteringAngle::getSRSA(elm.getAtomicNumber())),
      mBrowning(BrowningEmpiricalCrossSection::getBECS(elm.getAtomicNumber())),
      extrapolateBelowEnergy(method == 1 ? ToSI::eV(50.) : ToSI::eV(100.)),
      mSpwem(NISTMottScatteringAngle::getNISTMSA(elm.getAtomicNumber()).getSpwem()),
      mX1(NISTMottScatteringAngle::getNISTMSA(elm.getAtomicNumber()).getX1()),
      MottXSatMinEnergy(totalCrossSection(extrapolateBelowEnergy)),
      sfBrowning(MottXSatMinEnergy / mBrowning.totalCrossSection(extrapolateBelowEnergy)),
      name(StringT("CrossSection[NIST-Mott, ") + StringT(mElement.toAbbrev()) + "]")
   {
      //loadData(elm.getAtomicNumber());
      //if (!(method >= 1 && method <= 3))
      //   printf("NISTMottRS::setMethod: Invalid NISTMottRS method: method must = 1, 2, or 3.");
   }

   const char* NISTMottRS::toString()
   {
      return name.c_str();
   }

   __host__ __device__ const ElementT& NISTMottRS::getElement() const
   {
      return mElement;
   }

   __host__ __device__ double NISTMottRS::totalCrossSection(const double energy) const
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
         return scale * ULagrangeInterpolation::d1(mSpwem, DL50, PARAM, sigmaINTERPOLATIONORDER, ::logf(energy))[0];
      }
      else {
         return mRutherford.totalCrossSection(energy);
      }
   }

   __host__ __device__ double NISTMottRS::randomScatteringAngle(const double energy) const
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
         const float x0[] =  {
            DL50,
               0.
         };
         const float xinc[] = {
            PARAM,
            0.005
         };
         const float x[] = {
            ::log(energy),
            Random::random()
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

   const NISTMottRS* mScatter1[113];
   const NISTMottRS* mScatter2[113];
   const NISTMottRS* mScatter3[113];

   __device__ const NISTMottRS* d_mScatter1[113];
   __device__ const NISTMottRS* d_mScatter2[113];
   __device__ const NISTMottRS* d_mScatter3[113];

   __host__ __device__ const NISTMottRS& getNMRS1(int an)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      return *d_mScatter1[an];
#else
      return *mScatter1[an];
#endif
   }

   __host__ __device__ const NISTMottRS& getNMRS2(int an)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      return *d_mScatter2[an];
#else
      return *mScatter2[an];
#endif
   }

   __host__ __device__ const NISTMottRS& getNMRS3(int an)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      return *d_mScatter3[an];
#else
      return *mScatter3[an];
#endif
   }

   __host__ __device__ NISTMottRSFactory::NISTMottRSFactory(int method) :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      RandomizedScatterFactoryT("NIST Mott Inelastic Cross-Section", *Reference::d_NullReference), method(method >= 1 && method <= 3 ? method : 1)
#else
      RandomizedScatterFactoryT("NIST Mott Inelastic Cross-Section", mReferenceWebsite), method(method >= 1 && method <= 3 ? method : 1)
#endif
   {
   }

   __host__ __device__ const RandomizedScatterT& NISTMottRSFactory::get(const ElementT& elm) const
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

   const NISTMottRSFactory FactoryRef(1);
   const NISTMottRSFactory Factory100Ref(2);
   const NISTMottRSFactory Factory100LinRef(3);

   const RandomizedScatterFactoryT& Factory = FactoryRef;
   const RandomizedScatterFactoryT& Factory100 = Factory100Ref;
   const RandomizedScatterFactoryT& Factory100Lin = Factory100LinRef;

   __device__ const RandomizedScatterFactoryT* dFactory = nullptr;
   __device__ const RandomizedScatterFactoryT* dFactory100 = nullptr;
   __device__ const RandomizedScatterFactoryT* dFactory100Lin = nullptr;

   __global__ void initFactory()
   {
      dFactory = new NISTMottRSFactory(1);
      dFactory100 = new NISTMottRSFactory(2);
      dFactory100Lin = new NISTMottRSFactory(3);
   }

   void init()
   {
      mScatter1[1] = new NISTMottRS(Element::H, 1);
      mScatter1[2] = new NISTMottRS(Element::He, 1);
      mScatter1[3] = new NISTMottRS(Element::Li, 1);
      mScatter1[4] = new NISTMottRS(Element::Be, 1);
      mScatter1[5] = new NISTMottRS(Element::B, 1);
      mScatter1[6] = new NISTMottRS(Element::C, 1);
      mScatter1[7] = new NISTMottRS(Element::N, 1);
      mScatter1[8] = new NISTMottRS(Element::O, 1);
      mScatter1[9] = new NISTMottRS(Element::F, 1);
      mScatter1[10] = new NISTMottRS(Element::Ne, 1);
      mScatter1[11] = new NISTMottRS(Element::Na, 1);
      mScatter1[12] = new NISTMottRS(Element::Mg, 1);
      mScatter1[13] = new NISTMottRS(Element::Al, 1);
      mScatter1[14] = new NISTMottRS(Element::Si, 1);
      mScatter1[15] = new NISTMottRS(Element::P, 1);
      mScatter1[16] = new NISTMottRS(Element::S, 1);
      mScatter1[17] = new NISTMottRS(Element::Cl, 1);
      mScatter1[18] = new NISTMottRS(Element::Ar, 1);
      mScatter1[19] = new NISTMottRS(Element::K, 1);
      mScatter1[20] = new NISTMottRS(Element::Ca, 1);
      mScatter1[21] = new NISTMottRS(Element::Sc, 1);
      mScatter1[22] = new NISTMottRS(Element::Ti, 1);
      mScatter1[23] = new NISTMottRS(Element::V, 1);
      mScatter1[24] = new NISTMottRS(Element::Cr, 1);
      mScatter1[25] = new NISTMottRS(Element::Mn, 1);
      mScatter1[26] = new NISTMottRS(Element::Fe, 1);
      mScatter1[27] = new NISTMottRS(Element::Co, 1);
      mScatter1[28] = new NISTMottRS(Element::Ni, 1);
      mScatter1[29] = new NISTMottRS(Element::Cu, 1);
      mScatter1[30] = new NISTMottRS(Element::Zn, 1);
      mScatter1[31] = new NISTMottRS(Element::Ga, 1);
      mScatter1[32] = new NISTMottRS(Element::Ge, 1);
      mScatter1[33] = new NISTMottRS(Element::As, 1);
      mScatter1[34] = new NISTMottRS(Element::Se, 1);
      mScatter1[35] = new NISTMottRS(Element::Br, 1);
      mScatter1[36] = new NISTMottRS(Element::Kr, 1);
      mScatter1[37] = new NISTMottRS(Element::Rb, 1);
      mScatter1[38] = new NISTMottRS(Element::Sr, 1);
      mScatter1[39] = new NISTMottRS(Element::Y, 1);
      mScatter1[40] = new NISTMottRS(Element::Zr, 1);
      mScatter1[41] = new NISTMottRS(Element::Nb, 1);
      mScatter1[42] = new NISTMottRS(Element::Mo, 1);
      mScatter1[43] = new NISTMottRS(Element::Tc, 1);
      mScatter1[44] = new NISTMottRS(Element::Ru, 1);
      mScatter1[45] = new NISTMottRS(Element::Rh, 1);
      mScatter1[46] = new NISTMottRS(Element::Pd, 1);
      mScatter1[47] = new NISTMottRS(Element::Ag, 1);
      mScatter1[48] = new NISTMottRS(Element::Cd, 1);
      mScatter1[49] = new NISTMottRS(Element::In, 1);
      mScatter1[50] = new NISTMottRS(Element::Sn, 1);
      mScatter1[51] = new NISTMottRS(Element::Sb, 1);
      mScatter1[52] = new NISTMottRS(Element::Te, 1);
      mScatter1[53] = new NISTMottRS(Element::I, 1);
      mScatter1[54] = new NISTMottRS(Element::Xe, 1);
      mScatter1[55] = new NISTMottRS(Element::Cs, 1);
      mScatter1[56] = new NISTMottRS(Element::Ba, 1);
      mScatter1[57] = new NISTMottRS(Element::La, 1);
      mScatter1[58] = new NISTMottRS(Element::Ce, 1);
      mScatter1[59] = new NISTMottRS(Element::Pr, 1);
      mScatter1[60] = new NISTMottRS(Element::Nd, 1);
      mScatter1[61] = new NISTMottRS(Element::Pm, 1);
      mScatter1[62] = new NISTMottRS(Element::Sm, 1);
      mScatter1[63] = new NISTMottRS(Element::Eu, 1);
      mScatter1[64] = new NISTMottRS(Element::Gd, 1);
      mScatter1[65] = new NISTMottRS(Element::Tb, 1);
      mScatter1[66] = new NISTMottRS(Element::Dy, 1);
      mScatter1[67] = new NISTMottRS(Element::Ho, 1);
      mScatter1[68] = new NISTMottRS(Element::Er, 1);
      mScatter1[69] = new NISTMottRS(Element::Tm, 1);
      mScatter1[70] = new NISTMottRS(Element::Yb, 1);
      mScatter1[71] = new NISTMottRS(Element::Lu, 1);
      mScatter1[72] = new NISTMottRS(Element::Hf, 1);
      mScatter1[73] = new NISTMottRS(Element::Ta, 1);
      mScatter1[74] = new NISTMottRS(Element::W, 1);
      mScatter1[75] = new NISTMottRS(Element::Re, 1);
      mScatter1[76] = new NISTMottRS(Element::Os, 1);
      mScatter1[77] = new NISTMottRS(Element::Ir, 1);
      mScatter1[78] = new NISTMottRS(Element::Pt, 1);
      mScatter1[79] = new NISTMottRS(Element::Au, 1);
      mScatter1[80] = new NISTMottRS(Element::Hg, 1);
      mScatter1[81] = new NISTMottRS(Element::Tl, 1);
      mScatter1[82] = new NISTMottRS(Element::Pb, 1);
      mScatter1[83] = new NISTMottRS(Element::Bi, 1);
      mScatter1[84] = new NISTMottRS(Element::Po, 1);
      mScatter1[85] = new NISTMottRS(Element::At, 1);
      mScatter1[86] = new NISTMottRS(Element::Rn, 1);
      mScatter1[87] = new NISTMottRS(Element::Fr, 1);
      mScatter1[88] = new NISTMottRS(Element::Ra, 1);
      mScatter1[89] = new NISTMottRS(Element::Ac, 1);
      mScatter1[90] = new NISTMottRS(Element::Th, 1);
      mScatter1[91] = new NISTMottRS(Element::Pa, 1);
      mScatter1[92] = new NISTMottRS(Element::U, 1);
      mScatter1[93] = new NISTMottRS(Element::Np, 1);
      mScatter1[94] = new NISTMottRS(Element::Pu, 1);
      mScatter1[95] = new NISTMottRS(Element::Am, 1);
      mScatter1[96] = new NISTMottRS(Element::Cm, 1);

      mScatter2[1] = new NISTMottRS(Element::H, 2);
      mScatter2[2] = new NISTMottRS(Element::He, 2);
      mScatter2[3] = new NISTMottRS(Element::Li, 2);
      mScatter2[4] = new NISTMottRS(Element::Be, 2);
      mScatter2[5] = new NISTMottRS(Element::B, 2);
      mScatter2[6] = new NISTMottRS(Element::C, 2);
      mScatter2[7] = new NISTMottRS(Element::N, 2);
      mScatter2[8] = new NISTMottRS(Element::O, 2);
      mScatter2[9] = new NISTMottRS(Element::F, 2);
      mScatter2[10] = new NISTMottRS(Element::Ne, 2);
      mScatter2[11] = new NISTMottRS(Element::Na, 2);
      mScatter2[12] = new NISTMottRS(Element::Mg, 2);
      mScatter2[13] = new NISTMottRS(Element::Al, 2);
      mScatter2[14] = new NISTMottRS(Element::Si, 2);
      mScatter2[15] = new NISTMottRS(Element::P, 2);
      mScatter2[16] = new NISTMottRS(Element::S, 2);
      mScatter2[17] = new NISTMottRS(Element::Cl, 2);
      mScatter2[18] = new NISTMottRS(Element::Ar, 2);
      mScatter2[19] = new NISTMottRS(Element::K, 2);
      mScatter2[20] = new NISTMottRS(Element::Ca, 2);
      mScatter2[21] = new NISTMottRS(Element::Sc, 2);
      mScatter2[22] = new NISTMottRS(Element::Ti, 2);
      mScatter2[23] = new NISTMottRS(Element::V, 2);
      mScatter2[24] = new NISTMottRS(Element::Cr, 2);
      mScatter2[25] = new NISTMottRS(Element::Mn, 2);
      mScatter2[26] = new NISTMottRS(Element::Fe, 2);
      mScatter2[27] = new NISTMottRS(Element::Co, 2);
      mScatter2[28] = new NISTMottRS(Element::Ni, 2);
      mScatter2[29] = new NISTMottRS(Element::Cu, 2);
      mScatter2[30] = new NISTMottRS(Element::Zn, 2);
      mScatter2[31] = new NISTMottRS(Element::Ga, 2);
      mScatter2[32] = new NISTMottRS(Element::Ge, 2);
      mScatter2[33] = new NISTMottRS(Element::As, 2);
      mScatter2[34] = new NISTMottRS(Element::Se, 2);
      mScatter2[35] = new NISTMottRS(Element::Br, 2);
      mScatter2[36] = new NISTMottRS(Element::Kr, 2);
      mScatter2[37] = new NISTMottRS(Element::Rb, 2);
      mScatter2[38] = new NISTMottRS(Element::Sr, 2);
      mScatter2[39] = new NISTMottRS(Element::Y, 2);
      mScatter2[40] = new NISTMottRS(Element::Zr, 2);
      mScatter2[41] = new NISTMottRS(Element::Nb, 2);
      mScatter2[42] = new NISTMottRS(Element::Mo, 2);
      mScatter2[43] = new NISTMottRS(Element::Tc, 2);
      mScatter2[44] = new NISTMottRS(Element::Ru, 2);
      mScatter2[45] = new NISTMottRS(Element::Rh, 2);
      mScatter2[46] = new NISTMottRS(Element::Pd, 2);
      mScatter2[47] = new NISTMottRS(Element::Ag, 2);
      mScatter2[48] = new NISTMottRS(Element::Cd, 2);
      mScatter2[49] = new NISTMottRS(Element::In, 2);
      mScatter2[50] = new NISTMottRS(Element::Sn, 2);
      mScatter2[51] = new NISTMottRS(Element::Sb, 2);
      mScatter2[52] = new NISTMottRS(Element::Te, 2);
      mScatter2[53] = new NISTMottRS(Element::I, 2);
      mScatter2[54] = new NISTMottRS(Element::Xe, 2);
      mScatter2[55] = new NISTMottRS(Element::Cs, 2);
      mScatter2[56] = new NISTMottRS(Element::Ba, 2);
      mScatter2[57] = new NISTMottRS(Element::La, 2);
      mScatter2[58] = new NISTMottRS(Element::Ce, 2);
      mScatter2[59] = new NISTMottRS(Element::Pr, 2);
      mScatter2[60] = new NISTMottRS(Element::Nd, 2);
      mScatter2[61] = new NISTMottRS(Element::Pm, 2);
      mScatter2[62] = new NISTMottRS(Element::Sm, 2);
      mScatter2[63] = new NISTMottRS(Element::Eu, 2);
      mScatter2[64] = new NISTMottRS(Element::Gd, 2);
      mScatter2[65] = new NISTMottRS(Element::Tb, 2);
      mScatter2[66] = new NISTMottRS(Element::Dy, 2);
      mScatter2[67] = new NISTMottRS(Element::Ho, 2);
      mScatter2[68] = new NISTMottRS(Element::Er, 2);
      mScatter2[69] = new NISTMottRS(Element::Tm, 2);
      mScatter2[70] = new NISTMottRS(Element::Yb, 2);
      mScatter2[71] = new NISTMottRS(Element::Lu, 2);
      mScatter2[72] = new NISTMottRS(Element::Hf, 2);
      mScatter2[73] = new NISTMottRS(Element::Ta, 2);
      mScatter2[74] = new NISTMottRS(Element::W, 2);
      mScatter2[75] = new NISTMottRS(Element::Re, 2);
      mScatter2[76] = new NISTMottRS(Element::Os, 2);
      mScatter2[77] = new NISTMottRS(Element::Ir, 2);
      mScatter2[78] = new NISTMottRS(Element::Pt, 2);
      mScatter2[79] = new NISTMottRS(Element::Au, 2);
      mScatter2[80] = new NISTMottRS(Element::Hg, 2);
      mScatter2[81] = new NISTMottRS(Element::Tl, 2);
      mScatter2[82] = new NISTMottRS(Element::Pb, 2);
      mScatter2[83] = new NISTMottRS(Element::Bi, 2);
      mScatter2[84] = new NISTMottRS(Element::Po, 2);
      mScatter2[85] = new NISTMottRS(Element::At, 2);
      mScatter2[86] = new NISTMottRS(Element::Rn, 2);
      mScatter2[87] = new NISTMottRS(Element::Fr, 2);
      mScatter2[88] = new NISTMottRS(Element::Ra, 2);
      mScatter2[89] = new NISTMottRS(Element::Ac, 2);
      mScatter2[90] = new NISTMottRS(Element::Th, 2);
      mScatter2[91] = new NISTMottRS(Element::Pa, 2);
      mScatter2[92] = new NISTMottRS(Element::U, 2);
      mScatter2[93] = new NISTMottRS(Element::Np, 2);
      mScatter2[94] = new NISTMottRS(Element::Pu, 2);
      mScatter2[95] = new NISTMottRS(Element::Am, 2);
      mScatter2[96] = new NISTMottRS(Element::Cm, 2);

      mScatter3[1] = new NISTMottRS(Element::H, 3);
      mScatter3[2] = new NISTMottRS(Element::He, 3);
      mScatter3[3] = new NISTMottRS(Element::Li, 3);
      mScatter3[4] = new NISTMottRS(Element::Be, 3);
      mScatter3[5] = new NISTMottRS(Element::B, 3);
      mScatter3[6] = new NISTMottRS(Element::C, 3);
      mScatter3[7] = new NISTMottRS(Element::N, 3);
      mScatter3[8] = new NISTMottRS(Element::O, 3);
      mScatter3[9] = new NISTMottRS(Element::F, 3);
      mScatter3[10] = new NISTMottRS(Element::Ne, 3);
      mScatter3[11] = new NISTMottRS(Element::Na, 3);
      mScatter3[12] = new NISTMottRS(Element::Mg, 3);
      mScatter3[13] = new NISTMottRS(Element::Al, 3);
      mScatter3[14] = new NISTMottRS(Element::Si, 3);
      mScatter3[15] = new NISTMottRS(Element::P, 3);
      mScatter3[16] = new NISTMottRS(Element::S, 3);
      mScatter3[17] = new NISTMottRS(Element::Cl, 3);
      mScatter3[18] = new NISTMottRS(Element::Ar, 3);
      mScatter3[19] = new NISTMottRS(Element::K, 3);
      mScatter3[20] = new NISTMottRS(Element::Ca, 3);
      mScatter3[21] = new NISTMottRS(Element::Sc, 3);
      mScatter3[22] = new NISTMottRS(Element::Ti, 3);
      mScatter3[23] = new NISTMottRS(Element::V, 3);
      mScatter3[24] = new NISTMottRS(Element::Cr, 3);
      mScatter3[25] = new NISTMottRS(Element::Mn, 3);
      mScatter3[26] = new NISTMottRS(Element::Fe, 3);
      mScatter3[27] = new NISTMottRS(Element::Co, 3);
      mScatter3[28] = new NISTMottRS(Element::Ni, 3);
      mScatter3[29] = new NISTMottRS(Element::Cu, 3);
      mScatter3[30] = new NISTMottRS(Element::Zn, 3);
      mScatter3[31] = new NISTMottRS(Element::Ga, 3);
      mScatter3[32] = new NISTMottRS(Element::Ge, 3);
      mScatter3[33] = new NISTMottRS(Element::As, 3);
      mScatter3[34] = new NISTMottRS(Element::Se, 3);
      mScatter3[35] = new NISTMottRS(Element::Br, 3);
      mScatter3[36] = new NISTMottRS(Element::Kr, 3);
      mScatter3[37] = new NISTMottRS(Element::Rb, 3);
      mScatter3[38] = new NISTMottRS(Element::Sr, 3);
      mScatter3[39] = new NISTMottRS(Element::Y, 3);
      mScatter3[40] = new NISTMottRS(Element::Zr, 3);
      mScatter3[41] = new NISTMottRS(Element::Nb, 3);
      mScatter3[42] = new NISTMottRS(Element::Mo, 3);
      mScatter3[43] = new NISTMottRS(Element::Tc, 3);
      mScatter3[44] = new NISTMottRS(Element::Ru, 3);
      mScatter3[45] = new NISTMottRS(Element::Rh, 3);
      mScatter3[46] = new NISTMottRS(Element::Pd, 3);
      mScatter3[47] = new NISTMottRS(Element::Ag, 3);
      mScatter3[48] = new NISTMottRS(Element::Cd, 3);
      mScatter3[49] = new NISTMottRS(Element::In, 3);
      mScatter3[50] = new NISTMottRS(Element::Sn, 3);
      mScatter3[51] = new NISTMottRS(Element::Sb, 3);
      mScatter3[52] = new NISTMottRS(Element::Te, 3);
      mScatter3[53] = new NISTMottRS(Element::I, 3);
      mScatter3[54] = new NISTMottRS(Element::Xe, 3);
      mScatter3[55] = new NISTMottRS(Element::Cs, 3);
      mScatter3[56] = new NISTMottRS(Element::Ba, 3);
      mScatter3[57] = new NISTMottRS(Element::La, 3);
      mScatter3[58] = new NISTMottRS(Element::Ce, 3);
      mScatter3[59] = new NISTMottRS(Element::Pr, 3);
      mScatter3[60] = new NISTMottRS(Element::Nd, 3);
      mScatter3[61] = new NISTMottRS(Element::Pm, 3);
      mScatter3[62] = new NISTMottRS(Element::Sm, 3);
      mScatter3[63] = new NISTMottRS(Element::Eu, 3);
      mScatter3[64] = new NISTMottRS(Element::Gd, 3);
      mScatter3[65] = new NISTMottRS(Element::Tb, 3);
      mScatter3[66] = new NISTMottRS(Element::Dy, 3);
      mScatter3[67] = new NISTMottRS(Element::Ho, 3);
      mScatter3[68] = new NISTMottRS(Element::Er, 3);
      mScatter3[69] = new NISTMottRS(Element::Tm, 3);
      mScatter3[70] = new NISTMottRS(Element::Yb, 3);
      mScatter3[71] = new NISTMottRS(Element::Lu, 3);
      mScatter3[72] = new NISTMottRS(Element::Hf, 3);
      mScatter3[73] = new NISTMottRS(Element::Ta, 3);
      mScatter3[74] = new NISTMottRS(Element::W, 3);
      mScatter3[75] = new NISTMottRS(Element::Re, 3);
      mScatter3[76] = new NISTMottRS(Element::Os, 3);
      mScatter3[77] = new NISTMottRS(Element::Ir, 3);
      mScatter3[78] = new NISTMottRS(Element::Pt, 3);
      mScatter3[79] = new NISTMottRS(Element::Au, 3);
      mScatter3[80] = new NISTMottRS(Element::Hg, 3);
      mScatter3[81] = new NISTMottRS(Element::Tl, 3);
      mScatter3[82] = new NISTMottRS(Element::Pb, 3);
      mScatter3[83] = new NISTMottRS(Element::Bi, 3);
      mScatter3[84] = new NISTMottRS(Element::Po, 3);
      mScatter3[85] = new NISTMottRS(Element::At, 3);
      mScatter3[86] = new NISTMottRS(Element::Rn, 3);
      mScatter3[87] = new NISTMottRS(Element::Fr, 3);
      mScatter3[88] = new NISTMottRS(Element::Ra, 3);
      mScatter3[89] = new NISTMottRS(Element::Ac, 3);
      mScatter3[90] = new NISTMottRS(Element::Th, 3);
      mScatter3[91] = new NISTMottRS(Element::Pa, 3);
      mScatter3[92] = new NISTMottRS(Element::U, 3);
      mScatter3[93] = new NISTMottRS(Element::Np, 3);
      mScatter3[94] = new NISTMottRS(Element::Pu, 3);
      mScatter3[95] = new NISTMottRS(Element::Am, 3);
      mScatter3[96] = new NISTMottRS(Element::Cm, 3);
   }

   __global__ void initCuda()
   {
      d_mScatter1[1] = new NISTMottRS(*Element::dH, 1);
      d_mScatter1[2] = new NISTMottRS(*Element::dHe, 1);
      d_mScatter1[3] = new NISTMottRS(*Element::dLi, 1);
      d_mScatter1[4] = new NISTMottRS(*Element::dBe, 1);
      d_mScatter1[5] = new NISTMottRS(*Element::dB, 1);
      d_mScatter1[6] = new NISTMottRS(*Element::dC, 1);
      d_mScatter1[7] = new NISTMottRS(*Element::dN, 1);
      d_mScatter1[8] = new NISTMottRS(*Element::dO, 1);
      d_mScatter1[9] = new NISTMottRS(*Element::dF, 1);
      d_mScatter1[10] = new NISTMottRS(*Element::dNe, 1);
      d_mScatter1[11] = new NISTMottRS(*Element::dNa, 1);
      d_mScatter1[12] = new NISTMottRS(*Element::dMg, 1);
      d_mScatter1[13] = new NISTMottRS(*Element::dAl, 1);
      d_mScatter1[14] = new NISTMottRS(*Element::dSi, 1);
      d_mScatter1[15] = new NISTMottRS(*Element::dP, 1);
      d_mScatter1[16] = new NISTMottRS(*Element::dS, 1);
      d_mScatter1[17] = new NISTMottRS(*Element::dCl, 1);
      d_mScatter1[18] = new NISTMottRS(*Element::dAr, 1);
      d_mScatter1[19] = new NISTMottRS(*Element::dK, 1);
      d_mScatter1[20] = new NISTMottRS(*Element::dCa, 1);
      d_mScatter1[21] = new NISTMottRS(*Element::dSc, 1);
      d_mScatter1[22] = new NISTMottRS(*Element::dTi, 1);
      d_mScatter1[23] = new NISTMottRS(*Element::dV, 1);
      d_mScatter1[24] = new NISTMottRS(*Element::dCr, 1);
      d_mScatter1[25] = new NISTMottRS(*Element::dMn, 1);
      d_mScatter1[26] = new NISTMottRS(*Element::dFe, 1);
      d_mScatter1[27] = new NISTMottRS(*Element::dCo, 1);
      d_mScatter1[28] = new NISTMottRS(*Element::dNi, 1);
      d_mScatter1[29] = new NISTMottRS(*Element::dCu, 1);
      d_mScatter1[30] = new NISTMottRS(*Element::dZn, 1);
      d_mScatter1[31] = new NISTMottRS(*Element::dGa, 1);
      d_mScatter1[32] = new NISTMottRS(*Element::dGe, 1);
      d_mScatter1[33] = new NISTMottRS(*Element::dAs, 1);
      d_mScatter1[34] = new NISTMottRS(*Element::dSe, 1);
      d_mScatter1[35] = new NISTMottRS(*Element::dBr, 1);
      d_mScatter1[36] = new NISTMottRS(*Element::dKr, 1);
      d_mScatter1[37] = new NISTMottRS(*Element::dRb, 1);
      d_mScatter1[38] = new NISTMottRS(*Element::dSr, 1);
      d_mScatter1[39] = new NISTMottRS(*Element::dY, 1);
      d_mScatter1[40] = new NISTMottRS(*Element::dZr, 1);
      d_mScatter1[41] = new NISTMottRS(*Element::dNb, 1);
      d_mScatter1[42] = new NISTMottRS(*Element::dMo, 1);
      d_mScatter1[43] = new NISTMottRS(*Element::dTc, 1);
      d_mScatter1[44] = new NISTMottRS(*Element::dRu, 1);
      d_mScatter1[45] = new NISTMottRS(*Element::dRh, 1);
      d_mScatter1[46] = new NISTMottRS(*Element::dPd, 1);
      d_mScatter1[47] = new NISTMottRS(*Element::dAg, 1);
      d_mScatter1[48] = new NISTMottRS(*Element::dCd, 1);
      d_mScatter1[49] = new NISTMottRS(*Element::dIn, 1);
      d_mScatter1[50] = new NISTMottRS(*Element::dSn, 1);
      d_mScatter1[51] = new NISTMottRS(*Element::dSb, 1);
      d_mScatter1[52] = new NISTMottRS(*Element::dTe, 1);
      d_mScatter1[53] = new NISTMottRS(*Element::dI, 1);
      d_mScatter1[54] = new NISTMottRS(*Element::dXe, 1);
      d_mScatter1[55] = new NISTMottRS(*Element::dCs, 1);
      d_mScatter1[56] = new NISTMottRS(*Element::dBa, 1);
      d_mScatter1[57] = new NISTMottRS(*Element::dLa, 1);
      d_mScatter1[58] = new NISTMottRS(*Element::dCe, 1);
      d_mScatter1[59] = new NISTMottRS(*Element::dPr, 1);
      d_mScatter1[60] = new NISTMottRS(*Element::dNd, 1);
      d_mScatter1[61] = new NISTMottRS(*Element::dPm, 1);
      d_mScatter1[62] = new NISTMottRS(*Element::dSm, 1);
      d_mScatter1[63] = new NISTMottRS(*Element::dEu, 1);
      d_mScatter1[64] = new NISTMottRS(*Element::dGd, 1);
      d_mScatter1[65] = new NISTMottRS(*Element::dTb, 1);
      d_mScatter1[66] = new NISTMottRS(*Element::dDy, 1);
      d_mScatter1[67] = new NISTMottRS(*Element::dHo, 1);
      d_mScatter1[68] = new NISTMottRS(*Element::dEr, 1);
      d_mScatter1[69] = new NISTMottRS(*Element::dTm, 1);
      d_mScatter1[70] = new NISTMottRS(*Element::dYb, 1);
      d_mScatter1[71] = new NISTMottRS(*Element::dLu, 1);
      d_mScatter1[72] = new NISTMottRS(*Element::dHf, 1);
      d_mScatter1[73] = new NISTMottRS(*Element::dTa, 1);
      d_mScatter1[74] = new NISTMottRS(*Element::dW, 1);
      d_mScatter1[75] = new NISTMottRS(*Element::dRe, 1);
      d_mScatter1[76] = new NISTMottRS(*Element::dOs, 1);
      d_mScatter1[77] = new NISTMottRS(*Element::dIr, 1);
      d_mScatter1[78] = new NISTMottRS(*Element::dPt, 1);
      d_mScatter1[79] = new NISTMottRS(*Element::dAu, 1);
      d_mScatter1[80] = new NISTMottRS(*Element::dHg, 1);
      d_mScatter1[81] = new NISTMottRS(*Element::dTl, 1);
      d_mScatter1[82] = new NISTMottRS(*Element::dPb, 1);
      d_mScatter1[83] = new NISTMottRS(*Element::dBi, 1);
      d_mScatter1[84] = new NISTMottRS(*Element::dPo, 1);
      d_mScatter1[85] = new NISTMottRS(*Element::dAt, 1);
      d_mScatter1[86] = new NISTMottRS(*Element::dRn, 1);
      d_mScatter1[87] = new NISTMottRS(*Element::dFr, 1);
      d_mScatter1[88] = new NISTMottRS(*Element::dRa, 1);
      d_mScatter1[89] = new NISTMottRS(*Element::dAc, 1);
      d_mScatter1[90] = new NISTMottRS(*Element::dTh, 1);
      d_mScatter1[91] = new NISTMottRS(*Element::dPa, 1);
      d_mScatter1[92] = new NISTMottRS(*Element::dU, 1);
      d_mScatter1[93] = new NISTMottRS(*Element::dNp, 1);
      d_mScatter1[94] = new NISTMottRS(*Element::dPu, 1);
      d_mScatter1[95] = new NISTMottRS(*Element::dAm, 1);
      d_mScatter1[96] = new NISTMottRS(*Element::dCm, 1);

      d_mScatter2[1] = new NISTMottRS(*Element::dH, 2);
      d_mScatter2[2] = new NISTMottRS(*Element::dHe, 2);
      d_mScatter2[3] = new NISTMottRS(*Element::dLi, 2);
      d_mScatter2[4] = new NISTMottRS(*Element::dBe, 2);
      d_mScatter2[5] = new NISTMottRS(*Element::dB, 2);
      d_mScatter2[6] = new NISTMottRS(*Element::dC, 2);
      d_mScatter2[7] = new NISTMottRS(*Element::dN, 2);
      d_mScatter2[8] = new NISTMottRS(*Element::dO, 2);
      d_mScatter2[9] = new NISTMottRS(*Element::dF, 2);
      d_mScatter2[10] = new NISTMottRS(*Element::dNe, 2);
      d_mScatter2[11] = new NISTMottRS(*Element::dNa, 2);
      d_mScatter2[12] = new NISTMottRS(*Element::dMg, 2);
      d_mScatter2[13] = new NISTMottRS(*Element::dAl, 2);
      d_mScatter2[14] = new NISTMottRS(*Element::dSi, 2);
      d_mScatter2[15] = new NISTMottRS(*Element::dP, 2);
      d_mScatter2[16] = new NISTMottRS(*Element::dS, 2);
      d_mScatter2[17] = new NISTMottRS(*Element::dCl, 2);
      d_mScatter2[18] = new NISTMottRS(*Element::dAr, 2);
      d_mScatter2[19] = new NISTMottRS(*Element::dK, 2);
      d_mScatter2[20] = new NISTMottRS(*Element::dCa, 2);
      d_mScatter2[21] = new NISTMottRS(*Element::dSc, 2);
      d_mScatter2[22] = new NISTMottRS(*Element::dTi, 2);
      d_mScatter2[23] = new NISTMottRS(*Element::dV, 2);
      d_mScatter2[24] = new NISTMottRS(*Element::dCr, 2);
      d_mScatter2[25] = new NISTMottRS(*Element::dMn, 2);
      d_mScatter2[26] = new NISTMottRS(*Element::dFe, 2);
      d_mScatter2[27] = new NISTMottRS(*Element::dCo, 2);
      d_mScatter2[28] = new NISTMottRS(*Element::dNi, 2);
      d_mScatter2[29] = new NISTMottRS(*Element::dCu, 2);
      d_mScatter2[30] = new NISTMottRS(*Element::dZn, 2);
      d_mScatter2[31] = new NISTMottRS(*Element::dGa, 2);
      d_mScatter2[32] = new NISTMottRS(*Element::dGe, 2);
      d_mScatter2[33] = new NISTMottRS(*Element::dAs, 2);
      d_mScatter2[34] = new NISTMottRS(*Element::dSe, 2);
      d_mScatter2[35] = new NISTMottRS(*Element::dBr, 2);
      d_mScatter2[36] = new NISTMottRS(*Element::dKr, 2);
      d_mScatter2[37] = new NISTMottRS(*Element::dRb, 2);
      d_mScatter2[38] = new NISTMottRS(*Element::dSr, 2);
      d_mScatter2[39] = new NISTMottRS(*Element::dY, 2);
      d_mScatter2[40] = new NISTMottRS(*Element::dZr, 2);
      d_mScatter2[41] = new NISTMottRS(*Element::dNb, 2);
      d_mScatter2[42] = new NISTMottRS(*Element::dMo, 2);
      d_mScatter2[43] = new NISTMottRS(*Element::dTc, 2);
      d_mScatter2[44] = new NISTMottRS(*Element::dRu, 2);
      d_mScatter2[45] = new NISTMottRS(*Element::dRh, 2);
      d_mScatter2[46] = new NISTMottRS(*Element::dPd, 2);
      d_mScatter2[47] = new NISTMottRS(*Element::dAg, 2);
      d_mScatter2[48] = new NISTMottRS(*Element::dCd, 2);
      d_mScatter2[49] = new NISTMottRS(*Element::dIn, 2);
      d_mScatter2[50] = new NISTMottRS(*Element::dSn, 2);
      d_mScatter2[51] = new NISTMottRS(*Element::dSb, 2);
      d_mScatter2[52] = new NISTMottRS(*Element::dTe, 2);
      d_mScatter2[53] = new NISTMottRS(*Element::dI, 2);
      d_mScatter2[54] = new NISTMottRS(*Element::dXe, 2);
      d_mScatter2[55] = new NISTMottRS(*Element::dCs, 2);
      d_mScatter2[56] = new NISTMottRS(*Element::dBa, 2);
      d_mScatter2[57] = new NISTMottRS(*Element::dLa, 2);
      d_mScatter2[58] = new NISTMottRS(*Element::dCe, 2);
      d_mScatter2[59] = new NISTMottRS(*Element::dPr, 2);
      d_mScatter2[60] = new NISTMottRS(*Element::dNd, 2);
      d_mScatter2[61] = new NISTMottRS(*Element::dPm, 2);
      d_mScatter2[62] = new NISTMottRS(*Element::dSm, 2);
      d_mScatter2[63] = new NISTMottRS(*Element::dEu, 2);
      d_mScatter2[64] = new NISTMottRS(*Element::dGd, 2);
      d_mScatter2[65] = new NISTMottRS(*Element::dTb, 2);
      d_mScatter2[66] = new NISTMottRS(*Element::dDy, 2);
      d_mScatter2[67] = new NISTMottRS(*Element::dHo, 2);
      d_mScatter2[68] = new NISTMottRS(*Element::dEr, 2);
      d_mScatter2[69] = new NISTMottRS(*Element::dTm, 2);
      d_mScatter2[70] = new NISTMottRS(*Element::dYb, 2);
      d_mScatter2[71] = new NISTMottRS(*Element::dLu, 2);
      d_mScatter2[72] = new NISTMottRS(*Element::dHf, 2);
      d_mScatter2[73] = new NISTMottRS(*Element::dTa, 2);
      d_mScatter2[74] = new NISTMottRS(*Element::dW, 2);
      d_mScatter2[75] = new NISTMottRS(*Element::dRe, 2);
      d_mScatter2[76] = new NISTMottRS(*Element::dOs, 2);
      d_mScatter2[77] = new NISTMottRS(*Element::dIr, 2);
      d_mScatter2[78] = new NISTMottRS(*Element::dPt, 2);
      d_mScatter2[79] = new NISTMottRS(*Element::dAu, 2);
      d_mScatter2[80] = new NISTMottRS(*Element::dHg, 2);
      d_mScatter2[81] = new NISTMottRS(*Element::dTl, 2);
      d_mScatter2[82] = new NISTMottRS(*Element::dPb, 2);
      d_mScatter2[83] = new NISTMottRS(*Element::dBi, 2);
      d_mScatter2[84] = new NISTMottRS(*Element::dPo, 2);
      d_mScatter2[85] = new NISTMottRS(*Element::dAt, 2);
      d_mScatter2[86] = new NISTMottRS(*Element::dRn, 2);
      d_mScatter2[87] = new NISTMottRS(*Element::dFr, 2);
      d_mScatter2[88] = new NISTMottRS(*Element::dRa, 2);
      d_mScatter2[89] = new NISTMottRS(*Element::dAc, 2);
      d_mScatter2[90] = new NISTMottRS(*Element::dTh, 2);
      d_mScatter2[91] = new NISTMottRS(*Element::dPa, 2);
      d_mScatter2[92] = new NISTMottRS(*Element::dU, 2);
      d_mScatter2[93] = new NISTMottRS(*Element::dNp, 2);
      d_mScatter2[94] = new NISTMottRS(*Element::dPu, 2);
      d_mScatter2[95] = new NISTMottRS(*Element::dAm, 2);
      d_mScatter2[96] = new NISTMottRS(*Element::dCm, 2);

      d_mScatter3[1] = new NISTMottRS(*Element::dH, 3);
      d_mScatter3[2] = new NISTMottRS(*Element::dHe, 3);
      d_mScatter3[3] = new NISTMottRS(*Element::dLi, 3);
      d_mScatter3[4] = new NISTMottRS(*Element::dBe, 3);
      d_mScatter3[5] = new NISTMottRS(*Element::dB, 3);
      d_mScatter3[6] = new NISTMottRS(*Element::dC, 3);
      d_mScatter3[7] = new NISTMottRS(*Element::dN, 3);
      d_mScatter3[8] = new NISTMottRS(*Element::dO, 3);
      d_mScatter3[9] = new NISTMottRS(*Element::dF, 3);
      d_mScatter3[10] = new NISTMottRS(*Element::dNe, 3);
      d_mScatter3[11] = new NISTMottRS(*Element::dNa, 3);
      d_mScatter3[12] = new NISTMottRS(*Element::dMg, 3);
      d_mScatter3[13] = new NISTMottRS(*Element::dAl, 3);
      d_mScatter3[14] = new NISTMottRS(*Element::dSi, 3);
      d_mScatter3[15] = new NISTMottRS(*Element::dP, 3);
      d_mScatter3[16] = new NISTMottRS(*Element::dS, 3);
      d_mScatter3[17] = new NISTMottRS(*Element::dCl, 3);
      d_mScatter3[18] = new NISTMottRS(*Element::dAr, 3);
      d_mScatter3[19] = new NISTMottRS(*Element::dK, 3);
      d_mScatter3[20] = new NISTMottRS(*Element::dCa, 3);
      d_mScatter3[21] = new NISTMottRS(*Element::dSc, 3);
      d_mScatter3[22] = new NISTMottRS(*Element::dTi, 3);
      d_mScatter3[23] = new NISTMottRS(*Element::dV, 3);
      d_mScatter3[24] = new NISTMottRS(*Element::dCr, 3);
      d_mScatter3[25] = new NISTMottRS(*Element::dMn, 3);
      d_mScatter3[26] = new NISTMottRS(*Element::dFe, 3);
      d_mScatter3[27] = new NISTMottRS(*Element::dCo, 3);
      d_mScatter3[28] = new NISTMottRS(*Element::dNi, 3);
      d_mScatter3[29] = new NISTMottRS(*Element::dCu, 3);
      d_mScatter3[30] = new NISTMottRS(*Element::dZn, 3);
      d_mScatter3[31] = new NISTMottRS(*Element::dGa, 3);
      d_mScatter3[32] = new NISTMottRS(*Element::dGe, 3);
      d_mScatter3[33] = new NISTMottRS(*Element::dAs, 3);
      d_mScatter3[34] = new NISTMottRS(*Element::dSe, 3);
      d_mScatter3[35] = new NISTMottRS(*Element::dBr, 3);
      d_mScatter3[36] = new NISTMottRS(*Element::dKr, 3);
      d_mScatter3[37] = new NISTMottRS(*Element::dRb, 3);
      d_mScatter3[38] = new NISTMottRS(*Element::dSr, 3);
      d_mScatter3[39] = new NISTMottRS(*Element::dY, 3);
      d_mScatter3[40] = new NISTMottRS(*Element::dZr, 3);
      d_mScatter3[41] = new NISTMottRS(*Element::dNb, 3);
      d_mScatter3[42] = new NISTMottRS(*Element::dMo, 3);
      d_mScatter3[43] = new NISTMottRS(*Element::dTc, 3);
      d_mScatter3[44] = new NISTMottRS(*Element::dRu, 3);
      d_mScatter3[45] = new NISTMottRS(*Element::dRh, 3);
      d_mScatter3[46] = new NISTMottRS(*Element::dPd, 3);
      d_mScatter3[47] = new NISTMottRS(*Element::dAg, 3);
      d_mScatter3[48] = new NISTMottRS(*Element::dCd, 3);
      d_mScatter3[49] = new NISTMottRS(*Element::dIn, 3);
      d_mScatter3[50] = new NISTMottRS(*Element::dSn, 3);
      d_mScatter3[51] = new NISTMottRS(*Element::dSb, 3);
      d_mScatter3[52] = new NISTMottRS(*Element::dTe, 3);
      d_mScatter3[53] = new NISTMottRS(*Element::dI, 3);
      d_mScatter3[54] = new NISTMottRS(*Element::dXe, 3);
      d_mScatter3[55] = new NISTMottRS(*Element::dCs, 3);
      d_mScatter3[56] = new NISTMottRS(*Element::dBa, 3);
      d_mScatter3[57] = new NISTMottRS(*Element::dLa, 3);
      d_mScatter3[58] = new NISTMottRS(*Element::dCe, 3);
      d_mScatter3[59] = new NISTMottRS(*Element::dPr, 3);
      d_mScatter3[60] = new NISTMottRS(*Element::dNd, 3);
      d_mScatter3[61] = new NISTMottRS(*Element::dPm, 3);
      d_mScatter3[62] = new NISTMottRS(*Element::dSm, 3);
      d_mScatter3[63] = new NISTMottRS(*Element::dEu, 3);
      d_mScatter3[64] = new NISTMottRS(*Element::dGd, 3);
      d_mScatter3[65] = new NISTMottRS(*Element::dTb, 3);
      d_mScatter3[66] = new NISTMottRS(*Element::dDy, 3);
      d_mScatter3[67] = new NISTMottRS(*Element::dHo, 3);
      d_mScatter3[68] = new NISTMottRS(*Element::dEr, 3);
      d_mScatter3[69] = new NISTMottRS(*Element::dTm, 3);
      d_mScatter3[70] = new NISTMottRS(*Element::dYb, 3);
      d_mScatter3[71] = new NISTMottRS(*Element::dLu, 3);
      d_mScatter3[72] = new NISTMottRS(*Element::dHf, 3);
      d_mScatter3[73] = new NISTMottRS(*Element::dTa, 3);
      d_mScatter3[74] = new NISTMottRS(*Element::dW, 3);
      d_mScatter3[75] = new NISTMottRS(*Element::dRe, 3);
      d_mScatter3[76] = new NISTMottRS(*Element::dOs, 3);
      d_mScatter3[77] = new NISTMottRS(*Element::dIr, 3);
      d_mScatter3[78] = new NISTMottRS(*Element::dPt, 3);
      d_mScatter3[79] = new NISTMottRS(*Element::dAu, 3);
      d_mScatter3[80] = new NISTMottRS(*Element::dHg, 3);
      d_mScatter3[81] = new NISTMottRS(*Element::dTl, 3);
      d_mScatter3[82] = new NISTMottRS(*Element::dPb, 3);
      d_mScatter3[83] = new NISTMottRS(*Element::dBi, 3);
      d_mScatter3[84] = new NISTMottRS(*Element::dPo, 3);
      d_mScatter3[85] = new NISTMottRS(*Element::dAt, 3);
      d_mScatter3[86] = new NISTMottRS(*Element::dRn, 3);
      d_mScatter3[87] = new NISTMottRS(*Element::dFr, 3);
      d_mScatter3[88] = new NISTMottRS(*Element::dRa, 3);
      d_mScatter3[89] = new NISTMottRS(*Element::dAc, 3);
      d_mScatter3[90] = new NISTMottRS(*Element::dTh, 3);
      d_mScatter3[91] = new NISTMottRS(*Element::dPa, 3);
      d_mScatter3[92] = new NISTMottRS(*Element::dU, 3);
      d_mScatter3[93] = new NISTMottRS(*Element::dNp, 3);
      d_mScatter3[94] = new NISTMottRS(*Element::dPu, 3);
      d_mScatter3[95] = new NISTMottRS(*Element::dAm, 3);
      d_mScatter3[96] = new NISTMottRS(*Element::dCm, 3);
   }
}