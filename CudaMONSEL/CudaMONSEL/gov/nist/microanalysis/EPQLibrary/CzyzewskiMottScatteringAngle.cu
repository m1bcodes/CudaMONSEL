#include "gov\nist\microanalysis\EPQLibrary\CzyzewskiMottScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\CzyzewskiMottCrossSection.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Reference.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ScreenedRutherfordScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

#include "Amphibian\Algorithm.cuh"
#include "Amphibian\random.cuh"

#include "CudaUtil.h"

namespace CzyzewskiMottScatteringAngle
{
   static const Reference::Author* alRef[] = { &Reference::Czyzewski, &Reference::MacCallum, &Reference::DJoy };
   static const Reference::JournalArticle REFERENCE(Reference::JApplPhys, "68, No. 7", "", 1990, alRef, 3);

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const double MAX_CZYZEWSKI = 4.8065296e-15;
#else
   const double MAX_CZYZEWSKI = ToSI::keV(30.0);
#endif

   void CzyzewskiMottScatteringAngle::init(int an)
   {
      // Use the MottCrossSection then discard it
      CzyzewskiMottCrossSectionT mcs(an);
      mMeanFreePath.resize(CzyzewskiMottCrossSection::SpecialEnergyCount);
      mTotalCrossSection.resize(CzyzewskiMottCrossSection::SpecialEnergyCount);
      // Extract some useful stuff...
      for (int i = 0; i < CzyzewskiMottCrossSection::SpecialEnergyCount; ++i) {
         double e = CzyzewskiMottCrossSection::getSpecialEnergy(i);
         mMeanFreePath[i] = mcs.meanFreePath(e);
         mTotalCrossSection[i] = mcs.totalCrossSection(e);
      }
      // Calculate a normalized running sum that will be used to map a random
      // number between
      // [0,1.00) onto a scattering angle.
      mCummulativeDF.resize(CzyzewskiMottCrossSection::SpecialEnergyCount, VectorXd(CzyzewskiMottCrossSection::SpecialAngleCount, 0));
      for (int r = 0; r < CzyzewskiMottCrossSection::SpecialEnergyCount; ++r) {
         const double energy = CzyzewskiMottCrossSection::getSpecialEnergy(r);
         mCummulativeDF[r][0] = 0.0;
         for (int c = 1; c < CzyzewskiMottCrossSection::SpecialAngleCount; ++c) {
            const double cm = mcs.partialCrossSection(CzyzewskiMottCrossSection::getSpecialAngle(c - 1), energy);
            const double cp = mcs.partialCrossSection(CzyzewskiMottCrossSection::getSpecialAngle(c), energy);
            const double am = CzyzewskiMottCrossSection::getSpecialAngle(c - 1);
            const double ap = CzyzewskiMottCrossSection::getSpecialAngle(c);
            mCummulativeDF[r][c] = mCummulativeDF[r][c - 1] + ((cp + cm) * (ap - am)) / 2.0;
         }
         // Normalize to 1.0
         for (int c = 0; c < CzyzewskiMottCrossSection::SpecialAngleCount; ++c)
            mCummulativeDF[r][c] = mCummulativeDF[r][c] / mCummulativeDF[r][CzyzewskiMottCrossSection::SpecialAngleCount - 1];
      }
   }

   __host__ __device__ CzyzewskiMottScatteringAngle::CzyzewskiMottScatteringAngle(const ElementT& el) :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      RandomizedScatterT("Cyzewski", *Reference::dNullReference), mElement(el), mRutherford(ScreenedRutherfordScatteringAngle::getSRSA(el.getAtomicNumber()))
#else
      RandomizedScatterT("Cyzewski", REFERENCE), mElement(el), mRutherford(ScreenedRutherfordScatteringAngle::getSRSA(el.getAtomicNumber()))
#endif
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
#else
      init(el.getAtomicNumber());
#endif
   }

   StringT CzyzewskiMottScatteringAngle::toString() const
   {
      return "CrossSection[Czyzewski-Mott," + StringT(mElement.toAbbrev()) + "]";
   }


   const ElementT& CzyzewskiMottScatteringAngle::getElement() const
   {
      return mElement;
   }

   __host__ __device__ double CzyzewskiMottScatteringAngle::scatteringAngleForSpecialEnergy(int ei, double rand) const
   {
      VectorXd r(mCummulativeDF[ei]);
      int ai = Algorithm::binarySearch(r.data(), 0, r.size() - 1, rand);
      if (ai >= 0)
         return CzyzewskiMottCrossSection::getSpecialAngle(ai);
      else { // Interpolate between angles
         ai = -(ai + 1);
         if (!(ai >= 1)) printf("CzyzewskiMottScatteringAngle::scatteringAngleForSpecialEnergy 1: %d\n", ai);
         if (!(rand <= r[ai])) printf("CzyzewskiMottScatteringAngle::scatteringAngleForSpecialEnergy 2: %lf, %lf\n", rand, r[ai]);
         if (!(rand > r[ai - 1]))printf("CzyzewskiMottScatteringAngle::scatteringAngleForSpecialEnergy 3: %lf, %lf\n", rand, r[ai - 1]);;
         double am = CzyzewskiMottCrossSection::getSpecialAngle(ai - 1);
         return am + (((rand - r[ai - 1]) / (r[ai] - r[ai - 1])) * (CzyzewskiMottCrossSection::getSpecialAngle(ai) - am));
      }
   }

   double CzyzewskiMottScatteringAngle::randomScatteringAngle(double energy, double rand) const
   {
      if (!(rand >= 0)) printf("CzyzewskiMottScatteringAngle::randomScatteringAngle 1: %lf\n", rand);
      if (!(rand < 1.0)) printf("CzyzewskiMottScatteringAngle::randomScatteringAngle 2: %lf\n", rand);
      if (!(energy <= ToSI::keV(30.0))) printf("CzyzewskiMottScatteringAngle::randomScatteringAngle 3: %lf\n", energy);
      if (!(energy > 0.0)) printf("CzyzewskiMottScatteringAngle::randomScatteringAngle 4: %lf\n", energy);
      const int e = CzyzewskiMottCrossSection::getEnergyIndex(energy);
      const double e0 = CzyzewskiMottCrossSection::getSpecialEnergy(e);
      if (!(energy <= e0)) printf("CzyzewskiMottScatteringAngle::randomScatteringAngle 5: %lf\n", energy);
      const double a0 = scatteringAngleForSpecialEnergy(e, rand);
      if (energy == e0)
         return a0;
      else {
         const double a1 = scatteringAngleForSpecialEnergy(e - 1, rand);
         const double e1 = CzyzewskiMottCrossSection::getSpecialEnergy(e - 1);
         if (!(energy > e1)) printf("CzyzewskiMottScatteringAngle::randomScatteringAngle 6: %lf\n", energy);
         if (!(a1 >= 0.0)) printf("CzyzewskiMottScatteringAngle::randomScatteringAngle 7: %lf\n", energy);
         if (!(a0 >= 0.0)) printf("CzyzewskiMottScatteringAngle::randomScatteringAngle 8: %lf\n", energy);
         return a0 + ((energy - e0) / (e1 - e0)) * (a1 - a0);
      }
   }

   double CzyzewskiMottScatteringAngle::randomScatteringAngle(double energy) const
   {
      if (energy < MAX_CZYZEWSKI) {
         return randomScatteringAngle(energy, Random::random());
      }
      else {
         return mRutherford.randomScatteringAngle(energy);
      }
   }

   double CzyzewskiMottScatteringAngle::meanFreePath(double energy) const
   {
      int e = CzyzewskiMottCrossSection::getEnergyIndex(energy);
      if (e == 0)
         e = 1;
      const double e0 = CzyzewskiMottCrossSection::getSpecialEnergy(e - 1);
      const double e1 = CzyzewskiMottCrossSection::getSpecialEnergy(e);
      if (!(energy >= e0)) printf("CzyzewskiMottScatteringAngle::meanFreePath 1: %lf\n", energy);
      if (!(energy <= e1)) printf("CzyzewskiMottScatteringAngle::meanFreePath 2: %lf\n", energy);
      return mMeanFreePath[e - 1] + ((mMeanFreePath[e] - mMeanFreePath[e - 1]) * ((energy - e0) / (e1 - e0)));
   }

   __host__ __device__ double CzyzewskiMottScatteringAngle::totalCrossSection(double energy) const
   {
      if (energy < MAX_CZYZEWSKI) {
         int e = CzyzewskiMottCrossSection::getEnergyIndex(energy);
         if (e == 0)
            e = 1;
         const double e0 = CzyzewskiMottCrossSection::getSpecialEnergy(e - 1);
         const double e1 = CzyzewskiMottCrossSection::getSpecialEnergy(e);
         if (!(energy >= e0)) printf("CzyzewskiMottScatteringAngle::totalCrossSection 1: %lf\n", energy);
         if (!(energy <= e1)) printf("CzyzewskiMottScatteringAngle::totalCrossSection 2: %lf\n", energy);
         return mTotalCrossSection[e - 1] + ((mTotalCrossSection[e] - mTotalCrossSection[e - 1]) * ((energy - e0) / (e1 - e0)));
      }
      else {
         return mRutherford.totalCrossSection(energy);
      }
   }

   //const CzyzewskiMottScatteringAngle CMSA1(Element::H);
   //const CzyzewskiMottScatteringAngle CMSA2(Element::He);
   //const CzyzewskiMottScatteringAngle CMSA3(Element::Li);
   //const CzyzewskiMottScatteringAngle CMSA4(Element::Be);
   //const CzyzewskiMottScatteringAngle CMSA5(Element::B);
   //const CzyzewskiMottScatteringAngle CMSA6(Element::C);
   //const CzyzewskiMottScatteringAngle CMSA7(Element::N);
   //const CzyzewskiMottScatteringAngle CMSA8(Element::O);
   //const CzyzewskiMottScatteringAngle CMSA9(Element::F);
   //const CzyzewskiMottScatteringAngle CMSA10(Element::Ne);
   //const CzyzewskiMottScatteringAngle CMSA11(Element::Na);
   //const CzyzewskiMottScatteringAngle CMSA12(Element::Mg);
   //const CzyzewskiMottScatteringAngle CMSA13(Element::Al);
   //const CzyzewskiMottScatteringAngle CMSA14(Element::Si);
   //const CzyzewskiMottScatteringAngle CMSA15(Element::P);
   //const CzyzewskiMottScatteringAngle CMSA16(Element::S);
   //const CzyzewskiMottScatteringAngle CMSA17(Element::Cl);
   //const CzyzewskiMottScatteringAngle CMSA18(Element::Ar);
   //const CzyzewskiMottScatteringAngle CMSA19(Element::K);
   //const CzyzewskiMottScatteringAngle CMSA20(Element::Ca);
   //const CzyzewskiMottScatteringAngle CMSA21(Element::Sc);
   //const CzyzewskiMottScatteringAngle CMSA22(Element::Ti);
   //const CzyzewskiMottScatteringAngle CMSA23(Element::V);
   //const CzyzewskiMottScatteringAngle CMSA24(Element::Cr);
   //const CzyzewskiMottScatteringAngle CMSA25(Element::Mn);
   //const CzyzewskiMottScatteringAngle CMSA26(Element::Fe);
   //const CzyzewskiMottScatteringAngle CMSA27(Element::Co);
   //const CzyzewskiMottScatteringAngle CMSA28(Element::Ni);
   //const CzyzewskiMottScatteringAngle CMSA29(Element::Cu);
   //const CzyzewskiMottScatteringAngle CMSA30(Element::Zn);
   //const CzyzewskiMottScatteringAngle CMSA31(Element::Ga);
   //const CzyzewskiMottScatteringAngle CMSA32(Element::Ge);
   //const CzyzewskiMottScatteringAngle CMSA33(Element::As);
   //const CzyzewskiMottScatteringAngle CMSA34(Element::Se);
   //const CzyzewskiMottScatteringAngle CMSA35(Element::Br);
   //const CzyzewskiMottScatteringAngle CMSA36(Element::Kr);
   //const CzyzewskiMottScatteringAngle CMSA37(Element::Rb);
   //const CzyzewskiMottScatteringAngle CMSA38(Element::Sr);
   //const CzyzewskiMottScatteringAngle CMSA39(Element::Y);
   //const CzyzewskiMottScatteringAngle CMSA40(Element::Zr);
   //const CzyzewskiMottScatteringAngle CMSA41(Element::Nb);
   //const CzyzewskiMottScatteringAngle CMSA42(Element::Mo);
   //const CzyzewskiMottScatteringAngle CMSA43(Element::Tc);
   //const CzyzewskiMottScatteringAngle CMSA44(Element::Ru);
   //const CzyzewskiMottScatteringAngle CMSA45(Element::Rh);
   //const CzyzewskiMottScatteringAngle CMSA46(Element::Pd);
   //const CzyzewskiMottScatteringAngle CMSA47(Element::Ag);
   //const CzyzewskiMottScatteringAngle CMSA48(Element::Cd);
   //const CzyzewskiMottScatteringAngle CMSA49(Element::In);
   //const CzyzewskiMottScatteringAngle CMSA50(Element::Sn);
   //const CzyzewskiMottScatteringAngle CMSA51(Element::Sb);
   //const CzyzewskiMottScatteringAngle CMSA52(Element::Te);
   //const CzyzewskiMottScatteringAngle CMSA53(Element::I);
   //const CzyzewskiMottScatteringAngle CMSA54(Element::Xe);
   //const CzyzewskiMottScatteringAngle CMSA55(Element::Cs);
   //const CzyzewskiMottScatteringAngle CMSA56(Element::Ba);
   //const CzyzewskiMottScatteringAngle CMSA57(Element::La);
   //const CzyzewskiMottScatteringAngle CMSA58(Element::Ce);
   //const CzyzewskiMottScatteringAngle CMSA59(Element::Pr);
   //const CzyzewskiMottScatteringAngle CMSA60(Element::Nd);
   //const CzyzewskiMottScatteringAngle CMSA61(Element::Pm);
   //const CzyzewskiMottScatteringAngle CMSA62(Element::Sm);
   //const CzyzewskiMottScatteringAngle CMSA63(Element::Eu);
   //const CzyzewskiMottScatteringAngle CMSA64(Element::Gd);
   //const CzyzewskiMottScatteringAngle CMSA65(Element::Tb);
   //const CzyzewskiMottScatteringAngle CMSA66(Element::Dy);
   //const CzyzewskiMottScatteringAngle CMSA67(Element::Ho);
   //const CzyzewskiMottScatteringAngle CMSA68(Element::Er);
   //const CzyzewskiMottScatteringAngle CMSA69(Element::Tm);
   //const CzyzewskiMottScatteringAngle CMSA70(Element::Yb);
   //const CzyzewskiMottScatteringAngle CMSA71(Element::Lu);
   //const CzyzewskiMottScatteringAngle CMSA72(Element::Hf);
   //const CzyzewskiMottScatteringAngle CMSA73(Element::Ta);
   //const CzyzewskiMottScatteringAngle CMSA74(Element::W);
   //const CzyzewskiMottScatteringAngle CMSA75(Element::Re);
   //const CzyzewskiMottScatteringAngle CMSA76(Element::Os);
   //const CzyzewskiMottScatteringAngle CMSA77(Element::Ir);
   //const CzyzewskiMottScatteringAngle CMSA78(Element::Pt);
   //const CzyzewskiMottScatteringAngle CMSA79(Element::Au);
   //const CzyzewskiMottScatteringAngle CMSA80(Element::Hg);
   //const CzyzewskiMottScatteringAngle CMSA81(Element::Tl);
   //const CzyzewskiMottScatteringAngle CMSA82(Element::Pb);
   //const CzyzewskiMottScatteringAngle CMSA83(Element::Bi);
   //const CzyzewskiMottScatteringAngle CMSA84(Element::Po);
   //const CzyzewskiMottScatteringAngle CMSA85(Element::At);
   //const CzyzewskiMottScatteringAngle CMSA86(Element::Rn);
   //const CzyzewskiMottScatteringAngle CMSA87(Element::Fr);
   //const CzyzewskiMottScatteringAngle CMSA88(Element::Ra);
   //const CzyzewskiMottScatteringAngle CMSA89(Element::Ac);
   //const CzyzewskiMottScatteringAngle CMSA90(Element::Th);
   //const CzyzewskiMottScatteringAngle CMSA91(Element::Pa);
   //const CzyzewskiMottScatteringAngle CMSA92(Element::U);
   //const CzyzewskiMottScatteringAngle CMSA93(Element::Np);
   //const CzyzewskiMottScatteringAngle CMSA94(Element::Pu);

   CzyzewskiMottScatteringAngle const * mScatter[113];

   void init()
   {
      mScatter[1] = new CzyzewskiMottScatteringAngle(Element::H);
      mScatter[2] = new CzyzewskiMottScatteringAngle(Element::He);
      mScatter[3] = new CzyzewskiMottScatteringAngle(Element::Li);
      mScatter[4] = new CzyzewskiMottScatteringAngle(Element::Be);
      mScatter[5] = new CzyzewskiMottScatteringAngle(Element::B);
      mScatter[6] = new CzyzewskiMottScatteringAngle(Element::C);
      mScatter[7] = new CzyzewskiMottScatteringAngle(Element::N);
      mScatter[8] = new CzyzewskiMottScatteringAngle(Element::O);
      mScatter[9] = new CzyzewskiMottScatteringAngle(Element::F);
      mScatter[10] = new CzyzewskiMottScatteringAngle(Element::Ne);
      mScatter[11] = new CzyzewskiMottScatteringAngle(Element::Na);
      mScatter[12] = new CzyzewskiMottScatteringAngle(Element::Mg);
      mScatter[13] = new CzyzewskiMottScatteringAngle(Element::Al);
      mScatter[14] = new CzyzewskiMottScatteringAngle(Element::Si);
      mScatter[15] = new CzyzewskiMottScatteringAngle(Element::P);
      mScatter[16] = new CzyzewskiMottScatteringAngle(Element::S);
      mScatter[17] = new CzyzewskiMottScatteringAngle(Element::Cl);
      mScatter[18] = new CzyzewskiMottScatteringAngle(Element::Ar);
      mScatter[19] = new CzyzewskiMottScatteringAngle(Element::K);
      mScatter[20] = new CzyzewskiMottScatteringAngle(Element::Ca);
      mScatter[21] = new CzyzewskiMottScatteringAngle(Element::Sc);
      mScatter[22] = new CzyzewskiMottScatteringAngle(Element::Ti);
      mScatter[23] = new CzyzewskiMottScatteringAngle(Element::V);
      mScatter[24] = new CzyzewskiMottScatteringAngle(Element::Cr);
      mScatter[25] = new CzyzewskiMottScatteringAngle(Element::Mn);
      mScatter[26] = new CzyzewskiMottScatteringAngle(Element::Fe);
      mScatter[27] = new CzyzewskiMottScatteringAngle(Element::Co);
      mScatter[28] = new CzyzewskiMottScatteringAngle(Element::Ni);
      mScatter[29] = new CzyzewskiMottScatteringAngle(Element::Cu);
      mScatter[30] = new CzyzewskiMottScatteringAngle(Element::Zn);
      mScatter[31] = new CzyzewskiMottScatteringAngle(Element::Ga);
      mScatter[32] = new CzyzewskiMottScatteringAngle(Element::Ge);
      mScatter[33] = new CzyzewskiMottScatteringAngle(Element::As);
      mScatter[34] = new CzyzewskiMottScatteringAngle(Element::Se);
      mScatter[35] = new CzyzewskiMottScatteringAngle(Element::Br);
      mScatter[36] = new CzyzewskiMottScatteringAngle(Element::Kr);
      mScatter[37] = new CzyzewskiMottScatteringAngle(Element::Rb);
      mScatter[38] = new CzyzewskiMottScatteringAngle(Element::Sr);
      mScatter[39] = new CzyzewskiMottScatteringAngle(Element::Y);
      mScatter[40] = new CzyzewskiMottScatteringAngle(Element::Zr);
      mScatter[41] = new CzyzewskiMottScatteringAngle(Element::Nb);
      mScatter[42] = new CzyzewskiMottScatteringAngle(Element::Mo);
      mScatter[43] = new CzyzewskiMottScatteringAngle(Element::Tc);
      mScatter[44] = new CzyzewskiMottScatteringAngle(Element::Ru);
      mScatter[45] = new CzyzewskiMottScatteringAngle(Element::Rh);
      mScatter[46] = new CzyzewskiMottScatteringAngle(Element::Pd);
      mScatter[47] = new CzyzewskiMottScatteringAngle(Element::Ag);
      mScatter[48] = new CzyzewskiMottScatteringAngle(Element::Cd);
      mScatter[49] = new CzyzewskiMottScatteringAngle(Element::In);
      mScatter[50] = new CzyzewskiMottScatteringAngle(Element::Sn);
      mScatter[51] = new CzyzewskiMottScatteringAngle(Element::Sb);
      mScatter[52] = new CzyzewskiMottScatteringAngle(Element::Te);
      mScatter[53] = new CzyzewskiMottScatteringAngle(Element::I);
      mScatter[54] = new CzyzewskiMottScatteringAngle(Element::Xe);
      mScatter[55] = new CzyzewskiMottScatteringAngle(Element::Cs);
      mScatter[56] = new CzyzewskiMottScatteringAngle(Element::Ba);
      mScatter[57] = new CzyzewskiMottScatteringAngle(Element::La);
      mScatter[58] = new CzyzewskiMottScatteringAngle(Element::Ce);
      mScatter[59] = new CzyzewskiMottScatteringAngle(Element::Pr);
      mScatter[60] = new CzyzewskiMottScatteringAngle(Element::Nd);
      mScatter[61] = new CzyzewskiMottScatteringAngle(Element::Pm);
      mScatter[62] = new CzyzewskiMottScatteringAngle(Element::Sm);
      mScatter[63] = new CzyzewskiMottScatteringAngle(Element::Eu);
      mScatter[64] = new CzyzewskiMottScatteringAngle(Element::Gd);
      mScatter[65] = new CzyzewskiMottScatteringAngle(Element::Tb);
      mScatter[66] = new CzyzewskiMottScatteringAngle(Element::Dy);
      mScatter[67] = new CzyzewskiMottScatteringAngle(Element::Ho);
      mScatter[68] = new CzyzewskiMottScatteringAngle(Element::Er);
      mScatter[69] = new CzyzewskiMottScatteringAngle(Element::Tm);
      mScatter[70] = new CzyzewskiMottScatteringAngle(Element::Yb);
      mScatter[71] = new CzyzewskiMottScatteringAngle(Element::Lu);
      mScatter[72] = new CzyzewskiMottScatteringAngle(Element::Hf);
      mScatter[73] = new CzyzewskiMottScatteringAngle(Element::Ta);
      mScatter[74] = new CzyzewskiMottScatteringAngle(Element::W);
      mScatter[75] = new CzyzewskiMottScatteringAngle(Element::Re);
      mScatter[76] = new CzyzewskiMottScatteringAngle(Element::Os);
      mScatter[77] = new CzyzewskiMottScatteringAngle(Element::Ir);
      mScatter[78] = new CzyzewskiMottScatteringAngle(Element::Pt);
      mScatter[79] = new CzyzewskiMottScatteringAngle(Element::Au);
      mScatter[80] = new CzyzewskiMottScatteringAngle(Element::Hg);
      mScatter[81] = new CzyzewskiMottScatteringAngle(Element::Tl);
      mScatter[82] = new CzyzewskiMottScatteringAngle(Element::Pb);
      mScatter[83] = new CzyzewskiMottScatteringAngle(Element::Bi);
      mScatter[84] = new CzyzewskiMottScatteringAngle(Element::Po);
      mScatter[85] = new CzyzewskiMottScatteringAngle(Element::At);
      mScatter[86] = new CzyzewskiMottScatteringAngle(Element::Rn);
      mScatter[87] = new CzyzewskiMottScatteringAngle(Element::Fr);
      mScatter[88] = new CzyzewskiMottScatteringAngle(Element::Ra);
      mScatter[89] = new CzyzewskiMottScatteringAngle(Element::Ac);
      mScatter[90] = new CzyzewskiMottScatteringAngle(Element::Th);
      mScatter[91] = new CzyzewskiMottScatteringAngle(Element::Pa);
      mScatter[92] = new CzyzewskiMottScatteringAngle(Element::U);
      mScatter[93] = new CzyzewskiMottScatteringAngle(Element::Np);
      mScatter[94] = new CzyzewskiMottScatteringAngle(Element::Pu);
   }

   __device__ CzyzewskiMottScatteringAngle * dScatter[113];

   __global__ void initCuda()
   {
      dScatter[1] = new CzyzewskiMottScatteringAngle(*Element::dH);
      dScatter[2] = new CzyzewskiMottScatteringAngle(*Element::dHe);
      dScatter[3] = new CzyzewskiMottScatteringAngle(*Element::dLi);
      dScatter[4] = new CzyzewskiMottScatteringAngle(*Element::dBe);
      dScatter[5] = new CzyzewskiMottScatteringAngle(*Element::dB);
      dScatter[6] = new CzyzewskiMottScatteringAngle(*Element::dC);
      dScatter[7] = new CzyzewskiMottScatteringAngle(*Element::dN);
      dScatter[8] = new CzyzewskiMottScatteringAngle(*Element::dO);
      dScatter[9] = new CzyzewskiMottScatteringAngle(*Element::dF);
      dScatter[10] = new CzyzewskiMottScatteringAngle(*Element::dNe);
      dScatter[11] = new CzyzewskiMottScatteringAngle(*Element::dNa);
      dScatter[12] = new CzyzewskiMottScatteringAngle(*Element::dMg);
      dScatter[13] = new CzyzewskiMottScatteringAngle(*Element::dAl);
      dScatter[14] = new CzyzewskiMottScatteringAngle(*Element::dSi);
      dScatter[15] = new CzyzewskiMottScatteringAngle(*Element::dP);
      dScatter[16] = new CzyzewskiMottScatteringAngle(*Element::dS);
      dScatter[17] = new CzyzewskiMottScatteringAngle(*Element::dCl);
      dScatter[18] = new CzyzewskiMottScatteringAngle(*Element::dAr);
      dScatter[19] = new CzyzewskiMottScatteringAngle(*Element::dK);
      dScatter[20] = new CzyzewskiMottScatteringAngle(*Element::dCa);
      dScatter[21] = new CzyzewskiMottScatteringAngle(*Element::dSc);
      dScatter[22] = new CzyzewskiMottScatteringAngle(*Element::dTi);
      dScatter[23] = new CzyzewskiMottScatteringAngle(*Element::dV);
      dScatter[24] = new CzyzewskiMottScatteringAngle(*Element::dCr);
      dScatter[25] = new CzyzewskiMottScatteringAngle(*Element::dMn);
      dScatter[26] = new CzyzewskiMottScatteringAngle(*Element::dFe);
      dScatter[27] = new CzyzewskiMottScatteringAngle(*Element::dCo);
      dScatter[28] = new CzyzewskiMottScatteringAngle(*Element::dNi);
      dScatter[29] = new CzyzewskiMottScatteringAngle(*Element::dCu);
      dScatter[30] = new CzyzewskiMottScatteringAngle(*Element::dZn);
      dScatter[31] = new CzyzewskiMottScatteringAngle(*Element::dGa);
      dScatter[32] = new CzyzewskiMottScatteringAngle(*Element::dGe);
      dScatter[33] = new CzyzewskiMottScatteringAngle(*Element::dAs);
      dScatter[34] = new CzyzewskiMottScatteringAngle(*Element::dSe);
      dScatter[35] = new CzyzewskiMottScatteringAngle(*Element::dBr);
      dScatter[36] = new CzyzewskiMottScatteringAngle(*Element::dKr);
      dScatter[37] = new CzyzewskiMottScatteringAngle(*Element::dRb);
      dScatter[38] = new CzyzewskiMottScatteringAngle(*Element::dSr);
      dScatter[39] = new CzyzewskiMottScatteringAngle(*Element::dY);
      dScatter[40] = new CzyzewskiMottScatteringAngle(*Element::dZr);
      dScatter[41] = new CzyzewskiMottScatteringAngle(*Element::dNb);
      dScatter[42] = new CzyzewskiMottScatteringAngle(*Element::dMo);
      dScatter[43] = new CzyzewskiMottScatteringAngle(*Element::dTc);
      dScatter[44] = new CzyzewskiMottScatteringAngle(*Element::dRu);
      dScatter[45] = new CzyzewskiMottScatteringAngle(*Element::dRh);
      dScatter[46] = new CzyzewskiMottScatteringAngle(*Element::dPd);
      dScatter[47] = new CzyzewskiMottScatteringAngle(*Element::dAg);
      dScatter[48] = new CzyzewskiMottScatteringAngle(*Element::dCd);
      dScatter[49] = new CzyzewskiMottScatteringAngle(*Element::dIn);
      dScatter[50] = new CzyzewskiMottScatteringAngle(*Element::dSn);
      dScatter[51] = new CzyzewskiMottScatteringAngle(*Element::dSb);
      dScatter[52] = new CzyzewskiMottScatteringAngle(*Element::dTe);
      dScatter[53] = new CzyzewskiMottScatteringAngle(*Element::dI);
      dScatter[54] = new CzyzewskiMottScatteringAngle(*Element::dXe);
      dScatter[55] = new CzyzewskiMottScatteringAngle(*Element::dCs);
      dScatter[56] = new CzyzewskiMottScatteringAngle(*Element::dBa);
      dScatter[57] = new CzyzewskiMottScatteringAngle(*Element::dLa);
      dScatter[58] = new CzyzewskiMottScatteringAngle(*Element::dCe);
      dScatter[59] = new CzyzewskiMottScatteringAngle(*Element::dPr);
      dScatter[60] = new CzyzewskiMottScatteringAngle(*Element::dNd);
      dScatter[61] = new CzyzewskiMottScatteringAngle(*Element::dPm);
      dScatter[62] = new CzyzewskiMottScatteringAngle(*Element::dSm);
      dScatter[63] = new CzyzewskiMottScatteringAngle(*Element::dEu);
      dScatter[64] = new CzyzewskiMottScatteringAngle(*Element::dGd);
      dScatter[65] = new CzyzewskiMottScatteringAngle(*Element::dTb);
      dScatter[66] = new CzyzewskiMottScatteringAngle(*Element::dDy);
      dScatter[67] = new CzyzewskiMottScatteringAngle(*Element::dHo);
      dScatter[68] = new CzyzewskiMottScatteringAngle(*Element::dEr);
      dScatter[69] = new CzyzewskiMottScatteringAngle(*Element::dTm);
      dScatter[70] = new CzyzewskiMottScatteringAngle(*Element::dYb);
      dScatter[71] = new CzyzewskiMottScatteringAngle(*Element::dLu);
      dScatter[72] = new CzyzewskiMottScatteringAngle(*Element::dHf);
      dScatter[73] = new CzyzewskiMottScatteringAngle(*Element::dTa);
      dScatter[74] = new CzyzewskiMottScatteringAngle(*Element::dW);
      dScatter[75] = new CzyzewskiMottScatteringAngle(*Element::dRe);
      dScatter[76] = new CzyzewskiMottScatteringAngle(*Element::dOs);
      dScatter[77] = new CzyzewskiMottScatteringAngle(*Element::dIr);
      dScatter[78] = new CzyzewskiMottScatteringAngle(*Element::dPt);
      dScatter[79] = new CzyzewskiMottScatteringAngle(*Element::dAu);
      dScatter[80] = new CzyzewskiMottScatteringAngle(*Element::dHg);
      dScatter[81] = new CzyzewskiMottScatteringAngle(*Element::dTl);
      dScatter[82] = new CzyzewskiMottScatteringAngle(*Element::dPb);
      dScatter[83] = new CzyzewskiMottScatteringAngle(*Element::dBi);
      dScatter[84] = new CzyzewskiMottScatteringAngle(*Element::dPo);
      dScatter[85] = new CzyzewskiMottScatteringAngle(*Element::dAt);
      dScatter[86] = new CzyzewskiMottScatteringAngle(*Element::dRn);
      dScatter[87] = new CzyzewskiMottScatteringAngle(*Element::dFr);
      dScatter[88] = new CzyzewskiMottScatteringAngle(*Element::dRa);
      dScatter[89] = new CzyzewskiMottScatteringAngle(*Element::dAc);
      dScatter[90] = new CzyzewskiMottScatteringAngle(*Element::dTh);
      dScatter[91] = new CzyzewskiMottScatteringAngle(*Element::dPa);
      dScatter[92] = new CzyzewskiMottScatteringAngle(*Element::dU);
      dScatter[93] = new CzyzewskiMottScatteringAngle(*Element::dNp);
      dScatter[94] = new CzyzewskiMottScatteringAngle(*Element::dPu);
   }

   __host__ __device__ const VectorXd& CzyzewskiMottScatteringAngle::getMeanFreePath() const
   {
      return mMeanFreePath;
   }

   __host__ __device__ const VectorXd& CzyzewskiMottScatteringAngle::getTotalCrossSection() const
   {
      return mTotalCrossSection;
   }

   __host__ __device__ const MatrixXd& CzyzewskiMottScatteringAngle::getCummulativeDF() const
   {
      return mCummulativeDF;
   }

   __global__ void copyMeanFreePath(const unsigned int i, const double *data, const unsigned int len)
   {
      dScatter[i]->copyMeanFreePath<double>(data, len);
   }

   __global__ void copyTotalCrossSection(const unsigned int i, const double *data, const unsigned int len)
   {
      dScatter[i]->copyTotalCrossSection<double>(data, len);
   }

   __global__ void copyCummulativeDFRow(const unsigned int i, const unsigned int r, const double *data, const unsigned int len)
   {
      dScatter[i]->copyCummulativeDFRow<double>(r, data, len);
   }

   void copyDataToCuda()
   {
      for (int i = 0; i <= 94; ++i) {
         double *dmfp = nullptr;
         const VectorXd& mfp = mScatter[i]->getMeanFreePath();
         checkCudaErrors(cudaMalloc((void **)&dmfp, sizeof(double) * mfp.size()));
         checkCudaErrors(cudaMemcpy(dmfp, mfp.data(), sizeof(double) * mfp.size(), cudaMemcpyHostToDevice));
         copyMeanFreePath << <1, 1 >> >(i, dmfp, mfp.size());
         checkCudaErrors(cudaDeviceSynchronize());
         checkCudaErrors(cudaGetLastError());
         checkCudaErrors(cudaFree(dmfp));

         double *dtcs = nullptr;
         const VectorXd& tcs = mScatter[i]->getTotalCrossSection();
         checkCudaErrors(cudaMalloc((void **)&dtcs, sizeof(double) * tcs.size()));
         checkCudaErrors(cudaMemcpy(dtcs, tcs.data(), sizeof(double) * tcs.size(), cudaMemcpyHostToDevice));
         copyTotalCrossSection << <1, 1 >> >(i, dtcs, tcs.size());
         checkCudaErrors(cudaDeviceSynchronize());
         checkCudaErrors(cudaGetLastError());
         checkCudaErrors(cudaFree(dtcs));

         const MatrixXd& cdf = mScatter[i]->getCummulativeDF();
         for (int r = 0; r <= cdf.size(); ++r) {
            double *dcdf = nullptr;
            const VectorXd& cdf = mScatter[i]->getTotalCrossSection();
            checkCudaErrors(cudaMalloc((void **)&dcdf, sizeof(double) * cdf.size()));
            checkCudaErrors(cudaMemcpy(dcdf, cdf.data(), sizeof(double) * cdf.size(), cudaMemcpyHostToDevice));
            copyCummulativeDFRow << <1, 1 >> >(i, r, dcdf, cdf.size());
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaGetLastError());
            checkCudaErrors(cudaFree(dcdf));
         }
      }
   }

   __host__ __device__ const CzyzewskiMottScatteringAngle& getCMSA(int an)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      return *dScatter[an];
#else
      return *mScatter[an];
#endif
   }

   __host__ __device__ CzyzewskiMottRandomizedScatterFactory::CzyzewskiMottRandomizedScatterFactory() :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      RandomizedScatterFactoryT("Czyzewski Mott cross-section", *Reference::dNullReference)
#else
      RandomizedScatterFactoryT("Czyzewski Mott cross-section", REFERENCE)
#endif
   {
   }

   __host__ __device__ const RandomizedScatterT& CzyzewskiMottRandomizedScatterFactory::get(const Element::Element& elm) const
   {
      return getCMSA(elm.getAtomicNumber());
   }

   void CzyzewskiMottRandomizedScatterFactory::initializeDefaultStrategy()
   {
   }

   const CzyzewskiMottRandomizedScatterFactory CzyzewskiMottFactory;
   const RandomizedScatterFactoryT& Factory = CzyzewskiMottFactory;
}
