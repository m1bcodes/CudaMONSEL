#include "gov\nist\microanalysis\EPQLibrary\CzyzewskiMottScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\CzyzewskiMottCrossSection.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Reference.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ScreenedRutherfordScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

#include "Amphibian\Algorithm.cuh"

namespace CzyzewskiMottScatteringAngle
{
   static const Reference::Author* alRef[] = { &Reference::Czyzewski, &Reference::MacCallum, &Reference::DJoy };
   static const Reference::JournalArticle REFERENCE(Reference::JApplPhys, "68, No. 7", "", 1990, alRef, 3);

   const double MAX_CZYZEWSKI = ToSI::keV(30.0);

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
         double energy = CzyzewskiMottCrossSection::getSpecialEnergy(r);
         mCummulativeDF[r][0] = 0.0;
         for (int c = 1; c < CzyzewskiMottCrossSection::SpecialAngleCount; ++c) {
            double cm = mcs.partialCrossSection(CzyzewskiMottCrossSection::getSpecialAngle(c - 1), energy);
            double cp = mcs.partialCrossSection(CzyzewskiMottCrossSection::getSpecialAngle(c), energy);
            double am = CzyzewskiMottCrossSection::getSpecialAngle(c - 1);
            double ap = CzyzewskiMottCrossSection::getSpecialAngle(c);
            mCummulativeDF[r][c] = mCummulativeDF[r][c - 1] + ((cp + cm) * (ap - am)) / 2.0;
         }
         // Normalize to 1.0
         for (int c = 0; c < CzyzewskiMottCrossSection::SpecialAngleCount; ++c)
            mCummulativeDF[r][c] = mCummulativeDF[r][c] / mCummulativeDF[r][CzyzewskiMottCrossSection::SpecialAngleCount - 1];
      }

   }

   CzyzewskiMottScatteringAngle::CzyzewskiMottScatteringAngle(const ElementT& el) : RandomizedScatterT("Cyzewski", REFERENCE), mElement(el), mRutherford(ScreenedRutherfordScatteringAngle::getSRSA(el.getAtomicNumber()))
   {
      init(el.getAtomicNumber());
   }

   CzyzewskiMottScatteringAngle::CzyzewskiMottScatteringAngle(int an) : RandomizedScatterT("Cyzewski", REFERENCE), mElement(Element::byAtomicNumber(an)), mRutherford(ScreenedRutherfordScatteringAngle::getSRSA(an))
   {
      init(an);
   }

   CzyzewskiMottScatteringAngle::CzyzewskiMottScatteringAngle(const CzyzewskiMottScatteringAngle& other) : RandomizedScatterT("Cyzewski", REFERENCE), mElement(other.mElement), mRutherford(other.mRutherford), mMeanFreePath(other.mMeanFreePath), mTotalCrossSection(other.mTotalCrossSection), mCummulativeDF(other.mCummulativeDF)
   {
   }

   StringT CzyzewskiMottScatteringAngle::toString() const
   {
      return "CrossSection[Czyzewski-Mott," + StringT(mElement.toAbbrev()) + "]";
   }


   const ElementT& CzyzewskiMottScatteringAngle::getElement() const
   {
      return mElement;
   }

   double CzyzewskiMottScatteringAngle::scatteringAngleForSpecialEnergy(int ei, double rand) const
   {
      VectorXd r = mCummulativeDF[ei];
      int ai = Algorithm::binarySearch(r.data(), 0, r.size()-1, rand);
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
      int e = CzyzewskiMottCrossSection::getEnergyIndex(energy);
      double e0 = CzyzewskiMottCrossSection::getSpecialEnergy(e);
      if (!(energy <= e0)) printf("CzyzewskiMottScatteringAngle::randomScatteringAngle 5: %lf\n", energy);
      double a0 = scatteringAngleForSpecialEnergy(e, rand);
      if (energy == e0)
         return a0;
      else {
         double a1 = scatteringAngleForSpecialEnergy(e - 1, rand);
         double e1 = CzyzewskiMottCrossSection::getSpecialEnergy(e - 1);
         if (!(energy > e1)) printf("CzyzewskiMottScatteringAngle::randomScatteringAngle 6: %lf\n", energy);
         if (!(a1 >= 0.0)) printf("CzyzewskiMottScatteringAngle::randomScatteringAngle 7: %lf\n", energy);
         if (!(a0 >= 0.0)) printf("CzyzewskiMottScatteringAngle::randomScatteringAngle 8: %lf\n", energy);
         return a0 + ((energy - e0) / (e1 - e0)) * (a1 - a0);
      }
   }

   double CzyzewskiMottScatteringAngle::randomScatteringAngle(double energy) const
   {
      if (energy < MAX_CZYZEWSKI)
         return randomScatteringAngle(energy, Math2::random());
      else {
         return mRutherford.randomScatteringAngle(energy);
      }
   }

   double CzyzewskiMottScatteringAngle::meanFreePath(double energy) const
   {
      int e = CzyzewskiMottCrossSection::getEnergyIndex(energy);
      if (e == 0)
         e = 1;
      double e0 = CzyzewskiMottCrossSection::getSpecialEnergy(e - 1);
      double e1 = CzyzewskiMottCrossSection::getSpecialEnergy(e);
      if (!(energy >= e0)) printf("CzyzewskiMottScatteringAngle::meanFreePath 1: %lf\n", energy);
      if (!(energy <= e1)) printf("CzyzewskiMottScatteringAngle::meanFreePath 2: %lf\n", energy);
      return mMeanFreePath[e - 1] + ((mMeanFreePath[e] - mMeanFreePath[e - 1]) * ((energy - e0) / (e1 - e0)));
   }

   double CzyzewskiMottScatteringAngle::totalCrossSection(double energy) const
      {
      if (energy < MAX_CZYZEWSKI) {
         int e = CzyzewskiMottCrossSection::getEnergyIndex(energy);
         if (e == 0)
            e = 1;
         double e0 = CzyzewskiMottCrossSection::getSpecialEnergy(e - 1);
         double e1 = CzyzewskiMottCrossSection::getSpecialEnergy(e);
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

   const CzyzewskiMottScatteringAngle CMSA1(1);
   const CzyzewskiMottScatteringAngle CMSA2(2);
   const CzyzewskiMottScatteringAngle CMSA3(3);
   const CzyzewskiMottScatteringAngle CMSA4(4);
   const CzyzewskiMottScatteringAngle CMSA5(5);
   const CzyzewskiMottScatteringAngle CMSA6(6);
   const CzyzewskiMottScatteringAngle CMSA7(7);
   const CzyzewskiMottScatteringAngle CMSA8(8);
   const CzyzewskiMottScatteringAngle CMSA9(9);
   const CzyzewskiMottScatteringAngle CMSA10(10);
   const CzyzewskiMottScatteringAngle CMSA11(11);
   const CzyzewskiMottScatteringAngle CMSA12(12);
   const CzyzewskiMottScatteringAngle CMSA13(13);
   const CzyzewskiMottScatteringAngle CMSA14(14);
   const CzyzewskiMottScatteringAngle CMSA15(15);
   const CzyzewskiMottScatteringAngle CMSA16(16);
   const CzyzewskiMottScatteringAngle CMSA17(17);
   const CzyzewskiMottScatteringAngle CMSA18(18);
   const CzyzewskiMottScatteringAngle CMSA19(19);
   const CzyzewskiMottScatteringAngle CMSA20(20);
   const CzyzewskiMottScatteringAngle CMSA21(21);
   const CzyzewskiMottScatteringAngle CMSA22(22);
   const CzyzewskiMottScatteringAngle CMSA23(23);
   const CzyzewskiMottScatteringAngle CMSA24(24);
   const CzyzewskiMottScatteringAngle CMSA25(25);
   const CzyzewskiMottScatteringAngle CMSA26(26);
   const CzyzewskiMottScatteringAngle CMSA27(27);
   const CzyzewskiMottScatteringAngle CMSA28(28);
   const CzyzewskiMottScatteringAngle CMSA29(29);
   const CzyzewskiMottScatteringAngle CMSA30(30);
   const CzyzewskiMottScatteringAngle CMSA31(31);
   const CzyzewskiMottScatteringAngle CMSA32(32);
   const CzyzewskiMottScatteringAngle CMSA33(33);
   const CzyzewskiMottScatteringAngle CMSA34(34);
   const CzyzewskiMottScatteringAngle CMSA35(35);
   const CzyzewskiMottScatteringAngle CMSA36(36);
   const CzyzewskiMottScatteringAngle CMSA37(37);
   const CzyzewskiMottScatteringAngle CMSA38(38);
   const CzyzewskiMottScatteringAngle CMSA39(39);
   const CzyzewskiMottScatteringAngle CMSA40(40);
   const CzyzewskiMottScatteringAngle CMSA41(41);
   const CzyzewskiMottScatteringAngle CMSA42(42);
   const CzyzewskiMottScatteringAngle CMSA43(43);
   const CzyzewskiMottScatteringAngle CMSA44(44);
   const CzyzewskiMottScatteringAngle CMSA45(45);
   const CzyzewskiMottScatteringAngle CMSA46(46);
   const CzyzewskiMottScatteringAngle CMSA47(47);
   const CzyzewskiMottScatteringAngle CMSA48(48);
   const CzyzewskiMottScatteringAngle CMSA49(49);
   const CzyzewskiMottScatteringAngle CMSA50(50);
   const CzyzewskiMottScatteringAngle CMSA51(51);
   const CzyzewskiMottScatteringAngle CMSA52(52);
   const CzyzewskiMottScatteringAngle CMSA53(53);
   const CzyzewskiMottScatteringAngle CMSA54(54);
   const CzyzewskiMottScatteringAngle CMSA55(55);
   const CzyzewskiMottScatteringAngle CMSA56(56);
   const CzyzewskiMottScatteringAngle CMSA57(57);
   const CzyzewskiMottScatteringAngle CMSA58(58);
   const CzyzewskiMottScatteringAngle CMSA59(59);
   const CzyzewskiMottScatteringAngle CMSA60(60);
   const CzyzewskiMottScatteringAngle CMSA61(61);
   const CzyzewskiMottScatteringAngle CMSA62(62);
   const CzyzewskiMottScatteringAngle CMSA63(63);
   const CzyzewskiMottScatteringAngle CMSA64(64);
   const CzyzewskiMottScatteringAngle CMSA65(65);
   const CzyzewskiMottScatteringAngle CMSA66(66);
   const CzyzewskiMottScatteringAngle CMSA67(67);
   const CzyzewskiMottScatteringAngle CMSA68(68);
   const CzyzewskiMottScatteringAngle CMSA69(69);
   const CzyzewskiMottScatteringAngle CMSA70(70);
   const CzyzewskiMottScatteringAngle CMSA71(71);
   const CzyzewskiMottScatteringAngle CMSA72(72);
   const CzyzewskiMottScatteringAngle CMSA73(73);
   const CzyzewskiMottScatteringAngle CMSA74(74);
   const CzyzewskiMottScatteringAngle CMSA75(75);
   const CzyzewskiMottScatteringAngle CMSA76(76);
   const CzyzewskiMottScatteringAngle CMSA77(77);
   const CzyzewskiMottScatteringAngle CMSA78(78);
   const CzyzewskiMottScatteringAngle CMSA79(79);
   const CzyzewskiMottScatteringAngle CMSA80(80);
   const CzyzewskiMottScatteringAngle CMSA81(81);
   const CzyzewskiMottScatteringAngle CMSA82(82);
   const CzyzewskiMottScatteringAngle CMSA83(83);
   const CzyzewskiMottScatteringAngle CMSA84(84);
   const CzyzewskiMottScatteringAngle CMSA85(85);
   const CzyzewskiMottScatteringAngle CMSA86(86);
   const CzyzewskiMottScatteringAngle CMSA87(87);
   const CzyzewskiMottScatteringAngle CMSA88(88);
   const CzyzewskiMottScatteringAngle CMSA89(89);
   const CzyzewskiMottScatteringAngle CMSA90(90);
   const CzyzewskiMottScatteringAngle CMSA91(91);
   const CzyzewskiMottScatteringAngle CMSA92(92);
   const CzyzewskiMottScatteringAngle CMSA93(93);
   const CzyzewskiMottScatteringAngle CMSA94(94);

   CzyzewskiMottScatteringAngle const * mScatter[113] = {
      nullptr,
      &CMSA1,
      &CMSA2,
      &CMSA3,
      &CMSA4,
      &CMSA5,
      &CMSA6,
      &CMSA7,
      &CMSA8,
      &CMSA9,
      &CMSA10,
      &CMSA11,
      &CMSA12,
      &CMSA13,
      &CMSA14,
      &CMSA15,
      &CMSA16,
      &CMSA17,
      &CMSA18,
      &CMSA19,
      &CMSA20,
      &CMSA21,
      &CMSA22,
      &CMSA23,
      &CMSA24,
      &CMSA25,
      &CMSA26,
      &CMSA27,
      &CMSA28,
      &CMSA29,
      &CMSA30,
      &CMSA31,
      &CMSA32,
      &CMSA33,
      &CMSA34,
      &CMSA35,
      &CMSA36,
      &CMSA37,
      &CMSA38,
      &CMSA39,
      &CMSA40,
      &CMSA41,
      &CMSA42,
      &CMSA43,
      &CMSA44,
      &CMSA45,
      &CMSA46,
      &CMSA47,
      &CMSA48,
      &CMSA49,
      &CMSA50,
      &CMSA51,
      &CMSA52,
      &CMSA53,
      &CMSA54,
      &CMSA55,
      &CMSA56,
      &CMSA57,
      &CMSA58,
      &CMSA59,
      &CMSA60,
      &CMSA61,
      &CMSA62,
      &CMSA63,
      &CMSA64,
      &CMSA65,
      &CMSA66,
      &CMSA67,
      &CMSA68,
      &CMSA69,
      &CMSA70,
      &CMSA71,
      &CMSA72,
      &CMSA73,
      &CMSA74,
      &CMSA75,
      &CMSA76,
      &CMSA77,
      &CMSA78,
      &CMSA79,
      &CMSA80,
      &CMSA81,
      &CMSA82,
      &CMSA83,
      &CMSA84,
      &CMSA85,
      &CMSA86,
      &CMSA87,
      &CMSA88,
      &CMSA89,
      &CMSA90,
      &CMSA91,
      &CMSA92,
      &CMSA93,
      &CMSA94
   };

   const CzyzewskiMottScatteringAngle& getCMSA(int an)
   {
      return *mScatter[an];
   }

   CzyzewskiMottRandomizedScatterFactory::CzyzewskiMottRandomizedScatterFactory() : RandomizedScatterFactoryT("Czyzewski Mott cross-section", REFERENCE)
   {
   }

   const RandomizedScatterT& CzyzewskiMottRandomizedScatterFactory::get(const Element::Element& elm) const
   {
      return getCMSA(elm.getAtomicNumber());
   }

   void CzyzewskiMottRandomizedScatterFactory::initializeDefaultStrategy()
   {
   }

   const CzyzewskiMottRandomizedScatterFactory CzyzewskiMottFactory;
   const RandomizedScatterFactoryT& Factory = CzyzewskiMottFactory;
}
