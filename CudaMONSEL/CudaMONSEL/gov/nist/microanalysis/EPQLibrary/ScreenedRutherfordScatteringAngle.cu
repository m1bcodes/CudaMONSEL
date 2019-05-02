#include "gov\nist\microanalysis\EPQLibrary\ScreenedRutherfordScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

namespace ScreenedRutherfordScatteringAngle
{
   static Reference::CrudeReference REFERENCE(Reference::CrudeReference("NBSMONTE.FOR"));

   ScreenedRutherfordScatteringAngle::ScreenedRutherfordScatteringAngle(const ElementT& elm) : RandomizedScatterT("Screened Rutherford", REFERENCE), mElement(elm)
   {
   }

   ScreenedRutherfordScatteringAngle::ScreenedRutherfordScatteringAngle(int an) : RandomizedScatterT("Screened Rutherford", REFERENCE), mElement(Element::byAtomicNumber(an))
   {
   }
   
   ScreenedRutherfordScatteringAngle::ScreenedRutherfordScatteringAngle(const ScreenedRutherfordScatteringAngle& other) : RandomizedScatterT("Screened Rutherford", REFERENCE), mElement(other.mElement)
   {
   }

   StringT ScreenedRutherfordScatteringAngle::toString() const
   {
      return "CrossSection[Screened-Rutherford," + StringT(mElement.toAbbrev()) + "]";
   }

   
   const ElementT& ScreenedRutherfordScatteringAngle::getElement() const
   {
      return mElement;
   }

   double ScreenedRutherfordScatteringAngle::totalCrossSection(double energy) const
   {
      // Ref: Heinrich 1981 p 459 convert to SI units
      double z = mElement.getAtomicNumber();
      double zp = ::pow(z, 1.0 / 3.0);
      return (7.670843088080456e-38 * zp * (1.0 + z)) / ((energy + ((5.44967975966321e-19 * zp * zp))));
   }

   double ScreenedRutherfordScatteringAngle::randomScatteringAngle(double energy) const
   {
      // This method for calculating the scattering angle is taken from
      // NBSMONTE.FOR
      double alpha = (5.44968e-19 * ::pow(mElement.getAtomicNumber(), 2.0 / 3.0)) / energy;
      double r = Math2::random();
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

   const ScreenedRutherfordScatteringAngle SRSA1(1);
   const ScreenedRutherfordScatteringAngle SRSA2(2);
   const ScreenedRutherfordScatteringAngle SRSA3(3);
   const ScreenedRutherfordScatteringAngle SRSA4(4);
   const ScreenedRutherfordScatteringAngle SRSA5(5);
   const ScreenedRutherfordScatteringAngle SRSA6(6);
   const ScreenedRutherfordScatteringAngle SRSA7(7);
   const ScreenedRutherfordScatteringAngle SRSA8(8);
   const ScreenedRutherfordScatteringAngle SRSA9(9);
   const ScreenedRutherfordScatteringAngle SRSA10(10);
   const ScreenedRutherfordScatteringAngle SRSA11(11);
   const ScreenedRutherfordScatteringAngle SRSA12(12);
   const ScreenedRutherfordScatteringAngle SRSA13(13);
   const ScreenedRutherfordScatteringAngle SRSA14(14);
   const ScreenedRutherfordScatteringAngle SRSA15(15);
   const ScreenedRutherfordScatteringAngle SRSA16(16);
   const ScreenedRutherfordScatteringAngle SRSA17(17);
   const ScreenedRutherfordScatteringAngle SRSA18(18);
   const ScreenedRutherfordScatteringAngle SRSA19(19);
   const ScreenedRutherfordScatteringAngle SRSA20(20);
   const ScreenedRutherfordScatteringAngle SRSA21(21);
   const ScreenedRutherfordScatteringAngle SRSA22(22);
   const ScreenedRutherfordScatteringAngle SRSA23(23);
   const ScreenedRutherfordScatteringAngle SRSA24(24);
   const ScreenedRutherfordScatteringAngle SRSA25(25);
   const ScreenedRutherfordScatteringAngle SRSA26(26);
   const ScreenedRutherfordScatteringAngle SRSA27(27);
   const ScreenedRutherfordScatteringAngle SRSA28(28);
   const ScreenedRutherfordScatteringAngle SRSA29(29);
   const ScreenedRutherfordScatteringAngle SRSA30(30);
   const ScreenedRutherfordScatteringAngle SRSA31(31);
   const ScreenedRutherfordScatteringAngle SRSA32(32);
   const ScreenedRutherfordScatteringAngle SRSA33(33);
   const ScreenedRutherfordScatteringAngle SRSA34(34);
   const ScreenedRutherfordScatteringAngle SRSA35(35);
   const ScreenedRutherfordScatteringAngle SRSA36(36);
   const ScreenedRutherfordScatteringAngle SRSA37(37);
   const ScreenedRutherfordScatteringAngle SRSA38(38);
   const ScreenedRutherfordScatteringAngle SRSA39(39);
   const ScreenedRutherfordScatteringAngle SRSA40(40);
   const ScreenedRutherfordScatteringAngle SRSA41(41);
   const ScreenedRutherfordScatteringAngle SRSA42(42);
   const ScreenedRutherfordScatteringAngle SRSA43(43);
   const ScreenedRutherfordScatteringAngle SRSA44(44);
   const ScreenedRutherfordScatteringAngle SRSA45(45);
   const ScreenedRutherfordScatteringAngle SRSA46(46);
   const ScreenedRutherfordScatteringAngle SRSA47(47);
   const ScreenedRutherfordScatteringAngle SRSA48(48);
   const ScreenedRutherfordScatteringAngle SRSA49(49);
   const ScreenedRutherfordScatteringAngle SRSA50(50);
   const ScreenedRutherfordScatteringAngle SRSA51(51);
   const ScreenedRutherfordScatteringAngle SRSA52(52);
   const ScreenedRutherfordScatteringAngle SRSA53(53);
   const ScreenedRutherfordScatteringAngle SRSA54(54);
   const ScreenedRutherfordScatteringAngle SRSA55(55);
   const ScreenedRutherfordScatteringAngle SRSA56(56);
   const ScreenedRutherfordScatteringAngle SRSA57(57);
   const ScreenedRutherfordScatteringAngle SRSA58(58);
   const ScreenedRutherfordScatteringAngle SRSA59(59);
   const ScreenedRutherfordScatteringAngle SRSA60(60);
   const ScreenedRutherfordScatteringAngle SRSA61(61);
   const ScreenedRutherfordScatteringAngle SRSA62(62);
   const ScreenedRutherfordScatteringAngle SRSA63(63);
   const ScreenedRutherfordScatteringAngle SRSA64(64);
   const ScreenedRutherfordScatteringAngle SRSA65(65);
   const ScreenedRutherfordScatteringAngle SRSA66(66);
   const ScreenedRutherfordScatteringAngle SRSA67(67);
   const ScreenedRutherfordScatteringAngle SRSA68(68);
   const ScreenedRutherfordScatteringAngle SRSA69(69);
   const ScreenedRutherfordScatteringAngle SRSA70(70);
   const ScreenedRutherfordScatteringAngle SRSA71(71);
   const ScreenedRutherfordScatteringAngle SRSA72(72);
   const ScreenedRutherfordScatteringAngle SRSA73(73);
   const ScreenedRutherfordScatteringAngle SRSA74(74);
   const ScreenedRutherfordScatteringAngle SRSA75(75);
   const ScreenedRutherfordScatteringAngle SRSA76(76);
   const ScreenedRutherfordScatteringAngle SRSA77(77);
   const ScreenedRutherfordScatteringAngle SRSA78(78);
   const ScreenedRutherfordScatteringAngle SRSA79(79);
   const ScreenedRutherfordScatteringAngle SRSA80(80);
   const ScreenedRutherfordScatteringAngle SRSA81(81);
   const ScreenedRutherfordScatteringAngle SRSA82(82);
   const ScreenedRutherfordScatteringAngle SRSA83(83);
   const ScreenedRutherfordScatteringAngle SRSA84(84);
   const ScreenedRutherfordScatteringAngle SRSA85(85);
   const ScreenedRutherfordScatteringAngle SRSA86(86);
   const ScreenedRutherfordScatteringAngle SRSA87(87);
   const ScreenedRutherfordScatteringAngle SRSA88(88);
   const ScreenedRutherfordScatteringAngle SRSA89(89);
   const ScreenedRutherfordScatteringAngle SRSA90(90);
   const ScreenedRutherfordScatteringAngle SRSA91(91);
   const ScreenedRutherfordScatteringAngle SRSA92(92);
   const ScreenedRutherfordScatteringAngle SRSA93(93);
   const ScreenedRutherfordScatteringAngle SRSA94(94);
   const ScreenedRutherfordScatteringAngle SRSA95(95);
   const ScreenedRutherfordScatteringAngle SRSA96(96);

   ScreenedRutherfordScatteringAngle const * mScatter[113] = {
      &SRSA1,
      &SRSA2,
      &SRSA3,
      &SRSA4,
      &SRSA5,
      &SRSA6,
      &SRSA7,
      &SRSA8,
      &SRSA9,
      &SRSA10,
      &SRSA11,
      &SRSA12,
      &SRSA13,
      &SRSA14,
      &SRSA15,
      &SRSA16,
      &SRSA17,
      &SRSA18,
      &SRSA19,
      &SRSA20,
      &SRSA21,
      &SRSA22,
      &SRSA23,
      &SRSA24,
      &SRSA25,
      &SRSA26,
      &SRSA27,
      &SRSA28,
      &SRSA29,
      &SRSA30,
      &SRSA31,
      &SRSA32,
      &SRSA33,
      &SRSA34,
      &SRSA35,
      &SRSA36,
      &SRSA37,
      &SRSA38,
      &SRSA39,
      &SRSA40,
      &SRSA41,
      &SRSA42,
      &SRSA43,
      &SRSA44,
      &SRSA45,
      &SRSA46,
      &SRSA47,
      &SRSA48,
      &SRSA49,
      &SRSA50,
      &SRSA51,
      &SRSA52,
      &SRSA53,
      &SRSA54,
      &SRSA55,
      &SRSA56,
      &SRSA57,
      &SRSA58,
      &SRSA59,
      &SRSA60,
      &SRSA61,
      &SRSA62,
      &SRSA63,
      &SRSA64,
      &SRSA65,
      &SRSA66,
      &SRSA67,
      &SRSA68,
      &SRSA69,
      &SRSA70,
      &SRSA71,
      &SRSA72,
      &SRSA73,
      &SRSA74,
      &SRSA75,
      &SRSA76,
      &SRSA77,
      &SRSA78,
      &SRSA79,
      &SRSA80,
      &SRSA81,
      &SRSA82,
      &SRSA83,
      &SRSA84,
      &SRSA85,
      &SRSA86,
      &SRSA87,
      &SRSA88,
      &SRSA89,
      &SRSA90,
      &SRSA91,
      &SRSA92,
      &SRSA93,
      &SRSA94,
      &SRSA95,
      &SRSA96
   };

   ScreenedRutherfordRandomizedScatterFactory::ScreenedRutherfordRandomizedScatterFactory() : RandomizedScatterFactoryT("Screened Rutherford elastic cross-section", REFERENCE)
   {
   }

   const RandomizedScatterT& ScreenedRutherfordRandomizedScatterFactory::get(const ElementT& elm) const
   {
      return *mScatter[elm.getAtomicNumber()];
   }

   void ScreenedRutherfordRandomizedScatterFactory::initializeDefaultStrategy()
   {
   }

   const ScreenedRutherfordRandomizedScatterFactory FactoryScreenedRutherford;
   const RandomizedScatterFactoryT& Factory = FactoryScreenedRutherford;
}
