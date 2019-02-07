#ifndef _ELEMENT_CUH_
#define _ELEMENT_CUH_

#include <stdio.h>
#include <fstream>

#include "../Utility/CSVReader.h"
#include "ToSI.h"

namespace Element
{
   class Element
   {
   public:
      Element(int atomicNo);
      Element();

      static int atomicNumberForName(char* name);
      static Element byName(char* name);
      static Element byAtomicNumber(int an);
      static double getAtomicWeight(int atomicNo);
      static Element const * allElements();
      static Element* range(Element min, Element max);
      static double meanIonizationPotential(int atomicNo);
      int getAtomicNumber();
      double getAtomicWeight();
      double getMass();
      char const * toAbbrev();
      static char const * toAbbrev(int atomicNo);
      static char const * toString(int el);
      static bool isValid(int atomicNo);
      bool isValid();
      int compareTo(Element e);
      int hashCode();
      bool equals(Element el);
      char const * toString();
      double getIonizationEnergy();
      static char * const getListOfAbbreviations(Element minEl, Element maxEl);

   private:
      int mAtomicNumber;

      static void readAtomicWeights();
      Element readResolve();
   };

   const int elmNone = 0;
   const int elmH = 1;
   const int elmHe = 2;
   const int elmLi = 3;
   const int elmBe = 4;
   const int elmB = 5;
   const int elmC = 6;
   const int elmN = 7;
   const int elmO = 8;
   const int elmF = 9;
   const int elmNe = 10;
   const int elmNa = 11;
   const int elmMg = 12;
   const int elmAl = 13;
   const int elmSi = 14;
   const int elmP = 15;
   const int elmS = 16;
   const int elmCl = 17;
   const int elmAr = 18;
   const int elmK = 19;
   const int elmCa = 20;
   const int elmSc = 21;
   const int elmTi = 22;
   const int elmV = 23;
   const int elmCr = 24;
   const int elmMn = 25;
   const int elmFe = 26;
   const int elmCo = 27;
   const int elmNi = 28;
   const int elmCu = 29;
   const int elmZn = 30;
   const int elmGa = 31;
   const int elmGe = 32;
   const int elmAs = 33;
   const int elmSe = 34;
   const int elmBr = 35;
   const int elmKr = 36;
   const int elmRb = 37;
   const int elmSr = 38;
   const int elmY = 39;
   const int elmZr = 40;
   const int elmNb = 41;
   const int elmMo = 42;
   const int elmTc = 43;
   const int elmRu = 44;
   const int elmRh = 45;
   const int elmPd = 46;
   const int elmAg = 47;
   const int elmCd = 48;
   const int elmIn = 49;
   const int elmSn = 50;
   const int elmSb = 51;
   const int elmTe = 52;
   const int elmI = 53;
   const int elmXe = 54;
   const int elmCs = 55;
   const int elmBa = 56;
   const int elmLa = 57;
   const int elmCe = 58;
   const int elmPr = 59;
   const int elmNd = 60;
   const int elmPm = 61;
   const int elmSm = 62;
   const int elmEu = 63;
   const int elmGd = 64;
   const int elmTb = 65;
   const int elmDy = 66;
   const int elmHo = 67;
   const int elmEr = 68;
   const int elmTm = 69;
   const int elmYb = 70;
   const int elmLu = 71;
   const int elmHf = 72;
   const int elmTa = 73;
   const int elmW = 74;
   const int elmRe = 75;
   const int elmOs = 76;
   const int elmIr = 77;
   const int elmPt = 78;
   const int elmAu = 79;
   const int elmHg = 80;
   const int elmTl = 81;
   const int elmPb = 82;
   const int elmBi = 83;
   const int elmPo = 84;
   const int elmAt = 85;
   const int elmRn = 86;
   const int elmFr = 87;
   const int elmRa = 88;
   const int elmAc = 89;
   const int elmTh = 90;
   const int elmPa = 91;
   const int elmU = 92;
   const int elmNp = 93;
   const int elmPu = 94;
   const int elmAm = 95;
   const int elmCm = 96;
   const int elmBk = 97;
   const int elmCf = 98;
   const int elmEs = 99;
   const int elmFm = 100;
   const int elmMd = 101;
   const int elmNo = 102;
   const int elmLr = 103;
   const int elmRf = 104;
   const int elmDb = 105;
   const int elmSg = 106;
   const int elmBh = 107;
   const int elmHs = 108;
   const int elmMt = 109;
   const int elmUun = 110;
   const int elmUuu = 111;
   const int elmUub = 112;
   const int elmEndOfElements = 113;

   const long long serialVersionUID = 0x987360133793L;

   float mIonizationEnergy[104]; // Nominal in Joules, IonizationEnergies.csv
   float mAtomicWeight[112]; // nominal, in AMU, AtomicWeights.csv

   const Element None(0);
   const Element H(1);
   const Element He(2);
   const Element Li(3);
   const Element Be(4);
   const Element B(5);
   const Element C(6);
   const Element N(7);
   const Element O(8);
   const Element F(9);
   const Element Ne(10);
   const Element Na(11);
   const Element Mg(12);
   const Element Al(13);
   const Element Si(14);
   const Element P(15);
   const Element S(16);
   const Element Cl(17);
   const Element Ar(18);
   const Element K(19);
   const Element Ca(20);
   const Element Sc(21);
   const Element Ti(22);
   const Element V(23);
   const Element Cr(24);
   const Element Mn(25);
   const Element Fe(26);
   const Element Co(27);
   const Element Ni(28);
   const Element Cu(29);
   const Element Zn(30);
   const Element Ga(31);
   const Element Ge(32);
   const Element As(33);
   const Element Se(34);
   const Element Br(35);
   const Element Kr(36);
   const Element Rb(37);
   const Element Sr(38);
   const Element Y(39);
   const Element Zr(40);
   const Element Nb(41);
   const Element Mo(42);
   const Element Tc(43);
   const Element Ru(44);
   const Element Rh(45);
   const Element Pd(46);
   const Element Ag(47);
   const Element Cd(48);
   const Element In(49);
   const Element Sn(50);
   const Element Sb(51);
   const Element Te(52);
   const Element I(53);
   const Element Xe(54);
   const Element Cs(55);
   const Element Ba(56);
   const Element La(57);
   const Element Ce(58);
   const Element Pr(59);
   const Element Nd(60);
   const Element Pm(61);
   const Element Sm(62);
   const Element Eu(63);
   const Element Gd(64);
   const Element Tb(65);
   const Element Dy(66);
   const Element Ho(67);
   const Element Er(68);
   const Element Tm(69);
   const Element Yb(70);
   const Element Lu(71);
   const Element Hf(72);
   const Element Ta(73);
   const Element W(74);
   const Element Re(75);
   const Element Os(76);
   const Element Ir(77);
   const Element Pt(78);
   const Element Au(79);
   const Element Hg(80);
   const Element Tl(81);
   const Element Pb(82);
   const Element Bi(83);
   const Element Po(84);
   const Element At(85);
   const Element Rn(86);
   const Element Fr(87);
   const Element Ra(88);
   const Element Ac(89);
   const Element Th(90);
   const Element Pa(91);
   const Element U(92);
   const Element Np(93);
   const Element Pu(94);
   const Element Am(95);
   const Element Cm(96);
   const Element Bk(97);
   const Element Cf(98);
   const Element Es(99);
   const Element Fm(100);
   const Element Md(101);
   const Element No(102);
   const Element Lr(103);
   const Element Rf(104);
   const Element Db(105);
   const Element Sg(106);
   const Element Bh(107);
   const Element Hs(108);
   const Element Mt(109);
   const Element Uun(110);
   const Element Uuu(111);
   const Element Uub(112);

   const Element mAllElements[] = {
      H,
      He,
      Li,
      Be,
      B,
      C,
      N,
      O,
      F,
      Ne,
      Na,
      Mg,
      Al,
      Si,
      P,
      S,
      Cl,
      Ar,
      K,
      Ca,
      Sc,
      Ti,
      V,
      Cr,
      Mn,
      Fe,
      Co,
      Ni,
      Cu,
      Zn,
      Ga,
      Ge,
      As,
      Se,
      Br,
      Kr,
      Rb,
      Sr,
      Y,
      Zr,
      Nb,
      Mo,
      Tc,
      Ru,
      Rh,
      Pd,
      Ag,
      Cd,
      In,
      Sn,
      Sb,
      Te,
      I,
      Xe,
      Cs,
      Ba,
      La,
      Ce,
      Pr,
      Nd,
      Pm,
      Sm,
      Eu,
      Gd,
      Tb,
      Dy,
      Ho,
      Er,
      Tm,
      Yb,
      Lu,
      Hf,
      Ta,
      W,
      Re,
      Os,
      Ir,
      Pt,
      Au,
      Hg,
      Tl,
      Pb,
      Bi,
      Po,
      At,
      Rn,
      Fr,
      Ra,
      Ac,
      Th,
      Pa,
      U,
      Np,
      Pu,
      Am,
      Cm,
      Bk,
      Cf,
      Es,
      Fm,
      Md,
      No,
      Lr,
      Rf,
      Db,
      Sg,
      Bh,
      Hs,
      Mt,
      Uun,
      Uuu,
      Uub
   };

   char const * const mElementNames[] = {
      "None",
      "Hydrogen",
      "Helium",
      "Lithium",
      "Beryllium",
      "Boron",
      "Carbon",
      "Nitrogen",
      "Oxygen",
      "Fluorine",
      "Neon",
      "Sodium",
      "Magnesium",
      "Aluminum",
      "Silicon",
      "Phosphorus",
      "Sulfur",
      "Chlorine",
      "Argon",
      "Potassium",
      "Calcium",
      "Scandium",
      "Titanium",
      "Vanadium",
      "Chromium",
      "Manganese",
      "Iron",
      "Cobalt",
      "Nickel",
      "Copper",
      "Zinc",
      "Gallium",
      "Germanium",
      "Arsenic",
      "Selenium",
      "Bromine",
      "Krypton",
      "Rubidium",
      "Strontium",
      "Yttrium",
      "Zirconium",
      "Niobium",
      "Molybdenum",
      "Technetium",
      "Ruthenium",
      "Rhodium",
      "Palladium",
      "Silver",
      "Cadmium",
      "Indium",
      "Tin",
      "Antimony",
      "Tellurium",
      "Iodine",
      "Xenon",
      "Cesium",
      "Barium",
      "Lanthanum",
      "Cerium",
      "Praseodymium",
      "Neodymium",
      "Promethium",
      "Samarium",
      "Europium",
      "Gadolinium",
      "Terbium",
      "Dysprosium",
      "Holmium",
      "Erbium",
      "Thulium",
      "Ytterbium",
      "Lutetium",
      "Hafnium",
      "Tantalum",
      "Tungsten",
      "Rhenium",
      "Osmium",
      "Iridium",
      "Platinum",
      "Gold",
      "Mercury",
      "Thallium",
      "Lead",
      "Bismuth",
      "Polonium",
      "Astatine",
      "Radon",
      "Francium",
      "Radium",
      "Actinium",
      "Thorium",
      "Protactinium",
      "Uranium",
      "Neptunium",
      "Plutonium",
      "Americium",
      "Curium",
      "Berkelium",
      "Californium",
      "Einsteinium",
      "Fermium",
      "Mendelevium",
      "Nobelium",
      "Lawrencium",
      "Rutherfordium",
      "Dubnium",
      "Seaborgium",
      "Bohrium",
      "Hassium",
      "Meitnerium",
      "Ununnilium",
      "Unununium",
      "Ununbium",
      "End-of-elements"
   };

   char const * const mAbbreviations[] = {
      "",
      "H",
      "He",
      "Li",
      "Be",
      "B",
      "C",
      "N",
      "O",
      "F",
      "Ne",
      "Na",
      "Mg",
      "Al",
      "Si",
      "P",
      "S",
      "Cl",
      "Ar",
      "K",
      "Ca",
      "Sc",
      "Ti",
      "V",
      "Cr",
      "Mn",
      "Fe",
      "Co",
      "Ni",
      "Cu",
      "Zn",
      "Ga",
      "Ge",
      "As",
      "Se",
      "Br",
      "Kr",
      "Rb",
      "Sr",
      "Y",
      "Zr",
      "Nb",
      "Mo",
      "Tc",
      "Ru",
      "Rh",
      "Pd",
      "Ag",
      "Cd",
      "In",
      "Sn",
      "Sb",
      "Te",
      "I",
      "Xe",
      "Cs",
      "Ba",
      "La",
      "Ce",
      "Pr",
      "Nd",
      "Pm",
      "Sm",
      "Eu",
      "Gd",
      "Tb",
      "Dy",
      "Ho",
      "Er",
      "Tm",
      "Yb",
      "Lu",
      "Hf",
      "Ta",
      "W",
      "Re",
      "Os",
      "Ir",
      "Pt",
      "Au",
      "Hg",
      "Tl",
      "Pb",
      "Bi",
      "Po",
      "At",
      "Rn",
      "Fr",
      "Ra",
      "Ac",
      "Th",
      "Pa",
      "U",
      "Np",
      "Pu",
      "Am",
      "Cm",
      "Bk",
      "Cf",
      "Es",
      "Fm",
      "Md",
      "No",
      "Lr",
      "Rf",
      "Db",
      "Sg",
      "Bh",
      "Hs",
      "Mt",
      "Uun",
      "Uuu",
      "Uub",
      "EOE"
   };
}
#endif
