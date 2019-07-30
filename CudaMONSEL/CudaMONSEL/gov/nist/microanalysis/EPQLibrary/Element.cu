#include "Element.cuh"

#include <algorithm>

#include "CudaUtil.h"

namespace Element
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __device__ static const int numIonizationEnergy = 104;
   __device__ static const int numAtomicWeight = 112;

   __device__ const long long serialVersionUID = 0x987360133793L;

   __device__ const int elmNone = 0;
   __device__ const int elmH = 1;
   __device__ const int elmHe = 2;
   __device__ const int elmLi = 3;
   __device__ const int elmBe = 4;
   __device__ const int elmB = 5;
   __device__ const int elmC = 6;
   __device__ const int elmN = 7;
   __device__ const int elmO = 8;
   __device__ const int elmF = 9;
   __device__ const int elmNe = 10;
   __device__ const int elmNa = 11;
   __device__ const int elmMg = 12;
   __device__ const int elmAl = 13;
   __device__ const int elmSi = 14;
   __device__ const int elmP = 15;
   __device__ const int elmS = 16;
   __device__ const int elmCl = 17;
   __device__ const int elmAr = 18;
   __device__ const int elmK = 19;
   __device__ const int elmCa = 20;
   __device__ const int elmSc = 21;
   __device__ const int elmTi = 22;
   __device__ const int elmV = 23;
   __device__ const int elmCr = 24;
   __device__ const int elmMn = 25;
   __device__ const int elmFe = 26;
   __device__ const int elmCo = 27;
   __device__ const int elmNi = 28;
   __device__ const int elmCu = 29;
   __device__ const int elmZn = 30;
   __device__ const int elmGa = 31;
   __device__ const int elmGe = 32;
   __device__ const int elmAs = 33;
   __device__ const int elmSe = 34;
   __device__ const int elmBr = 35;
   __device__ const int elmKr = 36;
   __device__ const int elmRb = 37;
   __device__ const int elmSr = 38;
   __device__ const int elmY = 39;
   __device__ const int elmZr = 40;
   __device__ const int elmNb = 41;
   __device__ const int elmMo = 42;
   __device__ const int elmTc = 43;
   __device__ const int elmRu = 44;
   __device__ const int elmRh = 45;
   __device__ const int elmPd = 46;
   __device__ const int elmAg = 47;
   __device__ const int elmCd = 48;
   __device__ const int elmIn = 49;
   __device__ const int elmSn = 50;
   __device__ const int elmSb = 51;
   __device__ const int elmTe = 52;
   __device__ const int elmI = 53;
   __device__ const int elmXe = 54;
   __device__ const int elmCs = 55;
   __device__ const int elmBa = 56;
   __device__ const int elmLa = 57;
   __device__ const int elmCe = 58;
   __device__ const int elmPr = 59;
   __device__ const int elmNd = 60;
   __device__ const int elmPm = 61;
   __device__ const int elmSm = 62;
   __device__ const int elmEu = 63;
   __device__ const int elmGd = 64;
   __device__ const int elmTb = 65;
   __device__ const int elmDy = 66;
   __device__ const int elmHo = 67;
   __device__ const int elmEr = 68;
   __device__ const int elmTm = 69;
   __device__ const int elmYb = 70;
   __device__ const int elmLu = 71;
   __device__ const int elmHf = 72;
   __device__ const int elmTa = 73;
   __device__ const int elmW = 74;
   __device__ const int elmRe = 75;
   __device__ const int elmOs = 76;
   __device__ const int elmIr = 77;
   __device__ const int elmPt = 78;
   __device__ const int elmAu = 79;
   __device__ const int elmHg = 80;
   __device__ const int elmTl = 81;
   __device__ const int elmPb = 82;
   __device__ const int elmBi = 83;
   __device__ const int elmPo = 84;
   __device__ const int elmAt = 85;
   __device__ const int elmRn = 86;
   __device__ const int elmFr = 87;
   __device__ const int elmRa = 88;
   __device__ const int elmAc = 89;
   __device__ const int elmTh = 90;
   __device__ const int elmPa = 91;
   __device__ const int elmU = 92;
   __device__ const int elmNp = 93;
   __device__ const int elmPu = 94;
   __device__ const int elmAm = 95;
   __device__ const int elmCm = 96;
   __device__ const int elmBk = 97;
   __device__ const int elmCf = 98;
   __device__ const int elmEs = 99;
   __device__ const int elmFm = 100;
   __device__ const int elmMd = 101;
   __device__ const int elmNo = 102;
   __device__ const int elmLr = 103;
   __device__ const int elmRf = 104;
   __device__ const int elmDb = 105;
   __device__ const int elmSg = 106;
   __device__ const int elmBh = 107;
   __device__ const int elmHs = 108;
   __device__ const int elmMt = 109;
   __device__ const int elmUun = 110;
   __device__ const int elmUuu = 111;
   __device__ const int elmUub = 112;
   __device__ const int elmEndOfElements = 113;
#else
   static const int numIonizationEnergy = 104;
   static const int numAtomicWeight = 112;

   const long long serialVersionUID = 0x987360133793L;

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
#endif

   __device__ const Element *dNone = nullptr;
   __device__ const Element *dH = nullptr;
   __device__ const Element *dHe = nullptr;
   __device__ const Element *dLi = nullptr;
   __device__ const Element *dBe = nullptr;
   __device__ const Element *dB = nullptr;
   __device__ const Element *dC = nullptr;
   __device__ const Element *dN = nullptr;
   __device__ const Element *dO = nullptr;
   __device__ const Element *dF = nullptr;
   __device__ const Element *dNe = nullptr;
   __device__ const Element *dNa = nullptr;
   __device__ const Element *dMg = nullptr;
   __device__ const Element *dAl = nullptr;
   __device__ const Element *dSi = nullptr;
   __device__ const Element *dP = nullptr;
   __device__ const Element *dS = nullptr;
   __device__ const Element *dCl = nullptr;
   __device__ const Element *dAr = nullptr;
   __device__ const Element *dK = nullptr;
   __device__ const Element *dCa = nullptr;
   __device__ const Element *dSc = nullptr;
   __device__ const Element *dTi = nullptr;
   __device__ const Element *dV = nullptr;
   __device__ const Element *dCr = nullptr;
   __device__ const Element *dMn = nullptr;
   __device__ const Element *dFe = nullptr;
   __device__ const Element *dCo = nullptr;
   __device__ const Element *dNi = nullptr;
   __device__ const Element *dCu = nullptr;
   __device__ const Element *dZn = nullptr;
   __device__ const Element *dGa = nullptr;
   __device__ const Element *dGe = nullptr;
   __device__ const Element *dAs = nullptr;
   __device__ const Element *dSe = nullptr;
   __device__ const Element *dBr = nullptr;
   __device__ const Element *dKr = nullptr;
   __device__ const Element *dRb = nullptr;
   __device__ const Element *dSr = nullptr;
   __device__ const Element *dY = nullptr;
   __device__ const Element *dZr = nullptr;
   __device__ const Element *dNb = nullptr;
   __device__ const Element *dMo = nullptr;
   __device__ const Element *dTc = nullptr;
   __device__ const Element *dRu = nullptr;
   __device__ const Element *dRh = nullptr;
   __device__ const Element *dPd = nullptr;
   __device__ const Element *dAg = nullptr;
   __device__ const Element *dCd = nullptr;
   __device__ const Element *dIn = nullptr;
   __device__ const Element *dSn = nullptr;
   __device__ const Element *dSb = nullptr;
   __device__ const Element *dTe = nullptr;
   __device__ const Element *dI = nullptr;
   __device__ const Element *dXe = nullptr;
   __device__ const Element *dCs = nullptr;
   __device__ const Element *dBa = nullptr;
   __device__ const Element *dLa = nullptr;
   __device__ const Element *dCe = nullptr;
   __device__ const Element *dPr = nullptr;
   __device__ const Element *dNd = nullptr;
   __device__ const Element *dPm = nullptr;
   __device__ const Element *dSm = nullptr;
   __device__ const Element *dEu = nullptr;
   __device__ const Element *dGd = nullptr;
   __device__ const Element *dTb = nullptr;
   __device__ const Element *dDy = nullptr;
   __device__ const Element *dHo = nullptr;
   __device__ const Element *dEr = nullptr;
   __device__ const Element *dTm = nullptr;
   __device__ const Element *dYb = nullptr;
   __device__ const Element *dLu = nullptr;
   __device__ const Element *dHf = nullptr;
   __device__ const Element *dTa = nullptr;
   __device__ const Element *dW = nullptr;
   __device__ const Element *dRe = nullptr;
   __device__ const Element *dOs = nullptr;
   __device__ const Element *dIr = nullptr;
   __device__ const Element *dPt = nullptr;
   __device__ const Element *dAu = nullptr;
   __device__ const Element *dHg = nullptr;
   __device__ const Element *dTl = nullptr;
   __device__ const Element *dPb = nullptr;
   __device__ const Element *dBi = nullptr;
   __device__ const Element *dPo = nullptr;
   __device__ const Element *dAt = nullptr;
   __device__ const Element *dRn = nullptr;
   __device__ const Element *dFr = nullptr;
   __device__ const Element *dRa = nullptr;
   __device__ const Element *dAc = nullptr;
   __device__ const Element *dTh = nullptr;
   __device__ const Element *dPa = nullptr;
   __device__ const Element *dU = nullptr;
   __device__ const Element *dNp = nullptr;
   __device__ const Element *dPu = nullptr;
   __device__ const Element *dAm = nullptr;
   __device__ const Element *dCm = nullptr;
   __device__ const Element *dBk = nullptr;
   __device__ const Element *dCf = nullptr;
   __device__ const Element *dEs = nullptr;
   __device__ const Element *dFm = nullptr;
   __device__ const Element *dMd = nullptr;
   __device__ const Element *dNo = nullptr;
   __device__ const Element *dLr = nullptr;
   __device__ const Element *dRf = nullptr;
   __device__ const Element *dDb = nullptr;
   __device__ const Element *dSg = nullptr;
   __device__ const Element *dBh = nullptr;
   __device__ const Element *dHs = nullptr;
   __device__ const Element *dMt = nullptr;
   __device__ const Element *dUun = nullptr;
   __device__ const Element *dUuu = nullptr;
   __device__ const Element *dUub = nullptr;

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

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ char const * const mAbbreviations[] = {
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
#else
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
#endif

   Element const * mAllElements[numAtomicWeight] = {
      &H,
      &He,
      &Li,
      &Be,
      &B,
      &C,
      &N,
      &O,
      &F,
      &Ne,
      &Na,
      &Mg,
      &Al,
      &Si,
      &P,
      &S,
      &Cl,
      &Ar,
      &K,
      &Ca,
      &Sc,
      &Ti,
      &V,
      &Cr,
      &Mn,
      &Fe,
      &Co,
      &Ni,
      &Cu,
      &Zn,
      &Ga,
      &Ge,
      &As,
      &Se,
      &Br,
      &Kr,
      &Rb,
      &Sr,
      &Y,
      &Zr,
      &Nb,
      &Mo,
      &Tc,
      &Ru,
      &Rh,
      &Pd,
      &Ag,
      &Cd,
      &In,
      &Sn,
      &Sb,
      &Te,
      &I,
      &Xe,
      &Cs,
      &Ba,
      &La,
      &Ce,
      &Pr,
      &Nd,
      &Pm,
      &Sm,
      &Eu,
      &Gd,
      &Tb,
      &Dy,
      &Ho,
      &Er,
      &Tm,
      &Yb,
      &Lu,
      &Hf,
      &Ta,
      &W,
      &Re,
      &Os,
      &Ir,
      &Pt,
      &Au,
      &Hg,
      &Tl,
      &Pb,
      &Bi,
      &Po,
      &At,
      &Rn,
      &Fr,
      &Ra,
      &Ac,
      &Th,
      &Pa,
      &U,
      &Np,
      &Pu,
      &Am,
      &Cm,
      &Bk,
      &Cf,
      &Es,
      &Fm,
      &Md,
      &No,
      &Lr,
      &Rf,
      &Db,
      &Sg,
      &Bh,
      &Hs,
      &Mt,
      &Uun,
      &Uuu,
      &Uub
   };

   static std::vector<float> mAtomicWeight; // nominal, in AMU, AtomicWeights.csv
   static std::vector<float> mIonizationEnergy; // Nominal in Joules, IonizationEnergies.csv

   __constant__ static float d_AtomicWeight[112], d_IonizationEnergy[104];

   static void readAtomicWeights()
   {
      if (!mAtomicWeight.empty()) return;
      mAtomicWeight.resize(numAtomicWeight);

      try {
         char filepath[] = ".\\gov\\nist\\microanalysis\\EPQLibrary\\AtomicWeights.csv";
         std::ifstream file(filepath);
         printf("Reading: %s\n", filepath);
         int idx = 0;
         for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
            if ((*loop)[0][0] == '/' && (*loop)[0][1] == '/') { // check if the first line should be removed
               continue;
            }
            mAtomicWeight[idx] = std::stof((*loop)[0]);
            ++idx;
         }
         file.close();
      }
      catch (std::exception&) {
         printf("didnt see file AtomicWeights.csv\n");
         throw 0; //throw new EPQFatalException("Fatal error while attempting to load the atomic weights data file.");
      }
   }

   static void readIonizationEnergy()
   {
      if (!mIonizationEnergy.empty()) return;
      mIonizationEnergy.resize(numIonizationEnergy);

      char filepath[] = ".\\gov\\nist\\microanalysis\\EPQLibrary\\IonizationEnergies.csv";
      //float hIonizationEnergy[numIonizationEnergy];
      try {
         printf("Reading: %s\n", filepath);
         std::ifstream file(filepath);
         if (!file.good()) throw 0;
         int idx = 0;
         for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
            if ((*loop)[0][0] == '/' && (*loop)[0][1] == '/') { // check if the first line should be removed
               continue;
            }
            if (CSVIterator::IsNaN((*loop)[0])) {
               mIonizationEnergy[idx] = -1.0;
            }
            else {
               mIonizationEnergy[idx] = ToSI::eV(std::stof((*loop)[0]));
            }
            ++idx;
         }
         file.close();
      }
      catch (std::exception&) {
         printf("Fatal error while attempting to load the ionization data file: %s.\n", filepath);
      }

      //int bytes = sizeof(float) * mIonizationEnergy.size();
      //checkCudaErrors(cudaMalloc(&d_IonizationEnergy, bytes));
      //checkCudaErrors(cudaMemcpyToSymbol(d_IonizationEnergy, mIonizationEnergy.data(), bytes, cudaMemcpyHostToDevice));
      //memcpy(mIonizationEnergy, hIonizationEnergy, sizeof(hIonizationEnergy));
   }

   void copyDataToCuda()
   {
      checkCudaErrors(cudaMemcpyToSymbol(d_AtomicWeight, mAtomicWeight.data(), sizeof(float) * mAtomicWeight.size()));
      checkCudaErrors(cudaMemcpyToSymbol(d_IonizationEnergy, mIonizationEnergy.data(), sizeof(float) * mIonizationEnergy.size()));
   }

   //struct Initializer
   //{
   //   Initializer()
   //   {
   //      readAtomicWeights();
   //      readIonizationEnergy();

   //      printf("Element::Initializer() completed: %d bytes\n", sizeof(mAllElements));
   //   }
   //};

   //Initializer init;

   void init()
   {
      readAtomicWeights();
      readIonizationEnergy();

      printf("InitializeElements() completed: %d bytes\n", sizeof(mAllElements));
   }

   __host__ __device__ Element::Element(int atomicNo)
   {
      if (atomicNo >= elmNone && atomicNo < elmEndOfElements) {
         mAtomicNumber = atomicNo;
      }
      else {
         printf("Wrong atomic number %d\n", atomicNo);
         mAtomicNumber = elmNone;
      }
   }

   Element::Element(const Element& other)
   {
      if (*this == other) return;
      mAtomicNumber = other.mAtomicNumber;
   }

   Element::Element()
   {
      mAtomicNumber = elmNone;
   }

   __host__ __device__ bool Element::operator==(const Element& other) const
   {
      return mAtomicNumber == other.mAtomicNumber;
   }

   const Element& Element::operator=(const Element& other)
   {
      if (*this == other) return *this;
      return byAtomicNumber(other.mAtomicNumber);
   }

   bool Element::operator<(const Element& other) const
   {
      return mAtomicNumber < other.mAtomicNumber;
   }

   static int atomicNumberForName(char const * name)
   {
      std::string strName(name);
      std::transform(strName.begin(), strName.end(), strName.begin(), ::tolower);
      for (int i = 0; i < elmEndOfElements + 1; ++i) {
         std::string nameLong = mElementNames[i];
         std::string nameAbbrev = mAbbreviations[i];
         std::transform(nameLong.begin(), nameLong.end(), nameLong.begin(), ::tolower);
         std::transform(nameAbbrev.begin(), nameAbbrev.end(), nameAbbrev.begin(), ::tolower);
         if (strName == nameLong || strName == nameAbbrev) {
            return i;
         }
      }
      return atoi(name);
   }

   const Element& byName(char const * name)
   {
      int z = atomicNumberForName(name);
      return z == 0 ? None : *mAllElements[z - 1];
   }

   const Element& byAtomicNumber(int an)
   {
      return (an > 0) && (an <= numAtomicWeight) ? *mAllElements[an - 1] : None;
   }

   __host__ __device__ float getAtomicWeight(int atomicNo)
   {
      if (atomicNo <= 0 || atomicNo >= numAtomicWeight) {
         printf("invalid atmoic number: %d\n", atomicNo);
         return elmNone;
      }
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      else if (!d_AtomicWeight || !d_AtomicWeight[atomicNo - 1]) {
         printf("need to load d_AtomicWeight on to device first\n");
         return elmNone;
      }
      return d_AtomicWeight[atomicNo - 1];
#else
      else if (!mAtomicWeight.size() || !mAtomicWeight[atomicNo - 1]) {
         readAtomicWeights();
         printf("need to load mAtomicWeight array by calling readAtomicWeights first\n");
      }
      return mAtomicWeight[atomicNo - 1];
#endif
   }

   //__device__ float getAtomicWeightDevice(int atomicNo)
   //{
   //   if (atomicNo <= 0 || atomicNo >= numAtomicWeight) {
   //      printf("invalid atmoic number: %d\n", atomicNo);
   //      return -1;
   //   }
   //   return d_AtomicWeight[atomicNo - 1];
   //}

   Element const * const * allElements()
   {
      return mAllElements;
   }

   //std::vector<const Element*> range(const Element& min, const Element& max)
   //{
   //   if (min.getAtomicNumber() <= max.getAtomicNumber()) {
   //      printf("make sure min < max when calling range");
   //      return std::vector<const Element*>();
   //   }
   //   int numElem = max.getAtomicNumber() - min.getAtomicNumber() + 1;
   //   std::vector<Element const *> res(numElem);
   //   for (int k = 0; k < numElem; ++k) {
   //      res[k] = mAllElements[min.getAtomicNumber() - 1 + k];
   //   }
   //   return res;
   //}

   //float Element::meanIonizationPotential(int atomicNo) {
   //   try {
   //      return MeanIonizationPotential.Berger64.compute(Element.byAtomicNumber(atomicNo));
   //   }
   //   catch (std::exception& ex) {
   //      return MeanIonizationPotential.Sternheimer64.compute(Element.byAtomicNumber(atomicNo));
   //   }
   //}

   __host__ __device__ int Element::getAtomicNumber() const
   {
      return mAtomicNumber;
   }

   __host__ __device__ float Element::getAtomicWeight() const
   {
      return ::Element::getAtomicWeight(mAtomicNumber);
   }

   __host__ __device__ float Element::getMass() const
   {
      return ToSI::AMU(::Element::getAtomicWeight(mAtomicNumber));
   }

   __host__ __device__ char const * Element::toAbbrev() const
   {
      return mAbbreviations[mAtomicNumber];
   }

   __host__ __device__ char const * toAbbrev(int atomicNo)
   {
      return mAbbreviations[atomicNo];
   }

   char const * toString(int el)
   {
      return mElementNames[el];
   }

   //float meanIonizationPotential()
   //{
   //   try {
   //      return MeanIonizationPotential.Berger64.compute(this);
   //   }
   //   catch (final Exception ex) {
   //      return MeanIonizationPotential.Sternheimer64.compute(this);
   //   }
   //}

   //public float energyLoss(float eK) {
   //   return BetheElectronEnergyLoss.JoyLuo1989.compute(this, eK);
   //}
   //

   //public float massAbsorptionCoefficient(float energy) {
   //   return AlgorithmUser.getDefaultMAC().compute(this, energy);
   //}
   //

   //public float massAbsorptionCoefficient(XRayTransition xrt)
   //throws EPQException{
   //   return AlgorithmUser.getDefaultMAC().compute(this, xrt);
   //}
   //

   bool isValid(int atomicNo)
   {
      return (atomicNo >= elmH) && (atomicNo < elmEndOfElements);
   }

   bool Element::isValid() const
   {
      return (mAtomicNumber >= elmH) && (mAtomicNumber < elmEndOfElements);
   }

   int Element::compareTo(const Element& e) const
   {
      if (mAtomicNumber < e.mAtomicNumber) {
         return -1;
      }
      else {
         return mAtomicNumber == e.mAtomicNumber ? 0 : 1;
      }
   }

   __host__ __device__ unsigned int Element::hashCode() const
   {
      // mAtomicNumber is always less than 128 (1<<7). Int has 31 + 1 bits. 31-7
      // = 24
      return mAtomicNumber << 24;
   }

   bool Element::equals(const Element& el)
   {
      return *this == el;
   }

   char const * Element::toString() const
   {
      return mElementNames[mAtomicNumber];
   }

   __host__ __device__ float Element::getIonizationEnergy() const
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      if (d_IonizationEnergy[mAtomicNumber - 1] <= 0) {
         printf("Element::getIonizationEnergy: load mIonizationEnergy by calling readIonizationEnergy first\n");
         return -1.0;
      }

      const float res = (mAtomicNumber - 1 <= numIonizationEnergy) ? d_IonizationEnergy[mAtomicNumber - 1] : -1.0;
      if (res == -1.0) {
         printf("Element::getIonizationEnergy: The ionization energy is not available for %s\n", toAbbrev());
      }
      return res;
#else
      if (mIonizationEnergy.size() <= 0 || mIonizationEnergy[mAtomicNumber - 1] <= 0) {
         readIonizationEnergy();
         printf("Element::getIonizationEnergy: load mIonizationEnergy by calling readIonizationEnergy first\n");
      }

      const float res = (mAtomicNumber - 1 <= numIonizationEnergy) ? mIonizationEnergy[mAtomicNumber - 1] : -1.0;
      if (res == -1.0) {
         printf("Element::getIonizationEnergy: The ionization energy is not available for %s\n", toAbbrev());
      }
      return res;
#endif
   }

   const Element& Element::readResolve()
   {
      return ::Element::byAtomicNumber(mAtomicNumber);
   }

   //ElementNameVecT getListOfAbbreviations(const Element& minEl, const Element& maxEl)
   //{
   //   int numEl = maxEl.getAtomicNumber() - minEl.getAtomicNumber() + 1;
   //   std::vector<std::string> res(numEl, "");
   //   for (int z = minEl.getAtomicNumber(); z <= maxEl.getAtomicNumber(); ++z) {
   //      std::string abbrev(toAbbrev(z));
   //      res[z - minEl.getAtomicNumber()] = abbrev;
   //   }
   //   return res;
   //}

   //ElementNameVecT getListOfElements(const Element& minEl, const Element& maxEl)
   //{
   //   std::vector<std::string> res;
   //   for (int z = minEl.getAtomicNumber(); z <= maxEl.getAtomicNumber(); ++z) {
   //      res.push_back(toString(z));
   //   }
   //   return res;
   //}

   __global__ void initCuda()
   {
      dNone = new Element(0);
      dH = new Element(1);
      dHe = new Element(2);
      dLi = new Element(3);
      dBe = new Element(4);
      dB = new Element(5);
      dC = new Element(6);
      dN = new Element(7);
      dO = new Element(8);
      dF = new Element(9);
      dNe = new Element(10);
      dNa = new Element(11);
      dMg = new Element(12);
      dAl = new Element(13);
      dSi = new Element(14);
      dP = new Element(15);
      dS = new Element(16);
      dCl = new Element(17);
      dAr = new Element(18);
      dK = new Element(19);
      dCa = new Element(20);
      dSc = new Element(21);
      dTi = new Element(22);
      dV = new Element(23);
      dCr = new Element(24);
      dMn = new Element(25);
      dFe = new Element(26);
      dCo = new Element(27);
      dNi = new Element(28);
      dCu = new Element(29);
      dZn = new Element(30);
      dGa = new Element(31);
      dGe = new Element(32);
      dAs = new Element(33);
      dSe = new Element(34);
      dBr = new Element(35);
      dKr = new Element(36);
      dRb = new Element(37);
      dSr = new Element(38);
      dY = new Element(39);
      dZr = new Element(40);
      dNb = new Element(41);
      dMo = new Element(42);
      dTc = new Element(43);
      dRu = new Element(44);
      dRh = new Element(45);
      dPd = new Element(46);
      dAg = new Element(47);
      dCd = new Element(48);
      dIn = new Element(49);
      dSn = new Element(50);
      dSb = new Element(51);
      dTe = new Element(52);
      dI = new Element(53);
      dXe = new Element(54);
      dCs = new Element(55);
      dBa = new Element(56);
      dLa = new Element(57);
      dCe = new Element(58);
      dPr = new Element(59);
      dNd = new Element(60);
      dPm = new Element(61);
      dSm = new Element(62);
      dEu = new Element(63);
      dGd = new Element(64);
      dTb = new Element(65);
      dDy = new Element(66);
      dHo = new Element(67);
      dEr = new Element(68);
      dTm = new Element(69);
      dYb = new Element(70);
      dLu = new Element(71);
      dHf = new Element(72);
      dTa = new Element(73);
      dW = new Element(74);
      dRe = new Element(75);
      dOs = new Element(76);
      dIr = new Element(77);
      dPt = new Element(78);
      dAu = new Element(79);
      dHg = new Element(80);
      dTl = new Element(81);
      dPb = new Element(82);
      dBi = new Element(83);
      dPo = new Element(84);
      dAt = new Element(85);
      dRn = new Element(86);
      dFr = new Element(87);
      dRa = new Element(88);
      dAc = new Element(89);
      dTh = new Element(90);
      dPa = new Element(91);
      dU = new Element(92);
      dNp = new Element(93);
      dPu = new Element(94);
      dAm = new Element(95);
      dCm = new Element(96);
      dBk = new Element(97);
      dCf = new Element(98);
      dEs = new Element(99);
      dFm = new Element(100);
      dMd = new Element(101);
      dNo = new Element(102);
      dLr = new Element(103);
      dRf = new Element(104);
      dDb = new Element(105);
      dSg = new Element(106);
      dBh = new Element(107);
      dHs = new Element(108);
      dMt = new Element(109);
      dUun = new Element(110);
      dUuu = new Element(111);
      dUub = new Element(112);
   }
}
