#include "Element.cuh"
#include "..\..\..\..\CudaUtil.h"
#include "Amphibian\String.cuh"

namespace Element
{
   static const int numIonizationEnergy = 104;
   static const int numAtomicWeight = 112;

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

   __device__ const Element None;
   __device__ const Element H;
   __device__ const Element He;
   __device__ const Element Li;
   __device__ const Element Be;
   __device__ const Element B;
   __device__ const Element C;
   __device__ const Element N;
   __device__ const Element O;
   __device__ const Element F;
   __device__ const Element Ne;
   __device__ const Element Na;
   __device__ const Element Mg;
   __device__ const Element Al;
   __device__ const Element Si;
   __device__ const Element P;
   __device__ const Element S;
   __device__ const Element Cl;
   __device__ const Element Ar;
   __device__ const Element K;
   __device__ const Element Ca;
   __device__ const Element Sc;
   __device__ const Element Ti;
   __device__ const Element V;
   __device__ const Element Cr;
   __device__ const Element Mn;
   __device__ const Element Fe;
   __device__ const Element Co;
   __device__ const Element Ni;
   __device__ const Element Cu;
   __device__ const Element Zn;
   __device__ const Element Ga;
   __device__ const Element Ge;
   __device__ const Element As;
   __device__ const Element Se;
   __device__ const Element Br;
   __device__ const Element Kr;
   __device__ const Element Rb;
   __device__ const Element Sr;
   __device__ const Element Y;
   __device__ const Element Zr;
   __device__ const Element Nb;
   __device__ const Element Mo;
   __device__ const Element Tc;
   __device__ const Element Ru;
   __device__ const Element Rh;
   __device__ const Element Pd;
   __device__ const Element Ag;
   __device__ const Element Cd;
   __device__ const Element In;
   __device__ const Element Sn;
   __device__ const Element Sb;
   __device__ const Element Te;
   __device__ const Element I;
   __device__ const Element Xe;
   __device__ const Element Cs;
   __device__ const Element Ba;
   __device__ const Element La;
   __device__ const Element Ce;
   __device__ const Element Pr;
   __device__ const Element Nd;
   __device__ const Element Pm;
   __device__ const Element Sm;
   __device__ const Element Eu;
   __device__ const Element Gd;
   __device__ const Element Tb;
   __device__ const Element Dy;
   __device__ const Element Ho;
   __device__ const Element Er;
   __device__ const Element Tm;
   __device__ const Element Yb;
   __device__ const Element Lu;
   __device__ const Element Hf;
   __device__ const Element Ta;
   __device__ const Element W;
   __device__ const Element Re;
   __device__ const Element Os;
   __device__ const Element Ir;
   __device__ const Element Pt;
   __device__ const Element Au;
   __device__ const Element Hg;
   __device__ const Element Tl;
   __device__ const Element Pb;
   __device__ const Element Bi;
   __device__ const Element Po;
   __device__ const Element At;
   __device__ const Element Rn;
   __device__ const Element Fr;
   __device__ const Element Ra;
   __device__ const Element Ac;
   __device__ const Element Th;
   __device__ const Element Pa;
   __device__ const Element U;
   __device__ const Element Np;
   __device__ const Element Pu;
   __device__ const Element Am;
   __device__ const Element Cm;
   __device__ const Element Bk;
   __device__ const Element Cf;
   __device__ const Element Es;
   __device__ const Element Fm;
   __device__ const Element Md;
   __device__ const Element No;
   __device__ const Element Lr;
   __device__ const Element Rf;
   __device__ const Element Db;
   __device__ const Element Sg;
   __device__ const Element Bh;
   __device__ const Element Hs;
   __device__ const Element Mt;
   __device__ const Element Uun;
   __device__ const Element Uuu;
   __device__ const Element Uub;

   __device__ const Element mAllElements[112];
   /* = {
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
   };*/

   __device__ char const * const mElementNames[] = {
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

   __device__ char const * const mAbbreviations[] = {
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

   __device__ float mIonizationEnergy[104];
   __device__ float mAtomicWeight[112];

   __host__ __device__ Element::Element(int atomicNo)
   {
      if ((atomicNo >= elmNone) && (atomicNo < elmEndOfElements)) {
         mAtomicNumber = atomicNo;
      }
      else {
         printf("Wrong atomic number %d\n", atomicNo);
      }
   }

   __host__ __device__ Element::Element()
   {
      // mAtomicNumber = elmNone;
   }

   __host__ __device__ bool Element::operator==(const Element& other)
   {
      return mAtomicNumber == other.mAtomicNumber;
   }

   __host__ void readAtomicWeights()
   {
      float hAtomicWeight[112];
      try {
         std::ifstream file(".\\gov\\nist\\microanalysis\\EPQLibrary\\AtomicWeights.csv");
         int idx = 0;
         for (CSVIterator loop(file); loop != CSVIterator(); ++loop) { // TODO: check if the first line should be removed
            if ((*loop)[0][0] == '/' && (*loop)[0][1] == '/') {
               continue;
            }
            hAtomicWeight[idx] = std::stof((*loop)[0]);
            ++idx;
         }
      }
      catch (std::exception&) {
         printf("didnt see file AtomicWeights.csv\n");
         throw 0; //throw new EPQFatalException("Fatal error while attempting to load the atomic weights data file.");
      }
      checkCudaErrors(cudaMemcpyToSymbol(mAtomicWeight, &hAtomicWeight, sizeof(float) * 112));
   }

   __device__ int atomicNumberForName(char* name)
   {
      for (int i = 0; i < elmEndOfElements + 1; ++i) {
         if (String::AreEqual(mElementNames[i], name) || String::AreEqual(mAbbreviations[i], name)) { // TODO: make it case insensitive
            return i;
         }
      }
      return String::AToI(name);
   }

   __device__ Element byName(char* name)
   {
      int z = atomicNumberForName(name);
      return z == 0 ? None : mAllElements[z - 1];
   }

   __device__ Element byAtomicNumber(int an)
   {
      return (an >= 1) && (an < elmEndOfElements - 1) ? mAllElements[an - 1] : None;
   }

   __device__ double getAtomicWeight(int atomicNo)
   {
      if (mAtomicWeight == NULL) {
         printf("need to load mAtomicWeight array by calling readAtomicWeights first"); //readAtomicWeights();
      }
      return mAtomicWeight[atomicNo - 1];
   }

   __device__ Element const * allElements()
   {
      return mAllElements;
   }

   __device__ Element* range(Element min, Element max)
   {
      if (min.getAtomicNumber() <= max.getAtomicNumber()) {
         printf("make sure min < max when calling range");
         return NULL;
      }
      int numElem = max.getAtomicNumber() - min.getAtomicNumber() + 1;
      Element* res = new Element[numElem];
      for (int k = 0; k < numElem; ++k) {
         res[k] = mAllElements[min.getAtomicNumber() - 1 + k];
      }
      return res;
   }

   //double Element::meanIonizationPotential(int atomicNo) {
   //   try {
   //      return MeanIonizationPotential.Berger64.compute(Element.byAtomicNumber(atomicNo));
   //   }
   //   catch (std::exception& ex) {
   //      return MeanIonizationPotential.Sternheimer64.compute(Element.byAtomicNumber(atomicNo));
   //   }
   //}

   __device__ int Element::getAtomicNumber() {
      return mAtomicNumber;
   }

   __device__ double Element::getAtomicWeight()
   {
      return ::Element::getAtomicWeight(mAtomicNumber);
   }

   __device__ double Element::getMass()
   {
      return ToSI::AMU(::Element::getAtomicWeight(mAtomicNumber));
   }

   __device__ char const * Element::toAbbrev()
   {
      return mAbbreviations[mAtomicNumber];
   }

   __device__ char const * toAbbrev(int atomicNo)
   {
      return mAbbreviations[atomicNo];
   }

   __device__ char const * toString(int el)
   {
      return mElementNames[el];
   }

   //double meanIonizationPotential()
   //{
   //   try {
   //      return MeanIonizationPotential.Berger64.compute(this);
   //   }
   //   catch (final Exception ex) {
   //      return MeanIonizationPotential.Sternheimer64.compute(this);
   //   }
   //}

   //public double energyLoss(double eK) {
   //   return BetheElectronEnergyLoss.JoyLuo1989.compute(this, eK);
   //}
   //

   //public double massAbsorptionCoefficient(double energy) {
   //   return AlgorithmUser.getDefaultMAC().compute(this, energy);
   //}
   //

   //public double massAbsorptionCoefficient(XRayTransition xrt)
   //throws EPQException{
   //   return AlgorithmUser.getDefaultMAC().compute(this, xrt);
   //}
   //

   __device__ bool isValid(int atomicNo)
   {
      return (atomicNo >= elmH) && (atomicNo < elmEndOfElements);
   }

   __device__ bool Element::isValid()
   {
      return (mAtomicNumber >= elmH) && (mAtomicNumber < elmEndOfElements);
   }

   __device__ int Element::compareTo(Element e)
   {
      if (mAtomicNumber < e.mAtomicNumber) {
         return -1;
      }
      else {
         return mAtomicNumber == e.mAtomicNumber ? 0 : 1;
      }
   }

   __device__ int Element::hashCode()
   {
      // mAtomicNumber is always less than 128 (1<<7). Int has 31 + 1 bits. 31-7
      // = 24
      return mAtomicNumber << 24;
   }

   __device__ bool Element::equals(const Element& el)
   {
      return *this == el;
   }

   __device__ char const * Element::toString()
   {
      return mElementNames[mAtomicNumber];
   }

   __host__ void readIonizationEnergy()
   {
      float hIonizationEnergy[104];
      try {
         std::ifstream file(".\\gov\\nist\\microanalysis\\EPQLibrary\\IonizationEnergies.csv");
         int idx = 0;
         for (CSVIterator loop(file); loop != CSVIterator(); ++loop) { // TODO: check if the first line should be removed
            if ((*loop)[0][0] == '/' && (*loop)[0][1] == '/') {
               continue;
            }
            if (CSVIterator::IsNaN((*loop)[0])) {
               hIonizationEnergy[idx] = -1.0;
            }
            else {
               hIonizationEnergy[idx] = std::stof((*loop)[0]);
            }
         }
         ++idx;
      }
      catch (std::exception&) {
         printf("Fatal error while attempting to load the atomic weights data file.");
      }
      checkCudaErrors(cudaMemcpyToSymbol(mIonizationEnergy, &hIonizationEnergy, sizeof(float) * 104));
   }

   __device__ double Element::getIonizationEnergy()
   {
      if (mIonizationEnergy == NULL) {
         printf("load mIonizationEnergy by calling readIonizationEnergy first");
      }

      double res = (mAtomicNumber - 1 <= 104) ? mIonizationEnergy[mAtomicNumber - 1] : -1.0;
      if (res == -1.0) {
         printf("The ionization energy is not available for %s", toAbbrev()); // new EPQFatalException(");
      }
      return res;
   }

   __device__ Element Element::readResolve()
   {
      return ::Element::byAtomicNumber(mAtomicNumber);
   }

   //char * const Element::getListOfAbbreviations(Element minEl, Element maxEl)
   //{
   //   int numEl = maxEl.getAtomicNumber() - minEl.getAtomicNumber() + 1;
   //   char *res[] = new char*[numEl];
   //   for (int z = minEl.getAtomicNumber(); z <= maxEl.getAtomicNumber(); ++z)
   //      res.add(toAbbrev(z));
   //   return res;
   //}

   //static public final ArrayList<String> getListOfElements(Element minEl, Element maxEl)
   //{
   //   final ArrayList<String> res = new ArrayList<String>();
   //   for (int z = minEl.getAtomicNumber(); z <= maxEl.getAtomicNumber(); ++z)
   //      res.add(toString(z));
   //   return res;
   //}

   //__device__ bool AreEqual(Element& e1, Element& e2)
   //{
   //   if (&e1 == &e2) return true;
   //   return e1.getAtomicNumber() == e2.getAtomicNumber();
   //}
   
   __host__ void InitializeElements()
   {
      const Element hNone(0);
      const Element hH(1);
      const Element hHe(2);
      const Element hLi(3);
      const Element hBe(4);
      const Element hB(5);
      const Element hC(6);
      const Element hN(7);
      const Element hO(8);
      const Element hF(9);
      const Element hNe(10);
      const Element hNa(11);
      const Element hMg(12);
      const Element hAl(13);
      const Element hSi(14);
      const Element hP(15);
      const Element hS(16);
      const Element hCl(17);
      const Element hAr(18);
      const Element hK(19);
      const Element hCa(20);
      const Element hSc(21);
      const Element hTi(22);
      const Element hV(23);
      const Element hCr(24);
      const Element hMn(25);
      const Element hFe(26);
      const Element hCo(27);
      const Element hNi(28);
      const Element hCu(29);
      const Element hZn(30);
      const Element hGa(31);
      const Element hGe(32);
      const Element hAs(33);
      const Element hSe(34);
      const Element hBr(35);
      const Element hKr(36);
      const Element hRb(37);
      const Element hSr(38);
      const Element hY(39);
      const Element hZr(40);
      const Element hNb(41);
      const Element hMo(42);
      const Element hTc(43);
      const Element hRu(44);
      const Element hRh(45);
      const Element hPd(46);
      const Element hAg(47);
      const Element hCd(48);
      const Element hIn(49);
      const Element hSn(50);
      const Element hSb(51);
      const Element hTe(52);
      const Element hI(53);
      const Element hXe(54);
      const Element hCs(55);
      const Element hBa(56);
      const Element hLa(57);
      const Element hCe(58);
      const Element hPr(59);
      const Element hNd(60);
      const Element hPm(61);
      const Element hSm(62);
      const Element hEu(63);
      const Element hGd(64);
      const Element hTb(65);
      const Element hDy(66);
      const Element hHo(67);
      const Element hEr(68);
      const Element hTm(69);
      const Element hYb(70);
      const Element hLu(71);
      const Element hHf(72);
      const Element hTa(73);
      const Element hW(74);
      const Element hRe(75);
      const Element hOs(76);
      const Element hIr(77);
      const Element hPt(78);
      const Element hAu(79);
      const Element hHg(80);
      const Element hTl(81);
      const Element hPb(82);
      const Element hBi(83);
      const Element hPo(84);
      const Element hAt(85);
      const Element hRn(86);
      const Element hFr(87);
      const Element hRa(88);
      const Element hAc(89);
      const Element hTh(90);
      const Element hPa(91);
      const Element hU(92);
      const Element hNp(93);
      const Element hPu(94);
      const Element hAm(95);
      const Element hCm(96);
      const Element hBk(97);
      const Element hCf(98);
      const Element hEs(99);
      const Element hFm(100);
      const Element hMd(101);
      const Element hNo(102);
      const Element hLr(103);
      const Element hRf(104);
      const Element hDb(105);
      const Element hSg(106);
      const Element hBh(107);
      const Element hHs(108);
      const Element hMt(109);
      const Element hUun(110);
      const Element hUuu(111);
      const Element hUub(112);

      checkCudaErrors(cudaMemcpyToSymbol(None, &hNone, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(H, &hH, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(He, &hHe, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Li, &hLi, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Be, &hBe, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(B, &hB, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(C, &hC, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(N, &hN, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(O, &hO, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(F, &hF, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ne, &hNe, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Na, &hNa, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Mg, &hMg, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Al, &hAl, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Si, &hSi, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(P, &hP, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(S, &hS, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Cl, &hCl, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ar, &hAr, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(K, &hK, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ca, &hCa, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Sc, &hSc, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ti, &hTi, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(V, &hV, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Cr, &hCr, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Mn, &hMn, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Fe, &hFe, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Co, &hCo, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ni, &hNi, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Cu, &hCu, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Zn, &hZn, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ga, &hGa, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ge, &hGe, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(As, &hAs, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Se, &hSe, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Br, &hBr, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Kr, &hKr, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Rb, &hRb, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Sr, &hSr, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Y, &hY, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Zr, &hZr, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Nb, &hNb, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Mo, &hMo, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Tc, &hTc, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ru, &hRu, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Rh, &hRh, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Pd, &hPd, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ag, &hAg, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Cd, &hCd, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(In, &hIn, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Sn, &hSn, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Sb, &hSb, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Te, &hTe, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(I, &hI, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Xe, &hXe, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Cs, &hCs, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ba, &hBa, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(La, &hLa, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ce, &hCe, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Pr, &hPr, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Nd, &hNd, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Pm, &hPm, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Sm, &hSm, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Eu, &hEu, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Gd, &hGd, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Tb, &hTb, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Dy, &hDy, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ho, &hHo, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Er, &hEr, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Tm, &hTm, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Yb, &hYb, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Lu, &hLu, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Hf, &hHf, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ta, &hTa, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(W, &hW, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Re, &hRe, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Os, &hOs, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ir, &hIr, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Pt, &hPt, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Au, &hAu, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Hg, &hHg, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Tl, &hTl, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Pb, &hPb, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Bi, &hBi, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Po, &hPo, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(At, &hAt, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Rn, &hRn, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Fr, &hFr, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ra, &hRa, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Ac, &hAc, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Th, &hTh, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Pa, &hPa, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(U, &hU, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Np, &hNp, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Pu, &hPu, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Am, &hAm, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Cm, &hCm, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Bk, &hBk, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Cf, &hCf, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Es, &hEs, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Fm, &hFm, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Md, &hMd, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(No, &hNo, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Lr, &hLr, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Rf, &hRf, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Db, &hDb, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Sg, &hSg, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Bh, &hBh, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Hs, &hHs, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Mt, &hMt, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Uun, &hUun, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Uuu, &hUuu, sizeof(Element)));
      checkCudaErrors(cudaMemcpyToSymbol(Uub, &hUub, sizeof(Element)));

      Element hAllElements[] = {
         hH,
         hHe,
         hLi,
         hBe,
         hB,
         hC,
         hN,
         hO,
         hF,
         hNe,
         hNa,
         hMg,
         hAl,
         hSi,
         hP,
         hS,
         hCl,
         hAr,
         hK,
         hCa,
         hSc,
         hTi,
         hV,
         hCr,
         hMn,
         hFe,
         hCo,
         hNi,
         hCu,
         hZn,
         hGa,
         hGe,
         hAs,
         hSe,
         hBr,
         hKr,
         hRb,
         hSr,
         hY,
         hZr,
         hNb,
         hMo,
         hTc,
         hRu,
         hRh,
         hPd,
         hAg,
         hCd,
         hIn,
         hSn,
         hSb,
         hTe,
         hI,
         hXe,
         hCs,
         hBa,
         hLa,
         hCe,
         hPr,
         hNd,
         hPm,
         hSm,
         hEu,
         hGd,
         hTb,
         hDy,
         hHo,
         hEr,
         hTm,
         hYb,
         hLu,
         hHf,
         hTa,
         hW,
         hRe,
         hOs,
         hIr,
         hPt,
         hAu,
         hHg,
         hTl,
         hPb,
         hBi,
         hPo,
         hAt,
         hRn,
         hFr,
         hRa,
         hAc,
         hTh,
         hPa,
         hU,
         hNp,
         hPu,
         hAm,
         hCm,
         hBk,
         hCf,
         hEs,
         hFm,
         hMd,
         hNo,
         hLr,
         hRf,
         hDb,
         hSg,
         hBh,
         hHs,
         hMt,
         hUun,
         hUuu,
         hUub
      };
      checkCudaErrors(cudaMemcpyToSymbol(mAllElements, &hAllElements, sizeof(Element) * 112));
   }
}
