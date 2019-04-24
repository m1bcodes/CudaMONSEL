#include "Element.cuh"

#include <algorithm>

namespace Element
{
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

   const Element None = Element(0);
   const Element H = Element(1);
   const Element He = Element(2);
   const Element Li = Element(3);
   const Element Be = Element(4);
   const Element B = Element(5);
   const Element C = Element(6);
   const Element N = Element(7);
   const Element O = Element(8);
   const Element F = Element(9);
   const Element Ne = Element(10);
   const Element Na = Element(11);
   const Element Mg = Element(12);
   const Element Al = Element(13);
   const Element Si = Element(14);
   const Element P = Element(15);
   const Element S = Element(16);
   const Element Cl = Element(17);
   const Element Ar = Element(18);
   const Element K = Element(19);
   const Element Ca = Element(20);
   const Element Sc = Element(21);
   const Element Ti = Element(22);
   const Element V = Element(23);
   const Element Cr = Element(24);
   const Element Mn = Element(25);
   const Element Fe = Element(26);
   const Element Co = Element(27);
   const Element Ni = Element(28);
   const Element Cu = Element(29);
   const Element Zn = Element(30);
   const Element Ga = Element(31);
   const Element Ge = Element(32);
   const Element As = Element(33);
   const Element Se = Element(34);
   const Element Br = Element(35);
   const Element Kr = Element(36);
   const Element Rb = Element(37);
   const Element Sr = Element(38);
   const Element Y = Element(39);
   const Element Zr = Element(40);
   const Element Nb = Element(41);
   const Element Mo = Element(42);
   const Element Tc = Element(43);
   const Element Ru = Element(44);
   const Element Rh = Element(45);
   const Element Pd = Element(46);
   const Element Ag = Element(47);
   const Element Cd = Element(48);
   const Element In = Element(49);
   const Element Sn = Element(50);
   const Element Sb = Element(51);
   const Element Te = Element(52);
   const Element I = Element(53);
   const Element Xe = Element(54);
   const Element Cs = Element(55);
   const Element Ba = Element(56);
   const Element La = Element(57);
   const Element Ce = Element(58);
   const Element Pr = Element(59);
   const Element Nd = Element(60);
   const Element Pm = Element(61);
   const Element Sm = Element(62);
   const Element Eu = Element(63);
   const Element Gd = Element(64);
   const Element Tb = Element(65);
   const Element Dy = Element(66);
   const Element Ho = Element(67);
   const Element Er = Element(68);
   const Element Tm = Element(69);
   const Element Yb = Element(70);
   const Element Lu = Element(71);
   const Element Hf = Element(72);
   const Element Ta = Element(73);
   const Element W = Element(74);
   const Element Re = Element(75);
   const Element Os = Element(76);
   const Element Ir = Element(77);
   const Element Pt = Element(78);
   const Element Au = Element(79);
   const Element Hg = Element(80);
   const Element Tl = Element(81);
   const Element Pb = Element(82);
   const Element Bi = Element(83);
   const Element Po = Element(84);
   const Element At = Element(85);
   const Element Rn = Element(86);
   const Element Fr = Element(87);
   const Element Ra = Element(88);
   const Element Ac = Element(89);
   const Element Th = Element(90);
   const Element Pa = Element(91);
   const Element U = Element(92);
   const Element Np = Element(93);
   const Element Pu = Element(94);
   const Element Am = Element(95);
   const Element Cm = Element(96);
   const Element Bk = Element(97);
   const Element Cf = Element(98);
   const Element Es = Element(99);
   const Element Fm = Element(100);
   const Element Md = Element(101);
   const Element No = Element(102);
   const Element Lr = Element(103);
   const Element Rf = Element(104);
   const Element Db = Element(105);
   const Element Sg = Element(106);
   const Element Bh = Element(107);
   const Element Hs = Element(108);
   const Element Mt = Element(109);
   const Element Uun = Element(110);
   const Element Uuu = Element(111);
   const Element Uub = Element(112);

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

   Element mAllElements[numAtomicWeight];
   float mIonizationEnergy[numIonizationEnergy];
   float mAtomicWeight[numAtomicWeight];

   void InitializeElements()
   {
      readAtomicWeights();
      readIonizationEnergy();

      Element tmpAllElements[] = {
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
      memcpy(mAllElements, tmpAllElements, sizeof(tmpAllElements));
      printf("InitializeElements() completed: %d bytes\n", sizeof(mAllElements));
   }

   Element::Element(int atomicNo)
   {
      if (atomicNo >= elmNone && atomicNo < elmEndOfElements) {
         mAtomicNumber = atomicNo;
      }
      else {
         printf("Wrong atomic number %d\n", atomicNo);
      }
   }

   Element::Element()
   {
      mAtomicNumber = elmNone;
   }

   bool Element::operator==(const Element& other) const
   {
      return mAtomicNumber == other.mAtomicNumber;
   }

   Element::Element(const Element& other)
   {
      if (*this == other) return;
      mAtomicNumber = other.mAtomicNumber;
   }

   Element& Element::operator=(const Element& other)
   {
      if (*this == other) return *this;
      mAtomicNumber = other.mAtomicNumber;
   }

   bool Element::operator<(const Element& other) const
   {
      return mAtomicNumber < other.mAtomicNumber;
   }

   void readAtomicWeights()
   {
      float hAtomicWeight[numAtomicWeight];
      try {
         char filepath[] = ".\\gov\\nist\\microanalysis\\EPQLibrary\\AtomicWeights.csv";
         std::ifstream file(filepath);
         printf("Reading: %s\n", filepath);
         int idx = 0;
         for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
            if ((*loop)[0][0] == '/' && (*loop)[0][1] == '/') { // check if the first line should be removed
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
      memcpy(mAtomicWeight, hAtomicWeight, sizeof(hAtomicWeight));
      //checkCudaErrors(cudaMemcpyToSymbol(mAtomicWeight, &hAtomicWeight, sizeof(float) * 112));
   }

   int atomicNumberForName(char const * name)
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

   Element byName(char const * name)
   {
      int z = atomicNumberForName(name);
      return z == 0 ? None : mAllElements[z - 1];
   }

   Element byAtomicNumber(int an)
   {      
      return (an > 0) && (an < numAtomicWeight) ? mAllElements[an - 1] : None;
   }

   double getAtomicWeight(int atomicNo)
   {
      //printf("atomicNo: %d\n", atomicNo);
      if (atomicNo <= 0 || atomicNo >= numAtomicWeight) {
         //printf("invalid atmoic number: %d\n", atomicNo);
         return -1;
      }
      else if (!mAtomicWeight || !mAtomicWeight[atomicNo - 1]) {
         readAtomicWeights();
         printf("need to load mAtomicWeight array by calling readAtomicWeights first\n");
      }
      return mAtomicWeight[atomicNo - 1];
   }

   Element const * allElements()
   {
      return mAllElements;
   }

   std::vector<Element> range(const Element& min, const Element& max)
   {
      if (min.getAtomicNumber() <= max.getAtomicNumber()) {
         printf("make sure min < max when calling range");
         return std::vector<Element>();
      }
      int numElem = max.getAtomicNumber() - min.getAtomicNumber() + 1;
      std::vector<Element> res(numElem);
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

   int Element::getAtomicNumber() const
   {
      return mAtomicNumber;
   }

   double Element::getAtomicWeight() const
   {
      return ::Element::getAtomicWeight(mAtomicNumber);
   }

   double Element::getMass() const
   {
      return ToSI::AMU(::Element::getAtomicWeight(mAtomicNumber));
   }

   char const * Element::toAbbrev() const
   {
      return mAbbreviations[mAtomicNumber];
   }

   char const * toAbbrev(int atomicNo)
   {
      return mAbbreviations[atomicNo];
   }

   char const * toString(int el)
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

   unsigned int Element::hashCode() const
   {
      // mAtomicNumber is always less than 128 (1<<7). Int has 31 + 1 bits. 31-7
      // = 24
      return mAtomicNumber << 24;
   }

   bool Element::equals(const Element& el)
   {
      return *this == el;
   }

   char const * Element::toString()
   {
      return mElementNames[mAtomicNumber];
   }

   void readIonizationEnergy()
   {
      char filepath[] = ".\\gov\\nist\\microanalysis\\EPQLibrary\\IonizationEnergies.csv";
      float hIonizationEnergy[numIonizationEnergy];
      try {
         printf("Reading: %s\n", filepath);
         std::ifstream file(filepath);
         int idx = 0;
         for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
            if ((*loop)[0][0] == '/' && (*loop)[0][1] == '/') { // check if the first line should be removed
               continue;
            }
            if (CSVIterator::IsNaN((*loop)[0])) {
               hIonizationEnergy[idx] = -1.0;
            }
            else {
               hIonizationEnergy[idx] = std::stof((*loop)[0]);
            }
            ++idx;
         }
      }
      catch (std::exception&) {
         printf("Fatal error while attempting to load the ionization data file: %s.\n", filepath);
      }
      memcpy(mIonizationEnergy, hIonizationEnergy, sizeof(hIonizationEnergy));
   }

   double Element::getIonizationEnergy() const
   {
      if (mIonizationEnergy <= 0 || mIonizationEnergy[mAtomicNumber - 1] <= 0) {
         readIonizationEnergy();
         printf("load mIonizationEnergy by calling readIonizationEnergy first\n");
      }

      double res = (mAtomicNumber - 1 <= numIonizationEnergy) ? mIonizationEnergy[mAtomicNumber - 1] : -1.0;
      if (res == -1.0) {
         printf("The ionization energy is not available for %s\n", toAbbrev()); // new EPQFatalException(");
      }
      return res;
   }

   Element Element::readResolve()
   {
      return ::Element::byAtomicNumber(mAtomicNumber);
   }

   ElementNameVecT getListOfAbbreviations(const Element& minEl, const Element& maxEl)
   {
      int numEl = maxEl.getAtomicNumber() - minEl.getAtomicNumber() + 1;
      std::vector<std::string> res(numEl, "");
      for (int z = minEl.getAtomicNumber(); z <= maxEl.getAtomicNumber(); ++z) {
         std::string abbrev(toAbbrev(z));
         res[z - minEl.getAtomicNumber()] = abbrev;
      }
      return res;
   }

   ElementNameVecT getListOfElements(const Element& minEl, const Element& maxEl)
   {
      std::vector<std::string> res;
      for (int z = minEl.getAtomicNumber(); z <= maxEl.getAtomicNumber(); ++z) {
         res.push_back(toString(z));
      }
      return res;
   }
}
