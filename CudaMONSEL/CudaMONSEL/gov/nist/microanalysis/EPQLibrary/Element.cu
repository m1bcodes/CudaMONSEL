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

   static void readAtomicWeights()
   {
      if (!mAtomicWeight.empty()) return;
      mAtomicWeight.resize(numAtomicWeight);

      //float hAtomicWeight[numAtomicWeight];
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
      //memcpy(mAtomicWeight, hAtomicWeight, sizeof(hAtomicWeight));
      //checkCudaErrors(cudaMemcpyToSymbol(mAtomicWeight, &hAtomicWeight, sizeof(float) * 112));
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
      //memcpy(mIonizationEnergy, hIonizationEnergy, sizeof(hIonizationEnergy));
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

   double getAtomicWeight(int atomicNo)
   {
      //printf("atomicNo: %d\n", atomicNo);
      if (atomicNo <= 0 || atomicNo >= numAtomicWeight) {
         //printf("invalid atmoic number: %d\n", atomicNo);
         return -1;
      }
      else if (!mAtomicWeight.size() || !mAtomicWeight[atomicNo - 1]) {
         readAtomicWeights();
         printf("need to load mAtomicWeight array by calling readAtomicWeights first\n");
      }
      return mAtomicWeight[atomicNo - 1];
   }

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

   char const * Element::toString() const
   {
      return mElementNames[mAtomicNumber];
   }

   double Element::getIonizationEnergy() const
   {
      if (mIonizationEnergy.size() <= 0 || mIonizationEnergy[mAtomicNumber - 1] <= 0) {
         readIonizationEnergy();
         printf("Element::getIonizationEnergy: load mIonizationEnergy by calling readIonizationEnergy first\n");
      }

      double res = (mAtomicNumber - 1 <= numIonizationEnergy) ? mIonizationEnergy[mAtomicNumber - 1] : -1.0;
      if (res == -1.0) {
         printf("Element::getIonizationEnergy: The ionization energy is not available for %s\n", toAbbrev());
      }
      return res;
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
}
