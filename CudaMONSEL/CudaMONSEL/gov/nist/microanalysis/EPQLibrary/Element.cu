#include "Element.cuh"

const long long Element::serialVersionUID = 0x987360133793L;

float Element::mIonizationEnergy[104];
float Element::mAtomicWeight[112];

const Element Element::None(0);
const Element Element::H(1);
const Element Element::He(2);
const Element Element::Li(3);
const Element Element::Be(4);
const Element Element::B(5);
const Element Element::C(6);
const Element Element::N(7);
const Element Element::O(8);
const Element Element::F(9);
const Element Element::Ne(10);
const Element Element::Na(11);
const Element Element::Mg(12);
const Element Element::Al(13);
const Element Element::Si(14);
const Element Element::P(15);
const Element Element::S(16);
const Element Element::Cl(17);
const Element Element::Ar(18);
const Element Element::K(19);
const Element Element::Ca(20);
const Element Element::Sc(21);
const Element Element::Ti(22);
const Element Element::V(23);
const Element Element::Cr(24);
const Element Element::Mn(25);
const Element Element::Fe(26);
const Element Element::Co(27);
const Element Element::Ni(28);
const Element Element::Cu(29);
const Element Element::Zn(30);
const Element Element::Ga(31);
const Element Element::Ge(32);
const Element Element::As(33);
const Element Element::Se(34);
const Element Element::Br(35);
const Element Element::Kr(36);
const Element Element::Rb(37);
const Element Element::Sr(38);
const Element Element::Y(39);
const Element Element::Zr(40);
const Element Element::Nb(41);
const Element Element::Mo(42);
const Element Element::Tc(43);
const Element Element::Ru(44);
const Element Element::Rh(45);
const Element Element::Pd(46);
const Element Element::Ag(47);
const Element Element::Cd(48);
const Element Element::In(49);
const Element Element::Sn(50);
const Element Element::Sb(51);
const Element Element::Te(52);
const Element Element::I(53);
const Element Element::Xe(54);
const Element Element::Cs(55);
const Element Element::Ba(56);
const Element Element::La(57);
const Element Element::Ce(58);
const Element Element::Pr(59);
const Element Element::Nd(60);
const Element Element::Pm(61);
const Element Element::Sm(62);
const Element Element::Eu(63);
const Element Element::Gd(64);
const Element Element::Tb(65);
const Element Element::Dy(66);
const Element Element::Ho(67);
const Element Element::Er(68);
const Element Element::Tm(69);
const Element Element::Yb(70);
const Element Element::Lu(71);
const Element Element::Hf(72);
const Element Element::Ta(73);
const Element Element::W(74);
const Element Element::Re(75);
const Element Element::Os(76);
const Element Element::Ir(77);
const Element Element::Pt(78);
const Element Element::Au(79);
const Element Element::Hg(80);
const Element Element::Tl(81);
const Element Element::Pb(82);
const Element Element::Bi(83);
const Element Element::Po(84);
const Element Element::At(85);
const Element Element::Rn(86);
const Element Element::Fr(87);
const Element Element::Ra(88);
const Element Element::Ac(89);
const Element Element::Th(90);
const Element Element::Pa(91);
const Element Element::U(92);
const Element Element::Np(93);
const Element Element::Pu(94);
const Element Element::Am(95);
const Element Element::Cm(96);
const Element Element::Bk(97);
const Element Element::Cf(98);
const Element Element::Es(99);
const Element Element::Fm(100);
const Element Element::Md(101);
const Element Element::No(102);
const Element Element::Lr(103);
const Element Element::Rf(104);
const Element Element::Db(105);
const Element Element::Sg(106);
const Element Element::Bh(107);
const Element Element::Hs(108);
const Element Element::Mt(109);
const Element Element::Uun(110);
const Element Element::Uuu(111);
const Element Element::Uub(112);

const Element Element::mAllElements[] = {
   Element::H,
   Element::He,
   Element::Li,
   Element::Be,
   Element::B,
   Element::C,
   Element::N,
   Element::O,
   Element::F,
   Element::Ne,
   Element::Na,
   Element::Mg,
   Element::Al,
   Element::Si,
   Element::P,
   Element::S,
   Element::Cl,
   Element::Ar,
   Element::K,
   Element::Ca,
   Element::Sc,
   Element::Ti,
   Element::V,
   Element::Cr,
   Element::Mn,
   Element::Fe,
   Element::Co,
   Element::Ni,
   Element::Cu,
   Element::Zn,
   Element::Ga,
   Element::Ge,
   Element::As,
   Element::Se,
   Element::Br,
   Element::Kr,
   Element::Rb,
   Element::Sr,
   Element::Y,
   Element::Zr,
   Element::Nb,
   Element::Mo,
   Element::Tc,
   Element::Ru,
   Element::Rh,
   Element::Pd,
   Element::Ag,
   Element::Cd,
   Element::In,
   Element::Sn,
   Element::Sb,
   Element::Te,
   Element::I,
   Element::Xe,
   Element::Cs,
   Element::Ba,
   Element::La,
   Element::Ce,
   Element::Pr,
   Element::Nd,
   Element::Pm,
   Element::Sm,
   Element::Eu,
   Element::Gd,
   Element::Tb,
   Element::Dy,
   Element::Ho,
   Element::Er,
   Element::Tm,
   Element::Yb,
   Element::Lu,
   Element::Hf,
   Element::Ta,
   Element::W,
   Element::Re,
   Element::Os,
   Element::Ir,
   Element::Pt,
   Element::Au,
   Element::Hg,
   Element::Tl,
   Element::Pb,
   Element::Bi,
   Element::Po,
   Element::At,
   Element::Rn,
   Element::Fr,
   Element::Ra,
   Element::Ac,
   Element::Th,
   Element::Pa,
   Element::U,
   Element::Np,
   Element::Pu,
   Element::Am,
   Element::Cm,
   Element::Bk,
   Element::Cf,
   Element::Es,
   Element::Fm,
   Element::Md,
   Element::No,
   Element::Lr,
   Element::Rf,
   Element::Db,
   Element::Sg,
   Element::Bh,
   Element::Hs,
   Element::Mt,
   Element::Uun,
   Element::Uuu,
   Element::Uub
};

char const * const Element::mElementNames[] = {
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

char const * const Element::mAbbreviations[] = {
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

Element::Element(int atomicNo)
{
   if ((atomicNo >= elmNone) && (atomicNo < elmEndOfElements)) {
      mAtomicNumber = atomicNo;
   }
   else {
      printf("Wrong atomic number %d\n", atomicNo);
   }
}

Element::Element()
{
   mAtomicNumber = Element::elmNone;
}

void Element::readAtomicWeights()
{
   try {
      std::ifstream file("AtomicWeights.csv");

      int idx = 0;
      for (CSVIterator loop(file); loop != CSVIterator(); ++loop) { // TODO: check if the first line should be removed
         mAtomicWeight[idx] = std::stof((*loop)[0]);
      }
      ++idx;
   }
   catch (std::exception&) {
      throw 0; //throw new EPQFatalException("Fatal error while attempting to load the atomic weights data file.");
   }
}

int Element::atomicNumberForName(char* name)
{
   for (int i = 0; i < elmEndOfElements + 1; ++i) {
      if ((strcmp(Element::mElementNames[i], name) == 0) || (strcmp(Element::mAbbreviations[i], name) == 0)) { // TODO: make it case insensitive
         return i;
      }
   }
   try {
      return std::stoi(name);
   }
   catch (std::exception&) {
      return Element::elmNone;
   }
}

Element Element::byName(char* name)
{
   int z = atomicNumberForName(name);
   return z == 0 ? Element::None : Element::mAllElements[z - 1];
}

Element Element::byAtomicNumber(int an)
{
   return (an >= 1) && (an < elmEndOfElements - 1) ? mAllElements[an - 1] : Element::None;
}

double Element::getAtomicWeight(int atomicNo)
{
   if (mAtomicWeight == NULL) {
      readAtomicWeights();
   }
   return mAtomicWeight[atomicNo - 1];
}

Element const * Element::allElements()
{
   return mAllElements;
}

//Element* Element::range(Element min, Element max)
//{
//   if (min.getAtomicNumber() <= max.getAtomicNumber()) {
//      throw 0;
//   }
//   return Arrays.copyOfRange(mAllElements, min.getAtomicNumber() - 1, max.getAtomicNumber() - 1);
//}
//
//double Element::meanIonizationPotential(int atomicNo) {
//   try {
//      return MeanIonizationPotential.Berger64.compute(Element.byAtomicNumber(atomicNo));
//   }
//   catch (std::exception& ex) {
//      return MeanIonizationPotential.Sternheimer64.compute(Element.byAtomicNumber(atomicNo));
//   }
//}

int Element::getAtomicNumber() {
   return mAtomicNumber;
}

double Element::getAtomicWeight()
{
   return getAtomicWeight(mAtomicNumber);
}

double Element::getMass()
{
   return ToSI::AMU(getAtomicWeight(mAtomicNumber));
}

char const * Element::toAbbrev()
{
   return mAbbreviations[mAtomicNumber];
}

char const * Element::toAbbrev(int atomicNo)
{
   return mAbbreviations[atomicNo];
}

char const * Element::toString(int el)
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

bool Element::isValid(int atomicNo)
{
   return (atomicNo >= elmH) && (atomicNo < elmEndOfElements);
}

bool Element::isValid()
{
   return (mAtomicNumber >= elmH) && (mAtomicNumber < elmEndOfElements);
}

int Element::compareTo(Element e)
{
   if (mAtomicNumber < e.mAtomicNumber) {
      return -1;
   }
   else {
      return mAtomicNumber == e.mAtomicNumber ? 0 : 1;
   }
}

int Element::hashCode()
{
   // mAtomicNumber is always less than 128 (1<<7). Int has 31 + 1 bits. 31-7
   // = 24
   return mAtomicNumber << 24;
}

bool Element::equals(Element el)
{
   return el.mAtomicNumber == mAtomicNumber;
}

char const * Element::toString()
{
   return mElementNames[mAtomicNumber];
}

double Element::getIonizationEnergy()
{
   int idx = 0;
   try {
      std::ifstream file("IonizationEnergies.csv");

      for (CSVIterator loop(file); loop != CSVIterator(); ++loop) { // TODO: check if the first line should be removed
         if ((*loop)[0][0] == '/' && (*loop)[0][1] == '/') {
            continue;
         }
         if (CSVIterator::IsNaN((*loop)[0])) {
            mIonizationEnergy[idx] = -1.0;
         }
         else {
            mIonizationEnergy[idx] = std::stof((*loop)[0]);
         }
      }
      ++idx;
   }
   catch (std::exception&) {
      throw 0; // throw new EPQFatalException("Fatal error while attempting to load the atomic weights data file.");
   }

   double res = (mAtomicNumber - 1 <= 104) ? mIonizationEnergy[mAtomicNumber - 1] : -1.0;
   if (res == -1.0) {
      throw 0; // new EPQFatalException("The ionization energy is not available for " + toAbbrev());
   }
   return res;
}

Element Element::readResolve()
{
   return Element::byAtomicNumber(mAtomicNumber);
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


