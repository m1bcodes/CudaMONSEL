#include "Element.cuh"

namespace Element
{
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
      mAtomicNumber = elmNone;
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
         if ((strcmp(mElementNames[i], name) == 0) || (strcmp(mAbbreviations[i], name) == 0)) { // TODO: make it case insensitive
            return i;
         }
      }
      try {
         return std::stoi(name);
      }
      catch (std::exception&) {
         return elmNone;
      }
   }

   Element Element::byName(char* name)
   {
      int z = atomicNumberForName(name);
      return z == 0 ? None : mAllElements[z - 1];
   }

   Element Element::byAtomicNumber(int an)
   {
      return (an >= 1) && (an < elmEndOfElements - 1) ? mAllElements[an - 1] : None;
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

}
