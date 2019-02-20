#ifndef _ELEMENT_CUH_
#define _ELEMENT_CUH_

#include <stdio.h>
#include <fstream>

#include "../Utility/CSVReader.h"
#include "ToSI.cuh"

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

   float mIonizationEnergy[104]; // Nominal in Joules, IonizationEnergies.csv
   float mAtomicWeight[112]; // nominal, in AMU, AtomicWeights.csv

   extern const long long serialVersionUID;
}
#endif
