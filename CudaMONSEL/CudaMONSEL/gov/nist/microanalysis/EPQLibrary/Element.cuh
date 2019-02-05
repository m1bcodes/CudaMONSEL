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

   const Element None;
   const Element H;
   const Element He;
   const Element Li;
   const Element Be;
   const Element B;
   const Element C;
   const Element N;
   const Element O;
   const Element F;
   const Element Ne;
   const Element Na;
   const Element Mg;
   const Element Al;
   const Element Si;
   const Element P;
   const Element S;
   const Element Cl;
   const Element Ar;
   const Element K;
   const Element Ca;
   const Element Sc;
   const Element Ti;
   const Element V;
   const Element Cr;
   const Element Mn;
   const Element Fe;
   const Element Co;
   const Element Ni;
   const Element Cu;
   const Element Zn;
   const Element Ga;
   const Element Ge;
   const Element As;
   const Element Se;
   const Element Br;
   const Element Kr;
   const Element Rb;
   const Element Sr;
   const Element Y;
   const Element Zr;
   const Element Nb;
   const Element Mo;
   const Element Tc;
   const Element Ru;
   const Element Rh;
   const Element Pd;
   const Element Ag;
   const Element Cd;
   const Element In;
   const Element Sn;
   const Element Sb;
   const Element Te;
   const Element I;
   const Element Xe;
   const Element Cs;
   const Element Ba;
   const Element La;
   const Element Ce;
   const Element Pr;
   const Element Nd;
   const Element Pm;
   const Element Sm;
   const Element Eu;
   const Element Gd;
   const Element Tb;
   const Element Dy;
   const Element Ho;
   const Element Er;
   const Element Tm;
   const Element Yb;
   const Element Lu;
   const Element Hf;
   const Element Ta;
   const Element W;
   const Element Re;
   const Element Os;
   const Element Ir;
   const Element Pt;
   const Element Au;
   const Element Hg;
   const Element Tl;
   const Element Pb;
   const Element Bi;
   const Element Po;
   const Element At;
   const Element Rn;
   const Element Fr;
   const Element Ra;
   const Element Ac;
   const Element Th;
   const Element Pa;
   const Element U;
   const Element Np;
   const Element Pu;
   const Element Am;
   const Element Cm;
   const Element Bk;
   const Element Cf;
   const Element Es;
   const Element Fm;
   const Element Md;
   const Element No;
   const Element Lr;
   const Element Rf;
   const Element Db;
   const Element Sg;
   const Element Bh;
   const Element Hs;
   const Element Mt;
   const Element Uun;
   const Element Uuu;
   const Element Uub;

   const long long serialVersionUID;

   const Element mAllElements[elmEndOfElements + 1];
   char const * const mElementNames[elmEndOfElements + 1];
   char const * const mAbbreviations[elmEndOfElements + 1];

   float mIonizationEnergy[104]; // Nominal in Joules, IonizationEnergies.csv
   float mAtomicWeight[112]; // nominal, in AMU, AtomicWeights.csv

}
#endif
