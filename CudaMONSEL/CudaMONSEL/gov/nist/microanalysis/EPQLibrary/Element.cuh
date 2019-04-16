#ifndef _ELEMENT_CUH_
#define _ELEMENT_CUH_

#include <stdio.h>
#include <fstream>

#include "gov\nist\microanalysis\Utility\CSVReader.h"
#include "ToSI.cuh"

#include <vector>
#include <unordered_set>

namespace Element
{
   extern float mIonizationEnergy[104]; // Nominal in Joules, IonizationEnergies.csv
   extern float mAtomicWeight[112]; // nominal, in AMU, AtomicWeights.csv

   extern const long long serialVersionUID;

   extern const int elmNone;
   extern const int elmH;
   extern const int elmHe;
   extern const int elmLi;
   extern const int elmBe;
   extern const int elmB;
   extern const int elmC;
   extern const int elmN;
   extern const int elmO;
   extern const int elmF;
   extern const int elmNe;
   extern const int elmNa;
   extern const int elmMg;
   extern const int elmAl;
   extern const int elmSi;
   extern const int elmP;
   extern const int elmS;
   extern const int elmCl;
   extern const int elmAr;
   extern const int elmK;
   extern const int elmCa;
   extern const int elmSc;
   extern const int elmTi;
   extern const int elmV;
   extern const int elmCr;
   extern const int elmMn;
   extern const int elmFe;
   extern const int elmCo;
   extern const int elmNi;
   extern const int elmCu;
   extern const int elmZn;
   extern const int elmGa;
   extern const int elmGe;
   extern const int elmAs;
   extern const int elmSe;
   extern const int elmBr;
   extern const int elmKr;
   extern const int elmRb;
   extern const int elmSr;
   extern const int elmY;
   extern const int elmZr;
   extern const int elmNb;
   extern const int elmMo;
   extern const int elmTc;
   extern const int elmRu;
   extern const int elmRh;
   extern const int elmPd;
   extern const int elmAg;
   extern const int elmCd;
   extern const int elmIn;
   extern const int elmSn;
   extern const int elmSb;
   extern const int elmTe;
   extern const int elmI;
   extern const int elmXe;
   extern const int elmCs;
   extern const int elmBa;
   extern const int elmLa;
   extern const int elmCe;
   extern const int elmPr;
   extern const int elmNd;
   extern const int elmPm;
   extern const int elmSm;
   extern const int elmEu;
   extern const int elmGd;
   extern const int elmTb;
   extern const int elmDy;
   extern const int elmHo;
   extern const int elmEr;
   extern const int elmTm;
   extern const int elmYb;
   extern const int elmLu;
   extern const int elmHf;
   extern const int elmTa;
   extern const int elmW;
   extern const int elmRe;
   extern const int elmOs;
   extern const int elmIr;
   extern const int elmPt;
   extern const int elmAu;
   extern const int elmHg;
   extern const int elmTl;
   extern const int elmPb;
   extern const int elmBi;
   extern const int elmPo;
   extern const int elmAt;
   extern const int elmRn;
   extern const int elmFr;
   extern const int elmRa;
   extern const int elmAc;
   extern const int elmTh;
   extern const int elmPa;
   extern const int elmU;
   extern const int elmNp;
   extern const int elmPu;
   extern const int elmAm;
   extern const int elmCm;
   extern const int elmBk;
   extern const int elmCf;
   extern const int elmEs;
   extern const int elmFm;
   extern const int elmMd;
   extern const int elmNo;
   extern const int elmLr;
   extern const int elmRf;
   extern const int elmDb;
   extern const int elmSg;
   extern const int elmBh;
   extern const int elmHs;
   extern const int elmMt;
   extern const int elmUun;
   extern const int elmUuu;
   extern const int elmUub;
   extern const int elmEndOfElements;
   
   class Element
   {
   public:
      Element(int atomicNo);
      Element(const Element& atomicNo);
      Element();

      Element& operator=(const Element&);

      bool operator==(const Element&) const;
      bool operator<(const Element&) const;

      int getAtomicNumber() const;
      double getAtomicWeight() const;
      double getMass() const;
      char const * toAbbrev() const;

      bool isValid() const;
      int compareTo(const Element& e) const;
      unsigned int hashCode() const;
      bool equals(const Element& el);
      char const * toString();
      double getIonizationEnergy() const;

   private:
      Element readResolve();

      int mAtomicNumber;
   };

   struct HashFcn
   {
      inline unsigned int operator()(const Element& e) const
      {
         return e.hashCode();
      }
   };

   extern const Element None;
   extern const Element H;
   extern const Element He;
   extern const Element Li;
   extern const Element Be;
   extern const Element B;
   extern const Element C;
   extern const Element N;
   extern const Element O;
   extern const Element F;
   extern const Element Ne;
   extern const Element Na;
   extern const Element Mg;
   extern const Element Al;
   extern const Element Si;
   extern const Element P;
   extern const Element S;
   extern const Element Cl;
   extern const Element Ar;
   extern const Element K;
   extern const Element Ca;
   extern const Element Sc;
   extern const Element Ti;
   extern const Element V;
   extern const Element Cr;
   extern const Element Mn;
   extern const Element Fe;
   extern const Element Co;
   extern const Element Ni;
   extern const Element Cu;
   extern const Element Zn;
   extern const Element Ga;
   extern const Element Ge;
   extern const Element As;
   extern const Element Se;
   extern const Element Br;
   extern const Element Kr;
   extern const Element Rb;
   extern const Element Sr;
   extern const Element Y;
   extern const Element Zr;
   extern const Element Nb;
   extern const Element Mo;
   extern const Element Tc;
   extern const Element Ru;
   extern const Element Rh;
   extern const Element Pd;
   extern const Element Ag;
   extern const Element Cd;
   extern const Element In;
   extern const Element Sn;
   extern const Element Sb;
   extern const Element Te;
   extern const Element I;
   extern const Element Xe;
   extern const Element Cs;
   extern const Element Ba;
   extern const Element La;
   extern const Element Ce;
   extern const Element Pr;
   extern const Element Nd;
   extern const Element Pm;
   extern const Element Sm;
   extern const Element Eu;
   extern const Element Gd;
   extern const Element Tb;
   extern const Element Dy;
   extern const Element Ho;
   extern const Element Er;
   extern const Element Tm;
   extern const Element Yb;
   extern const Element Lu;
   extern const Element Hf;
   extern const Element Ta;
   extern const Element W;
   extern const Element Re;
   extern const Element Os;
   extern const Element Ir;
   extern const Element Pt;
   extern const Element Au;
   extern const Element Hg;
   extern const Element Tl;
   extern const Element Pb;
   extern const Element Bi;
   extern const Element Po;
   extern const Element At;
   extern const Element Rn;
   extern const Element Fr;
   extern const Element Ra;
   extern const Element Ac;
   extern const Element Th;
   extern const Element Pa;
   extern const Element U;
   extern const Element Np;
   extern const Element Pu;
   extern const Element Am;
   extern const Element Cm;
   extern const Element Bk;
   extern const Element Cf;
   extern const Element Es;
   extern const Element Fm;
   extern const Element Md;
   extern const Element No;
   extern const Element Lr;
   extern const Element Rf;
   extern const Element Db;
   extern const Element Sg;
   extern const Element Bh;
   extern const Element Hs;
   extern const Element Mt;
   extern const Element Uun;
   extern const Element Uuu;
   extern const Element Uub;

   void readAtomicWeights();
   void readIonizationEnergy();
   int atomicNumberForName(char const * name);
   Element byName(char const * name);
   Element byAtomicNumber(int an);
   double getAtomicWeight(int atomicNo);
   Element const * allElements();
   std::vector<Element> range(const Element& min, const Element& max);
   //double meanIonizationPotential(int atomicNo);

   char const * toAbbrev(int atomicNo);
   char const * toString(int el);
   bool isValid(int atomicNo);

   std::vector<std::string> getListOfAbbreviations(const Element& minEl, const Element& maxEl);
   std::vector<std::string> getListOfElements(const Element& minEl, const Element& maxEl);

   void InitializeElements();

   typedef std::unordered_set<Element, HashFcn> UnorderedSetT;
   typedef std::unordered_set<Element, HashFcn> OrderedSetT;
}
#endif
