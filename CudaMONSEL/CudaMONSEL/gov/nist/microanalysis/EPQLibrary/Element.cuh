#ifndef _ELEMENT_CUH_
#define _ELEMENT_CUH_

#include <stdio.h>
#include <fstream>

#include "../Utility/CSVReader.h"
#include "ToSI.cuh"

namespace Element
{
   extern __device__ float mIonizationEnergy[104]; // Nominal in Joules, IonizationEnergies.csv
   extern __device__ float mAtomicWeight[112]; // nominal, in AMU, AtomicWeights.csv

   extern __device__ const long long serialVersionUID;

   extern __device__ const int elmNone;
   extern __device__ const int elmH;
   extern __device__ const int elmHe;
   extern __device__ const int elmLi;
   extern __device__ const int elmBe;
   extern __device__ const int elmB;
   extern __device__ const int elmC;
   extern __device__ const int elmN;
   extern __device__ const int elmO;
   extern __device__ const int elmF;
   extern __device__ const int elmNe;
   extern __device__ const int elmNa;
   extern __device__ const int elmMg;
   extern __device__ const int elmAl;
   extern __device__ const int elmSi;
   extern __device__ const int elmP;
   extern __device__ const int elmS;
   extern __device__ const int elmCl;
   extern __device__ const int elmAr;
   extern __device__ const int elmK;
   extern __device__ const int elmCa;
   extern __device__ const int elmSc;
   extern __device__ const int elmTi;
   extern __device__ const int elmV;
   extern __device__ const int elmCr;
   extern __device__ const int elmMn;
   extern __device__ const int elmFe;
   extern __device__ const int elmCo;
   extern __device__ const int elmNi;
   extern __device__ const int elmCu;
   extern __device__ const int elmZn;
   extern __device__ const int elmGa;
   extern __device__ const int elmGe;
   extern __device__ const int elmAs;
   extern __device__ const int elmSe;
   extern __device__ const int elmBr;
   extern __device__ const int elmKr;
   extern __device__ const int elmRb;
   extern __device__ const int elmSr;
   extern __device__ const int elmY;
   extern __device__ const int elmZr;
   extern __device__ const int elmNb;
   extern __device__ const int elmMo;
   extern __device__ const int elmTc;
   extern __device__ const int elmRu;
   extern __device__ const int elmRh;
   extern __device__ const int elmPd;
   extern __device__ const int elmAg;
   extern __device__ const int elmCd;
   extern __device__ const int elmIn;
   extern __device__ const int elmSn;
   extern __device__ const int elmSb;
   extern __device__ const int elmTe;
   extern __device__ const int elmI;
   extern __device__ const int elmXe;
   extern __device__ const int elmCs;
   extern __device__ const int elmBa;
   extern __device__ const int elmLa;
   extern __device__ const int elmCe;
   extern __device__ const int elmPr;
   extern __device__ const int elmNd;
   extern __device__ const int elmPm;
   extern __device__ const int elmSm;
   extern __device__ const int elmEu;
   extern __device__ const int elmGd;
   extern __device__ const int elmTb;
   extern __device__ const int elmDy;
   extern __device__ const int elmHo;
   extern __device__ const int elmEr;
   extern __device__ const int elmTm;
   extern __device__ const int elmYb;
   extern __device__ const int elmLu;
   extern __device__ const int elmHf;
   extern __device__ const int elmTa;
   extern __device__ const int elmW;
   extern __device__ const int elmRe;
   extern __device__ const int elmOs;
   extern __device__ const int elmIr;
   extern __device__ const int elmPt;
   extern __device__ const int elmAu;
   extern __device__ const int elmHg;
   extern __device__ const int elmTl;
   extern __device__ const int elmPb;
   extern __device__ const int elmBi;
   extern __device__ const int elmPo;
   extern __device__ const int elmAt;
   extern __device__ const int elmRn;
   extern __device__ const int elmFr;
   extern __device__ const int elmRa;
   extern __device__ const int elmAc;
   extern __device__ const int elmTh;
   extern __device__ const int elmPa;
   extern __device__ const int elmU;
   extern __device__ const int elmNp;
   extern __device__ const int elmPu;
   extern __device__ const int elmAm;
   extern __device__ const int elmCm;
   extern __device__ const int elmBk;
   extern __device__ const int elmCf;
   extern __device__ const int elmEs;
   extern __device__ const int elmFm;
   extern __device__ const int elmMd;
   extern __device__ const int elmNo;
   extern __device__ const int elmLr;
   extern __device__ const int elmRf;
   extern __device__ const int elmDb;
   extern __device__ const int elmSg;
   extern __device__ const int elmBh;
   extern __device__ const int elmHs;
   extern __device__ const int elmMt;
   extern __device__ const int elmUun;
   extern __device__ const int elmUuu;
   extern __device__ const int elmUub;
   extern __device__ const int elmEndOfElements;
   
   class Element
   {
   public:
      __host__ __device__ Element(int atomicNo);
      __host__ __device__ Element();

      __host__ __device__ bool operator==(const Element&);

      __device__ int getAtomicNumber();
      __device__ double getAtomicWeight();
      __device__ double getMass();
      __device__ char const * toAbbrev();

      __device__ bool isValid();
      __device__ int compareTo(Element e);
      __device__ int hashCode();
      __device__ bool equals(const Element& el);
      __device__ char const * toString();
      __device__  double getIonizationEnergy();

   private:
      __device__ Element readResolve();

      int mAtomicNumber;
   };

   extern __device__ const Element None;
   extern __device__ const Element H;
   extern __device__ const Element He;
   extern __device__ const Element Li;
   extern __device__ const Element Be;
   extern __device__ const Element B;
   extern __device__ const Element C;
   extern __device__ const Element N;
   extern __device__ const Element O;
   extern __device__ const Element F;
   extern __device__ const Element Ne;
   extern __device__ const Element Na;
   extern __device__ const Element Mg;
   extern __device__ const Element Al;
   extern __device__ const Element Si;
   extern __device__ const Element P;
   extern __device__ const Element S;
   extern __device__ const Element Cl;
   extern __device__ const Element Ar;
   extern __device__ const Element K;
   extern __device__ const Element Ca;
   extern __device__ const Element Sc;
   extern __device__ const Element Ti;
   extern __device__ const Element V;
   extern __device__ const Element Cr;
   extern __device__ const Element Mn;
   extern __device__ const Element Fe;
   extern __device__ const Element Co;
   extern __device__ const Element Ni;
   extern __device__ const Element Cu;
   extern __device__ const Element Zn;
   extern __device__ const Element Ga;
   extern __device__ const Element Ge;
   extern __device__ const Element As;
   extern __device__ const Element Se;
   extern __device__ const Element Br;
   extern __device__ const Element Kr;
   extern __device__ const Element Rb;
   extern __device__ const Element Sr;
   extern __device__ const Element Y;
   extern __device__ const Element Zr;
   extern __device__ const Element Nb;
   extern __device__ const Element Mo;
   extern __device__ const Element Tc;
   extern __device__ const Element Ru;
   extern __device__ const Element Rh;
   extern __device__ const Element Pd;
   extern __device__ const Element Ag;
   extern __device__ const Element Cd;
   extern __device__ const Element In;
   extern __device__ const Element Sn;
   extern __device__ const Element Sb;
   extern __device__ const Element Te;
   extern __device__ const Element I;
   extern __device__ const Element Xe;
   extern __device__ const Element Cs;
   extern __device__ const Element Ba;
   extern __device__ const Element La;
   extern __device__ const Element Ce;
   extern __device__ const Element Pr;
   extern __device__ const Element Nd;
   extern __device__ const Element Pm;
   extern __device__ const Element Sm;
   extern __device__ const Element Eu;
   extern __device__ const Element Gd;
   extern __device__ const Element Tb;
   extern __device__ const Element Dy;
   extern __device__ const Element Ho;
   extern __device__ const Element Er;
   extern __device__ const Element Tm;
   extern __device__ const Element Yb;
   extern __device__ const Element Lu;
   extern __device__ const Element Hf;
   extern __device__ const Element Ta;
   extern __device__ const Element W;
   extern __device__ const Element Re;
   extern __device__ const Element Os;
   extern __device__ const Element Ir;
   extern __device__ const Element Pt;
   extern __device__ const Element Au;
   extern __device__ const Element Hg;
   extern __device__ const Element Tl;
   extern __device__ const Element Pb;
   extern __device__ const Element Bi;
   extern __device__ const Element Po;
   extern __device__ const Element At;
   extern __device__ const Element Rn;
   extern __device__ const Element Fr;
   extern __device__ const Element Ra;
   extern __device__ const Element Ac;
   extern __device__ const Element Th;
   extern __device__ const Element Pa;
   extern __device__ const Element U;
   extern __device__ const Element Np;
   extern __device__ const Element Pu;
   extern __device__ const Element Am;
   extern __device__ const Element Cm;
   extern __device__ const Element Bk;
   extern __device__ const Element Cf;
   extern __device__ const Element Es;
   extern __device__ const Element Fm;
   extern __device__ const Element Md;
   extern __device__ const Element No;
   extern __device__ const Element Lr;
   extern __device__ const Element Rf;
   extern __device__ const Element Db;
   extern __device__ const Element Sg;
   extern __device__ const Element Bh;
   extern __device__ const Element Hs;
   extern __device__ const Element Mt;
   extern __device__ const Element Uun;
   extern __device__ const Element Uuu;
   extern __device__ const Element Uub;

   __host__ void readAtomicWeights();
   __host__ void readIonizationEnergy();
   __device__ int atomicNumberForName(char* name);
   __device__ Element byName(char* name);
   __device__ Element byAtomicNumber(int an);
   __device__ double getAtomicWeight(int atomicNo);
   __device__ Element const * allElements();
   __device__ Element* range(Element min, Element max);
   //__device__ double meanIonizationPotential(int atomicNo);

   __device__ char const * toAbbrev(int atomicNo);
   __device__ char const * toString(int el);
   __device__ bool isValid(int atomicNo);

   //__device__ char* const getListOfAbbreviations(Element minEl, Element maxEl);

   //__device__ bool AreEqual(Element& e1, Element& e2);
   __host__ void InitializeElements();
}
#endif
