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

      __device__ int getAtomicNumber();
      __device__ double getAtomicWeight();
      __device__ double getMass();
      __device__ char const * toAbbrev();

      __device__ bool isValid();
      __device__ int compareTo(Element e);
      __device__ int hashCode();
      __device__ bool equals(Element el);
      __device__ char const * toString();
      __device__  double getIonizationEnergy();

   private:
      int mAtomicNumber;

      __device__ Element readResolve();
   };

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

   __host__ void InitializeElements();
}
#endif
