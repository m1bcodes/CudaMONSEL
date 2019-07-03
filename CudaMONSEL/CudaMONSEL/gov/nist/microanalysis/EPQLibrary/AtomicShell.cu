#include "AtomicShell.cuh"
#include "AlgorithmUser.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\EdgeEnergy.cuh"

#include "gov\nist\microanalysis\Utility\CSVReader.h"

namespace AtomicShell
{
   const int K = 0;
   const int LI = 1;
   const int LII = 2;
   const int LIII = 3;
   const int MI = 4;
   const int MII = 5;
   const int MIII = 6;
   const int MIV = 7;
   const int MV = 8;
   const int NI = 9;
   const int NII = 10;
   const int NIII = 11;
   const int NIV = 12;
   const int NV = 13;
   const int NVI = 14;
   const int NVII = 15;
   const int OI = 16;
   const int OII = 17;
   const int OIII = 18;
   const int OIV = 19;
   const int OV = 20;
   const int OVI = 21;
   const int OVII = 22;
   const int OVIII = 23;
   const int OIX = 24;
   const int PI = 25;
   const int PII = 26;
   const int PIII = 27;
   const int PIV = 28;
   const int PV = 29;
   const int PVI = 30;
   const int PVII = 31;
   const int PVIII = 32;
   const int PIX = 33;
   const int PX = 34;
   const int PXI = 35;
   const int QI = 36;
   const int QII = 37;
   const int QIII = 38;
   const int QIV = 39;
   const int QV = 40;
   const int QVI = 41;
   const int QVII = 42;
   const int QVIII = 43;
   const int QIX = 44;
   const int QX = 45;
   const int QXI = 46;
   const int QXII = 47;
   const int QXIII = 48;
   const int Last = 49;
   const int Continuum = 1000;
   const int NoShell = 1001;

   // Shell families
   const int NoFamily = 1999;
   const int KFamily = 2000;
   const int LFamily = 2001;
   const int MFamily = 2002;
   const int NFamily = 2003;
   const int OFamily = 2004;
   const int LastFamily = 2005;

   static char const * const mAtomicNames[] = {
      "1S",
      "2S",
      "2P1/2",
      "2P3/2",
      "3S",
      "3P1/2",
      "3P3/2",
      "3D3/2",
      "3D5/2",
      "4S",
      "4P1/2",
      "4P3/2",
      "4D3/2",
      "4D5/2",
      "4F5/2",
      "4F7/2",
      "5S",
      "5P1/2",
      "5P3/2",
      "5D3/2",
      "5D5/2",
      "5F5/2",
      "5F7/2",
      "5G7/2",
      "5G9/2",
      "6S",
      "6P1/2",
      "6P3/2",
      "6D3/2",
      "6D5/2",
      "6F5/2",
      "6F7/2",
      "6G7/2",
      "6G9/2",
      "6H9/2",
      "6H11/2",
      "7S",
      "7P1/2",
      "7P3/2",
      "7D3/2",
      "7D5/2",
      "7F5/2",
      "7F7/2",
      "7G7/2",
      "7G92/",
      "7H9/2",
      "7H11/2",
      "7I11/2",
      "7I13/2"
   };

   static char const * const mSiegbahnNames[] = {
      "K",
      "LI",
      "LII",
      "LIII",
      "MI",
      "MII",
      "MIII",
      "MIV",
      "MV",
      "NI",
      "NII",
      "NIII",
      "NIV",
      "NV",
      "NVI",
      "NVII",
      "OI",
      "OII",
      "OIII",
      "OIV",
      "OV",
      "OVI",
      "OVII",
      "OVIII",
      "OIX",
      "PI",
      "PII",
      "PIII",
      "PIV",
      "PV",
      "PVI",
      "PVII",
      "PVIII",
      "PIX",
      "PX",
      "PXI",
      "QI",
      "QII",
      "QIII",
      "QIV",
      "QV",
      "QVI",
      "QVII",
      "QVIII",
      "QIX",
      "QX",
      "QXI",
      "QXII",
      "QXIII"
   };

   static char const * const mIUPACNames[] = {
      "K",
      "L1",
      "L2",
      "L3",
      "M1",
      "M2",
      "M3",
      "M4",
      "M5",
      "N1",
      "N2",
      "N3",
      "N4",
      "N5",
      "N6",
      "N7",
      "O1",
      "O2",
      "O3",
      "O4",
      "O5",
      "O6",
      "O7",
      "O8",
      "O9",
      "P1",
      "P2",
      "P3",
      "P4",
      "P5",
      "P6",
      "P7",
      "P8",
      "P9",
      "P10",
      "P11",
      "Q1",
      "Q2",
      "Q3",
      "Q4",
      "Q5",
      "Q6",
      "Q7",
      "Q8",
      "Q9",
      "Q10",
      "Q11",
      "Q12",
      "Q13"
   };

   static const int mCapacity[] = {
      2,
      2,
      2,
      4,
      2,
      2,
      4,
      4,
      6,
      2,
      2,
      4,
      4,
      6,
      6,
      8,
      2,
      2,
      4,
      4,
      6,
      6,
      8,
      8,
      10,
      2,
      2,
      4,
      4,
      6,
      6,
      8,
      8,
      10,
      10,
      12,
      2,
      2,
      4,
      4,
      6,
      6,
      8,
      8,
      10,
      10,
      12,
      12,
      14
   };

   static const int mOrbitalAngularMomentum[] = {
      0, // N=1
      0,
      1,
      1, // N = 2
      0,
      1,
      1,
      2,
      2, // N = 3
      0,
      1,
      1,
      2,
      2,
      3,
      3, // N = 4
      0,
      1,
      1,
      2,
      2,
      3,
      3,
      4,
      4, // N = 5
      0,
      1,
      1,
      2,
      2,
      3,
      3,
      4,
      4,
      5,
      5, // N =6
      0,
      1,
      1,
      2,
      2,
      3,
      3,
      4,
      4,
      5,
      5,
      6,
      6
      // N =7
   };

   static double mTotalAngularMomentum[] = {
      0.5, // N=1
      0.5,
      0.5,
      1.5, // N = 2
      0.5,
      0.5,
      1.5,
      1.5,
      2.5, // N = 3
      0.5,
      0.5,
      1.5,
      1.5,
      2.5,
      2.5,
      3.5, // N = 4
      0.5,
      0.5,
      1.5,
      1.5,
      2.5,
      2.5,
      3.5,
      3.5,
      4.5, // N = 5
      0.5,
      0.5,
      1.5,
      1.5,
      2.5,
      2.5,
      3.5,
      3.5,
      4.5,
      4.5,
      5.5, // N =6
      0.5,
      0.5,
      1.5,
      1.5,
      2.5,
      2.5,
      3.5,
      3.5,
      4.5,
      4.5,
      5.5,
      5.5,
      6.5
      // N =7
   };

   static MatrixXi mOccupancy;
   static void loadGroundStateOccupancy()
   {
      if (mOccupancy.empty()) {
         mOccupancy.resize(102);
         char filepath[] = ".\\gov\\nist\\microanalysis\\EPQLibrary\\ElectronConfig.csv";
         try {
            printf("Reading: %s\n", filepath);
            std::ifstream file(filepath);
            if (!file.good()) throw 0;
            int idx = 0;
            for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
               if ((*loop)[0][0] == '#') { // check if the first line should be removed
                  continue;
               }
               int z = std::stoi((*loop)[0]);
               int sum = 0;
               mOccupancy[idx].reserve((*loop).size());
               for (int k = 1; k < (*loop).size(); ++k) {
                  int cur = std::stoi((*loop)[k]);
                  mOccupancy[idx].push_back(cur);
                  sum += cur;
               }
               if (sum != z) {
                  printf("getGroundStateOccupancy(): sum not equal to z");
               }
               ++idx;
            }
            file.close();
         }
         catch (std::exception&) {
            printf("Fatal error while attempting to load the ionization data file: %s.\n", filepath);
         }
      }
   }

   int AtomicShell::getGroundStateOccupancy() const
   {
      auto occ = mOccupancy[mElement.getAtomicNumber() - 1];
      return occ.size() > mShell ? occ[mShell] : 0;
   }

   bool isLineFamily(int f)
   {
      return (f == KFamily) || (f == LFamily) || (f == MFamily) || (f == NFamily) || (f == OFamily);
   }

   AtomicShell::AtomicShell(const ElementT& el, int shell) : mElement(el), mShell(shell)
   {
      if (!isValid(shell)) printf("Shell = %d\n", shell);
      if (sizeof(mAtomicNames) / sizeof(*mAtomicNames) != Last) printf("AtomicShell::AtomicShell: mAtomicNames wrong size %d\n", sizeof(mAtomicNames) / sizeof(*mAtomicNames));
      if (sizeof(mSiegbahnNames) / sizeof(*mSiegbahnNames) != Last) printf("AtomicShell::AtomicShell: mSiegbahnNames wrong size %d\n", sizeof(mSiegbahnNames) / sizeof(*mSiegbahnNames));
      if (sizeof(mCapacity) / sizeof(*mCapacity) != Last) printf("AtomicShell::AtomicShell: mCapacity wrong size %d\n", sizeof(mCapacity) / sizeof(*mCapacity));
      if (sizeof(mOrbitalAngularMomentum) / sizeof(*mOrbitalAngularMomentum) != Last) printf("AtomicShell::AtomicShell: mOrbitalAngularMomentum wrong size %d\n", sizeof(mOrbitalAngularMomentum) / sizeof(*mOrbitalAngularMomentum));
      loadGroundStateOccupancy();
   }

   char const * getAtomicName(int shell)
   {
      if ((shell >= K) && (shell < sizeof(mAtomicNames) / sizeof(*mAtomicNames))) {
         return mAtomicNames[shell];
      }
      else if (shell == Continuum) {
         return "Continuum";
      }
      else {
         return "Unknown";
      }
   }

   StringT AtomicShell::getAtomicName() const
   {
      StringT elm(mElement.toAbbrev());
      StringT shell(::AtomicShell::getAtomicName(mShell));
      return elm + " " + shell;
   }

   char const * getSiegbahnName(int shell)
   {
      if ((shell >= K) && (shell < Last)) {
         return mSiegbahnNames[shell];
      }
      else if (shell == Continuum) {
         return "Continuum";
      }
      else {
         return "Unknown";
      }
   }

   int parseSiegahnName(char const * str)
   {
      for (int sh = K; sh < Last; ++sh)
         if (str == mSiegbahnNames[sh])
            return sh;
      if (str == "Continuum")
         return Continuum;
      return NoShell;

   }

   StringT AtomicShell::getSiegbahnName() const
   {
      StringT elm(mElement.toAbbrev());
      StringT shell(::AtomicShell::getSiegbahnName(mShell));
      return elm + " " + shell;
   }

   char const * getIUPACName(int shell)
   {
      if ((shell >= K) && (shell < Last))
         return mIUPACNames[shell];
      else if (shell == Continuum)
         return "Continuum";
      else
         return "Unknown";
   }

   StringT AtomicShell::getIUPACName() const
   {
      StringT elm(mElement.toAbbrev());
      StringT shell(::AtomicShell::getIUPACName(mShell));
      return elm + " " + shell;
   }

   int parseIUPACName(char const * s)
   {
      std::string str(s);
      for (int i = 0; i < str.size(); i++) {
         str.at(i) = tolower(str.at(i));
      }
      for (int sh = 0; sh < Last; ++sh)
         if (str == mIUPACNames[sh])
            return sh;
      if (str == "Continuum")
         return Continuum;
      return NoShell;
   }

   double getEdgeEnergy(const ElementT& el, int shell)
   {
      return AlgorithmUser::getDefaultEdgeEnergy().compute(AtomicShell(el, shell));
   }

   double AtomicShell::getEdgeEnergy() const
   {
      return AlgorithmUser::getDefaultEdgeEnergy().compute(*this);
   }

   int getCapacity(int shell)
   {
      if (!isValid(shell) || (shell == Continuum)) {
         printf("AtomicShell::getCapacity wrong condition");
      }
      return mCapacity[shell];
   }

   int AtomicShell::getCapacity() const
   {
      return ::AtomicShell::getCapacity(mShell);
   }

   bool isValid(int shell)
   {
      return (shell >= K) && ((shell < Last) || (shell == Continuum));
   }

   int getFamily(int shell)
   {
      if (shell == K)
         return KFamily;
      else if (shell <= LIII)
         return LFamily;
      else if (shell <= MV)
         return MFamily;
      else if (shell <= NVII)
         return NFamily;
      else if (shell <= OIX)
         return OFamily;
      else
         return NoFamily;
   }

   int getFirstInFamily(int family)
   {
      switch (family) {
      case KFamily:
         return K;
      case LFamily:
         return LI;
      case MFamily:
         return MI;
      case NFamily:
         return NI;
      case OFamily:
         return OI;
      default:
         printf("Unknown family in getFirstInFamily");
      }
      return NoShell;
   }

   int getLastInFamily(int family)
   {
      switch (family) {
      case KFamily:
         return K;
      case LFamily:
         return LIII;
      case MFamily:
         return MV;
      case NFamily:
         return NVII;
      case OFamily:
         return OIX;
      default:
         printf("Unknown family in getFirstInFamily");
      }
      return NoShell;
   }

   char const * getFamilyName(int family)
   {
      switch (family) {
      case KFamily:
         return "K";
      case LFamily:
         return "L";
      case MFamily:
         return "M";
      case NFamily:
         return "N";
      case OFamily:
         return "O";
      }
      return "None";
   }

   int parseFamilyName(char const * s)
   {
      StringT str(s);
      for (int i = 0; i < str.size(); i++) {
         str[i] = toupper(str.at(i));
      }
      if (str == "K")
         return KFamily;
      else if (str == "L")
         return LFamily;
      else if (str == "M")
         return MFamily;
      else if (str == "N")
         return NFamily;
      else if (str == "O")
         return OFamily;
      return NoFamily;
   }

   int AtomicShell::getFamily() const
   {
      return ::AtomicShell::getFamily(mShell);
   }

   int getPrincipalQuantumNumber(int shell)
   {
      if (!isValid(shell)) {
         printf("AtomicShell::getPrincipalQuantumNumber invalid argument: %d\n", shell);
         return -1;
      }
      if (shell == K)
         return 1;
      else if (shell <= LIII)
         return 2;
      else if (shell <= MV)
         return 3;
      else if (shell <= NVII)
         return 4;
      else if (shell <= OIX)
         return 5;
      else if (shell <= PXI)
         return 6;
      else if (shell <= QXIII)
         return 7;
      printf("Unexpected shell in AtomicShell.getPrincipleQuantumNumber: %d" + shell);
      return -1;
   }

   int AtomicShell::getPrincipalQuantumNumber() const
   {
      return ::AtomicShell::getPrincipalQuantumNumber(mShell);
   }

   double AtomicShell::getEnergy() const
   {
      return getGroundStateOccupancy() > 0 ? getEdgeEnergy() + mElement.getIonizationEnergy() : NAN;
   }

   const ElementT& AtomicShell::getElement() const
   {
      return mElement;
   }

   int AtomicShell::getShell() const
   {
      return mShell;
   }

   StringT AtomicShell::toString() const
   {
      StringT elm(mElement.toAbbrev());
      StringT shell(::AtomicShell::getIUPACName(mShell));
      StringT ret = elm + " " + shell;
      return ret;
   }

   AtomicShell parseString(char const * const s)
   {
      std::string str(s);
      int p = str.find(' ');
      return AtomicShell(Element::byName(str.substr(0, p).c_str()), parseIUPACName(str.substr(p + 1).c_str()));
   }

   int AtomicShell::compareTo(const AtomicShell& otherShell) const
   {
      int an = mElement.getAtomicNumber();
      int ano = otherShell.mElement.getAtomicNumber();
      if (an < ano)
         return -1;
      else if (an == ano) {
         if (mShell < otherShell.mShell)
            return -1;
         else if (mShell == otherShell.mShell)
            return 0;
         else
            return 1;
      }
      else
         return 1;
   }

   AtomicShell AtomicShell::clone() const
   {
      return AtomicShell(mElement, mShell);
   }

   int AtomicShell::hashCode() const
   {
      // mShell < 1024 = 1<<10.
      return mElement.hashCode() ^ (mShell << 14);
   }

   bool AtomicShell::equals(const AtomicShell& sh) const
   {
      return (sh.mElement == mElement) && (sh.mShell == mShell);
   }

   bool exists(const Element::Element& elm, int shell)
   {
      return getEdgeEnergy(elm, shell) > 0.0;
   }

   bool AtomicShell::exists() const
   {
      return getEdgeEnergy() > 0.0;
   }

   int getOrbitalAngularMomentum(int shell)
   {
      return mOrbitalAngularMomentum[shell]; // S->0, P->1, D->2, F->3
   }

   int AtomicShell::getOrbitalAngularMomentum() const
   {
      return mOrbitalAngularMomentum[mShell]; // S->0, P->1, D->2, F->3
   }

   double AtomicShell::getTotalAngularMomentum() const
   {
      return mTotalAngularMomentum[mShell];
   }

   double getTotalAngularMomentum(int shell)
   {
      return mTotalAngularMomentum[shell];
   }

   bool electricDipolePermitted(int shell1, int shell2)
   {
      { // deltaJ=0,+1,-1 but no 0->0
         double deltaJ = abs(mTotalAngularMomentum[shell1] - mTotalAngularMomentum[shell2]);
         if (deltaJ > 1.0)
            return false;
         if (deltaJ != 0.0 && deltaJ != 1.0) {
            printf("AtomicShell::electricDipolePermitted wrong deltaJ: %lf\n", deltaJ);
         }
         // J=0->J=0
         // if((mTotalAngularMomentum[shell1] == 0.0) &&
         // (mTotalAngularMomentum[shell2] == 0.0))
         // return false;
      }
      // deltaL=+1,-1
      return abs(mOrbitalAngularMomentum[shell1] - mOrbitalAngularMomentum[shell2]) == 1.0;
   }

   bool electricQuadrupolePermitted(int shell1, int shell2)
   {
      { // deltaJ=0,+1,-1,+2,-2 but no 0->0
         double deltaJ = abs(mTotalAngularMomentum[shell1] - mTotalAngularMomentum[shell2]);
         if (deltaJ > 2.0)
            return false;
         if (deltaJ != 0.0 && deltaJ != 1.0 && deltaJ != 2.0) {
            printf("AtomicShell::electricQuadrupolePermitted wrong deltaJ: %lf\n", deltaJ);
         }
         if ((mTotalAngularMomentum[shell1] == 0.5) && (mTotalAngularMomentum[shell2] == 0.5))
            return false;
      }
      // deltaL=0,+2,-2
      double deltaL = abs(mOrbitalAngularMomentum[shell1] - mOrbitalAngularMomentum[shell2]);
      return (deltaL == 0.0) || (deltaL == 2.0);
   }
}