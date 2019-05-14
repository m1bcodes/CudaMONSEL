#ifndef _ATOMIC_SHELL_CUH_
#define _ATOMIC_SHELL_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace AtomicShell
{
   extern const int K;
   extern const int LI;
   extern const int LII;
   extern const int LIII;
   extern const int MI;
   extern const int MII;
   extern const int MIII;
   extern const int MIV;
   extern const int MV;
   extern const int NI;
   extern const int NII;
   extern const int NIII;
   extern const int NIV;
   extern const int NV;
   extern const int NVI;
   extern const int NVII;
   extern const int OI;
   extern const int OII;
   extern const int OIII;
   extern const int OIV;
   extern const int OV;
   extern const int OVI;
   extern const int OVII;
   extern const int OVIII;
   extern const int OIX;
   extern const int PI;
   extern const int PII;
   extern const int PIII;
   extern const int PIV;
   extern const int PV;
   extern const int PVI;
   extern const int PVII;
   extern const int PVIII;
   extern const int PIX;
   extern const int PX;
   extern const int PXI;
   extern const int QI;
   extern const int QII;
   extern const int QIII;
   extern const int QIV;
   extern const int QV;
   extern const int QVI;
   extern const int QVII;
   extern const int QVIII;
   extern const int QIX;
   extern const int QX;
   extern const int QXI;
   extern const int QXII;
   extern const int QXIII;
   extern const int Last;
   extern const int Continuum;
   extern const int NoShell;

   // Shell families
   extern const int NoFamily;
   extern const int KFamily;
   extern const int LFamily;
   extern const int MFamily;
   extern const int NFamily;
   extern const int OFamily;
   extern const int LastFamily;

   class AtomicShell final
   {
   public:
      AtomicShell(const ElementT& el, int shell);
      AtomicShell(const AtomicShell&);

      bool operator=(const AtomicShell&);

      int getGroundStateOccupancy() const;
      StringT getAtomicName() const;
      StringT getSiegbahnName() const;
      StringT getIUPACName() const;
      double getEdgeEnergy() const;
      int getCapacity() const;
      int getFamily() const;
      int getPrincipalQuantumNumber() const;
      double getEnergy() const;
      const ElementT& getElement() const;
      int getShell() const;
      StringT toString() const;
      int compareTo(const AtomicShell& otherShell) const;
      AtomicShell clone() const;
      int hashCode() const;
      bool equals(const AtomicShell& obj) const;
      bool exists() const;
      int getOrbitalAngularMomentum() const;
      double getTotalAngularMomentum() const;

   private:
      const ElementT& mElement;
      const int mShell;
   };

   char const * getAtomicName(int shell);
   char const * getSiegbahnName(int shell);
   int parseSiegahnName(char const * str);
   char const * getIUPACName(int shell);
   int parseIUPACName(char const * str);
   double getEdgeEnergy(const ElementT& el, int shell);
   int getCapacity(int shell);

   bool isValid(int shell);
   int getFamily(int shell);
   int getFirstInFamily(int family);
   char const * getFamilyName(int family);
   int parseFamilyName(char const * s);
   int getPrincipalQuantumNumber(int shell);

   bool exists(const ElementT& elm, int shell);
   int getOrbitalAngularMomentum(int shell);
   double getTotalAngularMomentum(int shell);
   bool electricDipolePermitted(int shell1, int shell2);
   bool electricQuadrupolePermitted(int shell1, int shell2);
}

#endif