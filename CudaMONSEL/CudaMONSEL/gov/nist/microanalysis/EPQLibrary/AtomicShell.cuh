#ifndef _ATOMIC_SHELL_CUH_
#define _ATOMIC_SHELL_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace AtomicShell
{
   class AtomicShell final
   {
   public:
      typedef std::vector<std::vector<int>> OccupancyT;

   public:
      AtomicShell(const ElementT& el, int shell);
      AtomicShell(const AtomicShell&);

      bool operator=(const AtomicShell&);

      int getGroundStateOccupancy() const;
      char const * getAtomicName() const;
      char const * getSiegbahnName() const;
      char const * getIUPACName() const;
      int getCapacity() const;
      int getFamily() const;
      int getPrincipalQuantumNumber() const;
      double getEnergy();
      const ElementT& getElement() const;
      int getShell() const;
      char const * const toString() const;
      int compareTo(const AtomicShell& otherShell) const;
      AtomicShell clone() const;
      int hashCode() const;
      bool equals(const AtomicShell& obj) const;
      int getOrbitalAngularMomentum() const;
      double getTotalAngularMomentum() const;

   private:
      void loadGroundStateOccupancy();
      const ElementT& mElement;
      const int mShell;

      OccupancyT mOccupancy;
   };

   char const * getAtomicName(int shell);
   char const * getSiegbahnName(int shell);
   int parseSiegahnName(char const * str);
   char const * getIUPACName(int shell);
   int parseIUPACName(char const * str);
   int getCapacity(int shell);

   bool isValid(int shell);
   int getFamily(int shell);
   int getFirstInFamily(int family);
   char const * getFamilyName(int family);
   int parseFamilyName(char const * s);

   bool exists(const ElementT& elm, int shell);
   bool exists();
   int getOrbitalAngularMomentum(int shell);
   double getTotalAngularMomentum(int shell);
   bool electricDipolePermitted(int shell1, int shell2);
   bool electricQuadrupolePermitted(int shell1, int shell2);
}

#endif