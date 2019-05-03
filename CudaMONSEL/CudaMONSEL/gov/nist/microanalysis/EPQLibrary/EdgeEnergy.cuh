#ifndef _EDGE_ENERGY_CUH_
#define _EDGE_ENERGY_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AlgorithmClass.cuh"

namespace EdgeEnergy
{
   class EdgeEnergy : public AlgorithmClassT
   {
   protected:
      void initializeDefaultStrategy() override;

      EdgeEnergy(StringT name, StringT ref);
      EdgeEnergy(StringT name, const ReferenceT& ref);

   public:
      AlgorithmClassT const * const * getAllImplementations() const override;

      virtual double compute(const AtomicShellT& shell) const = 0;
      virtual bool isSupported(const AtomicShellT& shell) const = 0;
   };

   class DiracHartreeSlaterIonizationEnergies : public EdgeEnergy
   {
   public:
      DiracHartreeSlaterIonizationEnergies();

      double compute(const AtomicShellT& shell) const override;
      bool isSupported(const AtomicShellT& shell) const override;
   };

   class NISTEdgeEnergy : public EdgeEnergy
   {
   public:
      NISTEdgeEnergy();
      double compute(const AtomicShellT& shell) const;
      bool isSupported(const AtomicShellT& shell) const;
   };

   class ChantlerEdgeEnergy : public EdgeEnergy
   {
   public:
      ChantlerEdgeEnergy();
      double compute(const AtomicShellT& shell) const;
      bool isSupported(const AtomicShellT& shell) const;
   };


   class WernishEdgeEnergy: public EdgeEnergy
   {
   public:
      WernishEdgeEnergy();
      double compute(const AtomicShellT& shell) const override;
      bool isSupported(const AtomicShellT& shell) const override;
   };

   class DTSAEdgeEnergy : public EdgeEnergy
   {
   public:
      DTSAEdgeEnergy();
      double compute(const AtomicShellT& shell) const override;
      bool isSupported(const AtomicShellT& shell) const override;
   };

   class SuperSetEdgeEnergy : public EdgeEnergy
   {
   public:
      SuperSetEdgeEnergy();
      double compute(const AtomicShellT& shell) const override;
      bool isSupported(const AtomicShellT& shell) const override;
   };

   extern const EdgeEnergy& DHSIonizationEnergy;
   extern const EdgeEnergy& NISTxrtdb;
   extern const EdgeEnergy& Chantler2005;
   extern const EdgeEnergy& Wernish84;
   extern const EdgeEnergy& DTSA;
   extern const EdgeEnergy& SuperSet;
}

#endif