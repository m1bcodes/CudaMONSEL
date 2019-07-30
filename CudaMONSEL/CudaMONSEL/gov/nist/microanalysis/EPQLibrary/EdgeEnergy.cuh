#ifndef _EDGE_ENERGY_CUH_
#define _EDGE_ENERGY_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AlgorithmClass.cuh"

namespace EdgeEnergy
{
   class EdgeEnergy : public AlgorithmClassT
   {
   protected:
      __host__ __device__ void initializeDefaultStrategy() override;

      EdgeEnergy(StringT name, StringT ref);
      EdgeEnergy(StringT name, const ReferenceT& ref);

   public:
      AlgorithmClassT const * const * getAllImplementations() const override;
      virtual float compute(const AtomicShellT& shell) const = 0;
      virtual bool isSupported(const AtomicShellT& shell) const = 0;
   };

   class DiracHartreeSlaterIonizationEnergies : public EdgeEnergy
   {
   public:
      DiracHartreeSlaterIonizationEnergies();
      float compute(const AtomicShellT& shell) const override;
      bool isSupported(const AtomicShellT& shell) const override;
      static void loadxionUis();

   private:
      static MatrixXd Uis; // nominally [100][9]
   };

   class NISTEdgeEnergy : public EdgeEnergy
   {
   public:
      NISTEdgeEnergy();
      float compute(const AtomicShellT& shell) const override;
      bool isSupported(const AtomicShellT& shell) const override;
      static void loadNISTxrtdb();

   private:
      static MatrixXd NISTEdgeEnergies;
   };

   class ChantlerEdgeEnergy : public EdgeEnergy
   {
   public:
      ChantlerEdgeEnergy();
      float compute(const AtomicShellT& shell) const override;
      bool isSupported(const AtomicShellT& shell) const override;
      static void loadFFastEdgeDB();

   private:
      static MatrixXd ChantlerEdgeEnergies;
   };

   class WernishEdgeEnergy: public EdgeEnergy
   {
   public:
      WernishEdgeEnergy();
      float compute(const AtomicShellT& shell) const override;
      bool isSupported(const AtomicShellT& shell) const override;
   };

   class DTSAEdgeEnergy : public EdgeEnergy
   {
   public:
      DTSAEdgeEnergy();
      float compute(const AtomicShellT& shell) const override;
      bool isSupported(const AtomicShellT& shell) const override;
      static void loadEdgeEnergies();

   private:
      static MatrixXd DTSAEdgeEnergies;
   };

   class SuperSetEdgeEnergy : public EdgeEnergy
   {
   public:
      SuperSetEdgeEnergy();
      float compute(const AtomicShellT& shell) const override;
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