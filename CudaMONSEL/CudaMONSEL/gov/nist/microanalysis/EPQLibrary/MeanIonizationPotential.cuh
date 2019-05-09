#ifndef _MEAN_IONIZATION_POTENTIAL_CUH_
#define _MEAN_IONIZATION_POTENTIAL_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AlgorithmClass.cuh"

namespace MeanIonizationPotential
{
   class MeanIonizationPotential : public AlgorithmClassT
   {
   public:
      void initializeDefaultStrategy() override;
      AlgorithmClassT const * const * getAllImplementations() const;

      double computeLn(const CompositionT& comp) const;
      virtual double compute(const ElementT& el) const = 0;

   protected:
      MeanIonizationPotential(StringT name, const ReferenceT& reference);
   };

   class Sternheimer64MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      Sternheimer64MeanIonizationPotential();
      double compute(const ElementT& el) const override;
   };

   class BergerAndSeltzerCITZAFMeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      BergerAndSeltzerCITZAFMeanIonizationPotential();
      double compute(const ElementT& el) const override;
   };

   class Bloch33MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      Bloch33MeanIonizationPotential();
      double compute(const ElementT& el) const override;
   };

   class Wilson41MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      Wilson41MeanIonizationPotential();
      double compute(const ElementT& el) const override;
   };

   class Springer67MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      Springer67MeanIonizationPotential();
      double compute(const ElementT& el) const override;
   };

   class Heinrich70MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      Heinrich70MeanIonizationPotential();
      double compute(const ElementT& el) const override;
   };

   class Duncumb69MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      Duncumb69MeanIonizationPotential();
      double compute(const ElementT& el) const override;
   };

   class Zeller75MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      Zeller75MeanIonizationPotential();
      double compute(const ElementT& el) const override;
   };

   class Berger64MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      Berger64MeanIonizationPotential();
      double compute(const ElementT& el) const override;

   private:
      void readTabulatedValues();
      static VectorXd mMeasured; // nominal, in Joules
   };

   class Berger83MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      Berger83MeanIonizationPotential();
      double compute(const ElementT& el) const override;

   private:
      void readTabulatedValues();
      static VectorXd mMeasured; // nominal, in Joules
   };

   extern const MeanIonizationPotential& Berger64;
   extern const MeanIonizationPotential& Berger83;
   extern const MeanIonizationPotential& Bloch33;
   extern const MeanIonizationPotential& Duncumb69;
   extern const MeanIonizationPotential& BergerAndSeltzerCITZAF;
   extern const MeanIonizationPotential& Heinrich70;
   extern const MeanIonizationPotential& Springer67;
   extern const MeanIonizationPotential& Sternheimer64;
   extern const MeanIonizationPotential& Wilson41;
   extern const MeanIonizationPotential& Zeller75;
}

#endif