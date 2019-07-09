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
      __host__ __device__ virtual double compute(const ElementT& el) const = 0;

   protected:
      __host__ __device__ MeanIonizationPotential(StringT name, const ReferenceT& reference);
   };

   class Sternheimer64MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Sternheimer64MeanIonizationPotential();
      __host__ __device__ double compute(const ElementT& el) const override;
   };
   extern __host__ __device__ double computeSternheimer64(const ElementT& el);

   class BergerAndSeltzerCITZAFMeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ BergerAndSeltzerCITZAFMeanIonizationPotential();
      __host__ __device__ double compute(const ElementT& el) const override;
   };

   class Bloch33MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Bloch33MeanIonizationPotential();
      __host__ __device__ double compute(const ElementT& el) const override;
   };

   class Wilson41MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Wilson41MeanIonizationPotential();
      __host__ __device__ double compute(const ElementT& el) const override;
   };
   extern __host__ __device__ double computeWilson41(const ElementT& el);

   class Springer67MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Springer67MeanIonizationPotential();
      __host__ __device__ double compute(const ElementT& el) const override;
   };

   class Heinrich70MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Heinrich70MeanIonizationPotential();
      __host__ __device__ double compute(const ElementT& el) const override;
   };

   class Duncumb69MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Duncumb69MeanIonizationPotential();
      __host__ __device__ double compute(const ElementT& el) const override;
   };

   class Zeller75MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Zeller75MeanIonizationPotential();
      __host__ __device__ double compute(const ElementT& el) const override;
   };

   class Berger64MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Berger64MeanIonizationPotential();
      __host__ __device__ double compute(const ElementT& el) const override;
      void readTabulatedValues();
      //__device__ void copyData(double *data);

   private:
      VectorXd mMeasured; // nominal, in Joules
   };

   class Berger83MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Berger83MeanIonizationPotential();
      __host__ __device__ double compute(const ElementT& el) const override;
      void readTabulatedValues();
      //__device__ void copyData(double *data);

   private:
      VectorXd mMeasured; // nominal, in Joules
   };

   extern Berger64MeanIonizationPotential &Berger64;
   extern Berger83MeanIonizationPotential &Berger83;
   extern const MeanIonizationPotential &Bloch33;
   extern const MeanIonizationPotential &Duncumb69;
   extern const MeanIonizationPotential &BergerAndSeltzerCITZAF;
   extern const MeanIonizationPotential &Heinrich70;
   extern const MeanIonizationPotential &Springer67;
   extern const MeanIonizationPotential &Sternheimer64;
   extern const MeanIonizationPotential &Wilson41;
   extern const MeanIonizationPotential &Zeller75;
}

#endif