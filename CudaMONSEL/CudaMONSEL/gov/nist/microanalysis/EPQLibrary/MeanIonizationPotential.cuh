#ifndef _MEAN_IONIZATION_POTENTIAL_CUH_
#define _MEAN_IONIZATION_POTENTIAL_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AlgorithmClass.cuh"

namespace MeanIonizationPotential
{
   class MeanIonizationPotential : public AlgorithmClassT
   {
      __host__ __device__ void initializeDefaultStrategy() override;

   public:
      AlgorithmClassT const * const * getAllImplementations() const;

      float computeLn(const CompositionT& comp) const;
      __host__ __device__ virtual float compute(const ElementT& el) const = 0;

   protected:
      __host__ __device__ MeanIonizationPotential(StringT name, const ReferenceT& reference);
   };

   class Sternheimer64MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Sternheimer64MeanIonizationPotential();
      __host__ __device__ float compute(const ElementT& el) const override;
   };
   extern __host__ __device__ float computeSternheimer64(const ElementT& el);

   class BergerAndSeltzerCITZAFMeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ BergerAndSeltzerCITZAFMeanIonizationPotential();
      __host__ __device__ float compute(const ElementT& el) const override;
   };

   class Bloch33MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Bloch33MeanIonizationPotential();
      __host__ __device__ float compute(const ElementT& el) const override;
   };

   class Wilson41MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Wilson41MeanIonizationPotential();
      __host__ __device__ float compute(const ElementT& el) const override;
   };
   extern __host__ __device__ float computeWilson41(const ElementT& el);

   class Springer67MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Springer67MeanIonizationPotential();
      __host__ __device__ float compute(const ElementT& el) const override;
   };

   class Heinrich70MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Heinrich70MeanIonizationPotential();
      __host__ __device__ float compute(const ElementT& el) const override;
   };

   class Duncumb69MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Duncumb69MeanIonizationPotential();
      __host__ __device__ float compute(const ElementT& el) const override;
   };

   class Zeller75MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Zeller75MeanIonizationPotential();
      __host__ __device__ float compute(const ElementT& el) const override;
   };

   class Berger64MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Berger64MeanIonizationPotential();
      __host__ __device__ float compute(const ElementT& el) const override;
      void readTabulatedValues();
      __host__ __device__ const VectorXf& getData() const;
      template<typename T>
      __device__ void assignData(T* data, const unsigned int len)
      {
         mMeasured.set_data(data, len);
      }

   private:
      VectorXf mMeasured; // nominal, in Joules
   };

   class Berger83MeanIonizationPotential : public MeanIonizationPotential
   {
   public:
      __host__ __device__ Berger83MeanIonizationPotential();
      __host__ __device__ float compute(const ElementT& el) const override;
      void readTabulatedValues();
      __host__ __device__ const VectorXf& getData() const;
      template<typename T>
      __device__ void assignData(T* data, const unsigned int len)
      {
         mMeasured.set_data(data, len);
      }

   private:
      VectorXf mMeasured; // nominal, in Joules
   };

   extern Berger64MeanIonizationPotential& Berger64;
   extern Berger83MeanIonizationPotential& Berger83;
   extern const MeanIonizationPotential& Bloch33;
   extern const MeanIonizationPotential& Duncumb69;
   extern const MeanIonizationPotential& BergerAndSeltzerCITZAF;
   extern const MeanIonizationPotential& Heinrich70;
   extern const MeanIonizationPotential& Springer67;
   extern const MeanIonizationPotential& Sternheimer64;
   extern const MeanIonizationPotential& Wilson41;
   extern const MeanIonizationPotential& Zeller75;

   extern __device__ Berger64MeanIonizationPotential* d_Berger64;
   extern __device__ Berger83MeanIonizationPotential* d_Berger83;
   extern __device__ const MeanIonizationPotential* d_Bloch33;
   extern __device__ const MeanIonizationPotential* d_Duncumb69;
   extern __device__ const MeanIonizationPotential* d_BergerAndSeltzerCITZAF;
   extern __device__ const MeanIonizationPotential* d_Heinrich70;
   extern __device__ const MeanIonizationPotential* d_Springer67;
   extern __device__ const MeanIonizationPotential* d_Sternheimer64;
   extern __device__ const MeanIonizationPotential* d_Wilson41;
   extern __device__ const MeanIonizationPotential* d_Zeller75;

   extern __global__ void initCuda();
   extern void transferDataToCuda();
}

#endif