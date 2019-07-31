#ifndef _TABULATED_INELASTIC_SM_CUH_
#define _TABULATED_INELASTIC_SM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ScatterMechanism.cuh"

namespace TabulatedInelasticSM
{
   typedef ::ScatterMechanism::data_type data_type;

   class TabulatedInelasticSM : public ScatterMechanismT
   {
   public:
      __host__ __device__ TabulatedInelasticSM(const SEmaterialT& mat, int methodSE, const char* tables[], data_type energyOffset);
      __host__ __device__ TabulatedInelasticSM(const SEmaterialT& mat, int methodSE, const char* tables[]);

      __host__ __device__ data_type scatterRate(const ElectronT* pe) override;
      __host__ __device__ ElectronT* scatter(ElectronT* pe) override;
      __host__ __device__ void setMaterial(const MaterialT* mat) override;

      void setMinEgenSE(data_type minEgenSE);
      data_type getMinEgenSE() const;
      int getMethodSE() const;
      void setBranchingRatios();
      void setBranchingRatios(data_type ratios[], int);
      void setEnergyGap(data_type energyGap);
      void setRateMult(data_type rateMult);
      bool isE0fromDispersion() const;
      void setE0fromDispersion(bool e0fromDispersion);

      StringT toString() const;

      __host__ __device__ const NUTableInterpolationT* gettableIIMFP() const;
      __host__ __device__ const NUTableInterpolationT* gettableReducedDeltaE() const;
      __host__ __device__ const NUTableInterpolationT* gettableTheta() const;
      __host__ __device__ const NUTableInterpolationT* gettableSEE0() const;

   private:
      __host__ __device__ void simESEf(data_type Eq, data_type deltaE, data_type r, data_type[2]);
      __host__ __device__ data_type pickBE(data_type Eq, data_type deltaE);

      const int methodSE;
      data_type energyOffset;

      const NUTableInterpolationT* tableIIMFP;
      const NUTableInterpolationT* tableReducedDeltaE;
      const NUTableInterpolationT* tableTheta;
      const NUTableInterpolationT* tableSEE0;

      data_type offsetFermiEnergy;
      data_type energyCBbottom;
      data_type minEgenSE;
      data_type workfunction;
      data_type bandgap;
      data_type energyGap;
      bool defaultRatios;
      MatrixXf cumulativeBranchingProbabilities;

      /*
      * bEref is the energy (relative to conduction band bottom) to which core
      * level binding energies are referenced. This is generally the Fermi energy
      * for metals and 0 for insulators or semiconductors.
      */
      data_type bEref;
      VectorXf kEa; // For convenience, because 1-d
      // tables
      // still require an array for input
      VectorXf coreEnergies;
      VectorXf interpInput;

      // Allowed energy ranges for interpolation table inputs
      VectorXf tableEiDomain;
      VectorXf tableIIMFPEiDomain;
      // Range of allowed SE initial energies on output
      VectorXf energyRangeSE0;

      // temp? IIMFP multiplier
      data_type rateMult;

      /*
      * E0 method. E0 is the energy available to ionize an inner shell. My
      * original method based on the description in Ding & Shimizu's SCANNING
      * article was to assume that deltaE, i.e., all energy transferred in an
      * inelastic event, is available for ionization. Ding et al later modified
      * this procedure. By the plasmon dispersion (I use Penn's form) deltaE can
      * be broken down into a q=0 part and a q != 0 part. The q = 0 part is E0.
      * When E0fromDispersion = false I used deltaE. When it is true I use Ding's
      * later method.
      */
      bool E0fromDispersion;
   };
}

#endif