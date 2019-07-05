#ifndef _TABULATED_INELASTIC_SM_CUH_
#define _TABULATED_INELASTIC_SM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ScatterMechanism.cuh"

namespace TabulatedInelasticSM
{
   class TabulatedInelasticSM : public ScatterMechanismT
   {
   public:
      TabulatedInelasticSM(const SEmaterialT& mat, int methodSE, const char* tables[], double energyOffset);
      TabulatedInelasticSM(const SEmaterialT& mat, int methodSE, const char* tables[]);

      double scatterRate(const ElectronT* pe) override;
      ElectronT* scatter(ElectronT* pe) override;
      void setMaterial(const MaterialT* mat) override;

      void setMinEgenSE(double minEgenSE);
      double getMinEgenSE() const;
      int getMethodSE() const;
      void setBranchingRatios();
      void setBranchingRatios(double ratios[], int);
      void setEnergyGap(double energyGap);
      void setRateMult(double rateMult);
      bool isE0fromDispersion() const;
      void setE0fromDispersion(bool e0fromDispersion);

      StringT toString() const;

   private:
      void simESEf(double Eq, double deltaE, double r, double[2]);
      double pickBE(double Eq, double deltaE);

      const int methodSE;
      double energyOffset;

      const NUTableInterpolationT* tableIIMFP;
      const NUTableInterpolationT* tableReducedDeltaE;
      const NUTableInterpolationT* tableTheta;
      const NUTableInterpolationT* tableSEE0;

      double offsetFermiEnergy;
      double energyCBbottom;
      double minEgenSE;
      double workfunction;
      double bandgap;
      double energyGap;
      bool defaultRatios;
      MatrixXd cumulativeBranchingProbabilities;

      /*
      * bEref is the energy (relative to conduction band bottom) to which core
      * level binding energies are referenced. This is generally the Fermi energy
      * for metals and 0 for insulators or semiconductors.
      */
      double bEref;
      VectorXd kEa; // For convenience, because 1-d
      // tables
      // still require an array for input
      VectorXd coreEnergies;
      VectorXd interpInput;

      // Allowed energy ranges for interpolation table inputs
      VectorXd tableEiDomain;
      VectorXd tableIIMFPEiDomain;
      // Range of allowed SE initial energies on output
      VectorXd energyRangeSE0;

      // temp? IIMFP multiplier
      double rateMult;

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