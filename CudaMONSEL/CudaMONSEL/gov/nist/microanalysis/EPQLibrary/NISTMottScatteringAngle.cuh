// file: gov\nist\microanalysis\EPQLibrary\NISTMottScatteringAngle.cuh

#ifndef _NIST_MOTT_SCATTERING_ANGLE_CUH_
#define _NIST_MOTT_SCATTERING_ANGLE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatter.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatterFactory.cuh"

namespace NISTMottScatteringAngle
{
   extern const int SPWEM_LEN;
   extern const int X1_LEN;
   extern const double DL50;
   extern const double PARAM;

   class NISTMottScatteringAngle : public RandomizedScatterT
   {
   public:
      NISTMottScatteringAngle(const ElementT& elm);
      StringT toString() const;

      const ElementT& getElement() const override;
      double totalCrossSection(double energy) const override;
      double randomScatteringAngle(double energy) const override;

   private:
      const ElementT& mElement;
      VectorXd mSpwem;
      MatrixXd mX1;
      const ScreenedRutherfordScatteringAngleT& mRutherford;
   };

   extern const double MAX_NISTMOTT;

   extern const NISTMottScatteringAngle* mScatter[113];

   class NISTMottRandomizedScatterFactory : public RandomizedScatterFactoryT
   {
   public:
      NISTMottRandomizedScatterFactory();

      const RandomizedScatterT& get(const ElementT& elm) const override;

   protected:
      void initializeDefaultStrategy() override;
   };

   extern const RandomizedScatterFactoryT& Factory;
}

#endif