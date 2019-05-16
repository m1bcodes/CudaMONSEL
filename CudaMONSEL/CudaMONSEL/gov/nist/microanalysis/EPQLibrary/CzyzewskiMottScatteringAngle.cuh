// file: gov\nist\microanalysis\EPQLibrary\CzyzewskiMottScatteringAngle.cuh
#ifndef _CZYZEWSKI_MOTT_SCATTERING_ANGLE_CUH_
#define _CZYZEWSKI_MOTT_SCATTERING_ANGLE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatter.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatterFactory.cuh"

namespace CzyzewskiMottScatteringAngle
{
   class CzyzewskiMottScatteringAngle : public RandomizedScatterT
   {
   public:
      CzyzewskiMottScatteringAngle(const ElementT& el);

      void init(int an);

      StringT toString() const;

      double randomScatteringAngle(double energy, double rand) const;

      double randomScatteringAngle(double energy) const override;
      double totalCrossSection(double energy) const override;
      const ElementT& getElement() const override;

      double meanFreePath(double energy) const;

   private:
      double scatteringAngleForSpecialEnergy(int ei, double rand) const;

   private:
      const ElementT& mElement;
      VectorXd mMeanFreePath;
      VectorXd mTotalCrossSection;
      MatrixXd mCummulativeDF;
      const ScreenedRutherfordScatteringAngleT& mRutherford;
   };
   
   class CzyzewskiMottRandomizedScatterFactory : public RandomizedScatterFactoryT
   {
   public:
      CzyzewskiMottRandomizedScatterFactory();
      const RandomizedScatterT& get(const ElementT& elm) const override;

   protected:
      void initializeDefaultStrategy() override;
   };

   extern const RandomizedScatterFactoryT& Factory;

   extern const CzyzewskiMottScatteringAngle& getCMSA(int an);

   extern void init();
}

#endif