// file: gov\nist\microanalysis\EPQLibrary\CzyzewskiMottScatteringAngle.cuh
#ifndef _CZYZEWSKI_MOTT_SCATTERING_ANGLE_CUH_
#define _CZYZEWSKI_MOTT_SCATTERING_ANGLE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatter.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatterFactory.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

namespace CzyzewskiMottScatteringAngle
{
   class CzyzewskiMottScatteringAngle : public RandomizedScatterT
   {
   public:
      CzyzewskiMottScatteringAngle(const ElementT& el);

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

   extern const CzyzewskiMottScatteringAngle CMSA1;
   extern const CzyzewskiMottScatteringAngle CMSA2;
   extern const CzyzewskiMottScatteringAngle CMSA3;
   extern const CzyzewskiMottScatteringAngle CMSA4;
   extern const CzyzewskiMottScatteringAngle CMSA5;
   extern const CzyzewskiMottScatteringAngle CMSA6;
   extern const CzyzewskiMottScatteringAngle CMSA7;
   extern const CzyzewskiMottScatteringAngle CMSA8;
   extern const CzyzewskiMottScatteringAngle CMSA9;
   extern const CzyzewskiMottScatteringAngle CMSA10;
   extern const CzyzewskiMottScatteringAngle CMSA11;
   extern const CzyzewskiMottScatteringAngle CMSA12;
   extern const CzyzewskiMottScatteringAngle CMSA13;
   extern const CzyzewskiMottScatteringAngle CMSA14;
   extern const CzyzewskiMottScatteringAngle CMSA15;
   extern const CzyzewskiMottScatteringAngle CMSA16;
   extern const CzyzewskiMottScatteringAngle CMSA17;
   extern const CzyzewskiMottScatteringAngle CMSA18;
   extern const CzyzewskiMottScatteringAngle CMSA19;
   extern const CzyzewskiMottScatteringAngle CMSA20;
   extern const CzyzewskiMottScatteringAngle CMSA21;
   extern const CzyzewskiMottScatteringAngle CMSA22;
   extern const CzyzewskiMottScatteringAngle CMSA23;
   extern const CzyzewskiMottScatteringAngle CMSA24;
   extern const CzyzewskiMottScatteringAngle CMSA25;
   extern const CzyzewskiMottScatteringAngle CMSA26;
   extern const CzyzewskiMottScatteringAngle CMSA27;
   extern const CzyzewskiMottScatteringAngle CMSA28;
   extern const CzyzewskiMottScatteringAngle CMSA29;
   extern const CzyzewskiMottScatteringAngle CMSA30;
   extern const CzyzewskiMottScatteringAngle CMSA31;
   extern const CzyzewskiMottScatteringAngle CMSA32;
   extern const CzyzewskiMottScatteringAngle CMSA33;
   extern const CzyzewskiMottScatteringAngle CMSA34;
   extern const CzyzewskiMottScatteringAngle CMSA35;
   extern const CzyzewskiMottScatteringAngle CMSA36;
   extern const CzyzewskiMottScatteringAngle CMSA37;
   extern const CzyzewskiMottScatteringAngle CMSA38;
   extern const CzyzewskiMottScatteringAngle CMSA39;
   extern const CzyzewskiMottScatteringAngle CMSA40;
   extern const CzyzewskiMottScatteringAngle CMSA41;
   extern const CzyzewskiMottScatteringAngle CMSA42;
   extern const CzyzewskiMottScatteringAngle CMSA43;
   extern const CzyzewskiMottScatteringAngle CMSA44;
   extern const CzyzewskiMottScatteringAngle CMSA45;
   extern const CzyzewskiMottScatteringAngle CMSA46;
   extern const CzyzewskiMottScatteringAngle CMSA47;
   extern const CzyzewskiMottScatteringAngle CMSA48;
   extern const CzyzewskiMottScatteringAngle CMSA49;
   extern const CzyzewskiMottScatteringAngle CMSA50;
   extern const CzyzewskiMottScatteringAngle CMSA51;
   extern const CzyzewskiMottScatteringAngle CMSA52;
   extern const CzyzewskiMottScatteringAngle CMSA53;
   extern const CzyzewskiMottScatteringAngle CMSA54;
   extern const CzyzewskiMottScatteringAngle CMSA55;
   extern const CzyzewskiMottScatteringAngle CMSA56;
   extern const CzyzewskiMottScatteringAngle CMSA57;
   extern const CzyzewskiMottScatteringAngle CMSA58;
   extern const CzyzewskiMottScatteringAngle CMSA59;
   extern const CzyzewskiMottScatteringAngle CMSA60;
   extern const CzyzewskiMottScatteringAngle CMSA61;
   extern const CzyzewskiMottScatteringAngle CMSA62;
   extern const CzyzewskiMottScatteringAngle CMSA63;
   extern const CzyzewskiMottScatteringAngle CMSA64;
   extern const CzyzewskiMottScatteringAngle CMSA65;
   extern const CzyzewskiMottScatteringAngle CMSA66;
   extern const CzyzewskiMottScatteringAngle CMSA67;
   extern const CzyzewskiMottScatteringAngle CMSA68;
   extern const CzyzewskiMottScatteringAngle CMSA69;
   extern const CzyzewskiMottScatteringAngle CMSA70;
   extern const CzyzewskiMottScatteringAngle CMSA71;
   extern const CzyzewskiMottScatteringAngle CMSA72;
   extern const CzyzewskiMottScatteringAngle CMSA73;
   extern const CzyzewskiMottScatteringAngle CMSA74;
   extern const CzyzewskiMottScatteringAngle CMSA75;
   extern const CzyzewskiMottScatteringAngle CMSA76;
   extern const CzyzewskiMottScatteringAngle CMSA77;
   extern const CzyzewskiMottScatteringAngle CMSA78;
   extern const CzyzewskiMottScatteringAngle CMSA79;
   extern const CzyzewskiMottScatteringAngle CMSA80;
   extern const CzyzewskiMottScatteringAngle CMSA81;
   extern const CzyzewskiMottScatteringAngle CMSA82;
   extern const CzyzewskiMottScatteringAngle CMSA83;
   extern const CzyzewskiMottScatteringAngle CMSA84;
   extern const CzyzewskiMottScatteringAngle CMSA85;
   extern const CzyzewskiMottScatteringAngle CMSA86;
   extern const CzyzewskiMottScatteringAngle CMSA87;
   extern const CzyzewskiMottScatteringAngle CMSA88;
   extern const CzyzewskiMottScatteringAngle CMSA89;
   extern const CzyzewskiMottScatteringAngle CMSA90;
   extern const CzyzewskiMottScatteringAngle CMSA91;
   extern const CzyzewskiMottScatteringAngle CMSA92;
   extern const CzyzewskiMottScatteringAngle CMSA93;
   extern const CzyzewskiMottScatteringAngle CMSA94;

   extern CzyzewskiMottScatteringAngle const * mScatter[113];

   class CzyzewskiMottRandomizedScatterFactory : public RandomizedScatterFactoryT
   {
   public:
      CzyzewskiMottRandomizedScatterFactory();
      const RandomizedScatterT& get(const ElementT& elm) const override;

   protected:
      void initializeDefaultStrategy() override;
   };

   extern const RandomizedScatterFactoryT& Factory;
}

#endif