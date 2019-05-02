// file: gov\nist\microanalysis\EPQLibrary\ScreenedRutherfordScatteringAngle.cuh

#ifndef _SCREENED_RUTHERFORD_SCATTERING_ANGLE_CUH_
#define _SCREENED_RUTHERFORD_SCATTERING_ANGLE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatter.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatterFactory.cuh"

namespace ScreenedRutherfordScatteringAngle
{
   class ScreenedRutherfordScatteringAngle : public RandomizedScatterT
   {
   public:
      ScreenedRutherfordScatteringAngle(const ElementT& elm);
      ScreenedRutherfordScatteringAngle(int an);
      ScreenedRutherfordScatteringAngle(const ScreenedRutherfordScatteringAngle& other);

      StringT toString() const;

      const ElementT& getElement() const override;
      double totalCrossSection(double energy) const override;
      double randomScatteringAngle(double energy) const override;

   private:
      const ElementT& mElement;
   };

   extern const ScreenedRutherfordScatteringAngle SRSA1;
   extern const ScreenedRutherfordScatteringAngle SRSA2;
   extern const ScreenedRutherfordScatteringAngle SRSA3;
   extern const ScreenedRutherfordScatteringAngle SRSA4;
   extern const ScreenedRutherfordScatteringAngle SRSA5;
   extern const ScreenedRutherfordScatteringAngle SRSA6;
   extern const ScreenedRutherfordScatteringAngle SRSA7;
   extern const ScreenedRutherfordScatteringAngle SRSA8;
   extern const ScreenedRutherfordScatteringAngle SRSA9;
   extern const ScreenedRutherfordScatteringAngle SRSA10;
   extern const ScreenedRutherfordScatteringAngle SRSA11;
   extern const ScreenedRutherfordScatteringAngle SRSA12;
   extern const ScreenedRutherfordScatteringAngle SRSA13;
   extern const ScreenedRutherfordScatteringAngle SRSA14;
   extern const ScreenedRutherfordScatteringAngle SRSA15;
   extern const ScreenedRutherfordScatteringAngle SRSA16;
   extern const ScreenedRutherfordScatteringAngle SRSA17;
   extern const ScreenedRutherfordScatteringAngle SRSA18;
   extern const ScreenedRutherfordScatteringAngle SRSA19;
   extern const ScreenedRutherfordScatteringAngle SRSA20;
   extern const ScreenedRutherfordScatteringAngle SRSA21;
   extern const ScreenedRutherfordScatteringAngle SRSA22;
   extern const ScreenedRutherfordScatteringAngle SRSA23;
   extern const ScreenedRutherfordScatteringAngle SRSA24;
   extern const ScreenedRutherfordScatteringAngle SRSA25;
   extern const ScreenedRutherfordScatteringAngle SRSA26;
   extern const ScreenedRutherfordScatteringAngle SRSA27;
   extern const ScreenedRutherfordScatteringAngle SRSA28;
   extern const ScreenedRutherfordScatteringAngle SRSA29;
   extern const ScreenedRutherfordScatteringAngle SRSA30;
   extern const ScreenedRutherfordScatteringAngle SRSA31;
   extern const ScreenedRutherfordScatteringAngle SRSA32;
   extern const ScreenedRutherfordScatteringAngle SRSA33;
   extern const ScreenedRutherfordScatteringAngle SRSA34;
   extern const ScreenedRutherfordScatteringAngle SRSA35;
   extern const ScreenedRutherfordScatteringAngle SRSA36;
   extern const ScreenedRutherfordScatteringAngle SRSA37;
   extern const ScreenedRutherfordScatteringAngle SRSA38;
   extern const ScreenedRutherfordScatteringAngle SRSA39;
   extern const ScreenedRutherfordScatteringAngle SRSA40;
   extern const ScreenedRutherfordScatteringAngle SRSA41;
   extern const ScreenedRutherfordScatteringAngle SRSA42;
   extern const ScreenedRutherfordScatteringAngle SRSA43;
   extern const ScreenedRutherfordScatteringAngle SRSA44;
   extern const ScreenedRutherfordScatteringAngle SRSA45;
   extern const ScreenedRutherfordScatteringAngle SRSA46;
   extern const ScreenedRutherfordScatteringAngle SRSA47;
   extern const ScreenedRutherfordScatteringAngle SRSA48;
   extern const ScreenedRutherfordScatteringAngle SRSA49;
   extern const ScreenedRutherfordScatteringAngle SRSA50;
   extern const ScreenedRutherfordScatteringAngle SRSA51;
   extern const ScreenedRutherfordScatteringAngle SRSA52;
   extern const ScreenedRutherfordScatteringAngle SRSA53;
   extern const ScreenedRutherfordScatteringAngle SRSA54;
   extern const ScreenedRutherfordScatteringAngle SRSA55;
   extern const ScreenedRutherfordScatteringAngle SRSA56;
   extern const ScreenedRutherfordScatteringAngle SRSA57;
   extern const ScreenedRutherfordScatteringAngle SRSA58;
   extern const ScreenedRutherfordScatteringAngle SRSA59;
   extern const ScreenedRutherfordScatteringAngle SRSA60;
   extern const ScreenedRutherfordScatteringAngle SRSA61;
   extern const ScreenedRutherfordScatteringAngle SRSA62;
   extern const ScreenedRutherfordScatteringAngle SRSA63;
   extern const ScreenedRutherfordScatteringAngle SRSA64;
   extern const ScreenedRutherfordScatteringAngle SRSA65;
   extern const ScreenedRutherfordScatteringAngle SRSA66;
   extern const ScreenedRutherfordScatteringAngle SRSA67;
   extern const ScreenedRutherfordScatteringAngle SRSA68;
   extern const ScreenedRutherfordScatteringAngle SRSA69;
   extern const ScreenedRutherfordScatteringAngle SRSA70;
   extern const ScreenedRutherfordScatteringAngle SRSA71;
   extern const ScreenedRutherfordScatteringAngle SRSA72;
   extern const ScreenedRutherfordScatteringAngle SRSA73;
   extern const ScreenedRutherfordScatteringAngle SRSA74;
   extern const ScreenedRutherfordScatteringAngle SRSA75;
   extern const ScreenedRutherfordScatteringAngle SRSA76;
   extern const ScreenedRutherfordScatteringAngle SRSA77;
   extern const ScreenedRutherfordScatteringAngle SRSA78;
   extern const ScreenedRutherfordScatteringAngle SRSA79;
   extern const ScreenedRutherfordScatteringAngle SRSA80;
   extern const ScreenedRutherfordScatteringAngle SRSA81;
   extern const ScreenedRutherfordScatteringAngle SRSA82;
   extern const ScreenedRutherfordScatteringAngle SRSA83;
   extern const ScreenedRutherfordScatteringAngle SRSA84;
   extern const ScreenedRutherfordScatteringAngle SRSA85;
   extern const ScreenedRutherfordScatteringAngle SRSA86;
   extern const ScreenedRutherfordScatteringAngle SRSA87;
   extern const ScreenedRutherfordScatteringAngle SRSA88;
   extern const ScreenedRutherfordScatteringAngle SRSA89;
   extern const ScreenedRutherfordScatteringAngle SRSA90;
   extern const ScreenedRutherfordScatteringAngle SRSA91;
   extern const ScreenedRutherfordScatteringAngle SRSA92;
   extern const ScreenedRutherfordScatteringAngle SRSA93;
   extern const ScreenedRutherfordScatteringAngle SRSA94;
   extern const ScreenedRutherfordScatteringAngle SRSA95;
   extern const ScreenedRutherfordScatteringAngle SRSA96;

   extern ScreenedRutherfordScatteringAngle const * mScatter[113];

   class ScreenedRutherfordRandomizedScatterFactory : public RandomizedScatterFactoryT
   {
   public:
      ScreenedRutherfordRandomizedScatterFactory();

      const RandomizedScatterT& get(const ElementT& elm) const override;

   protected:
      void initializeDefaultStrategy() override;
   };

   extern const RandomizedScatterFactoryT& Factory;
}

#endif