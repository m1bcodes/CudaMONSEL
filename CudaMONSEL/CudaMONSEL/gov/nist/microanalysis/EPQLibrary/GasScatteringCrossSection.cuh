#ifndef _GAS_SCATTERING_CROSS_SECTION_CUH_
#define _GAS_SCATTERING_CROSS_SECTION_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatter.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatterFactory.cuh"

namespace GasScatteringCrossSection
{
   class GasScatteringCrossSection : public RandomizedScatterT
   {
   public:
      GasScatteringCrossSection(const ElementT& elm);
      GasScatteringCrossSection(int an);
      GasScatteringCrossSection(const GasScatteringCrossSection& gscs);

      const RandomizedScatterT& getElasticModel();
      double ratioInelasticOverElastic() const;

      const ElementT& getElement() const override;
      double totalCrossSection(double energy) const override;
      double randomScatteringAngle(double energy) const override;

   private:
      const ElementT& mElement;
      const ScreenedRutherfordScatteringAngleT& mElastic;
   };

   extern const GasScatteringCrossSection GSCS1;
   extern const GasScatteringCrossSection GSCS2;
   extern const GasScatteringCrossSection GSCS3;
   extern const GasScatteringCrossSection GSCS4;
   extern const GasScatteringCrossSection GSCS5;
   extern const GasScatteringCrossSection GSCS6;
   extern const GasScatteringCrossSection GSCS7;
   extern const GasScatteringCrossSection GSCS8;
   extern const GasScatteringCrossSection GSCS9;
   extern const GasScatteringCrossSection GSCS10;
   extern const GasScatteringCrossSection GSCS11;
   extern const GasScatteringCrossSection GSCS12;
   extern const GasScatteringCrossSection GSCS13;
   extern const GasScatteringCrossSection GSCS14;
   extern const GasScatteringCrossSection GSCS15;
   extern const GasScatteringCrossSection GSCS16;
   extern const GasScatteringCrossSection GSCS17;
   extern const GasScatteringCrossSection GSCS18;
   extern const GasScatteringCrossSection GSCS19;
   extern const GasScatteringCrossSection GSCS20;
   extern const GasScatteringCrossSection GSCS21;
   extern const GasScatteringCrossSection GSCS22;
   extern const GasScatteringCrossSection GSCS23;
   extern const GasScatteringCrossSection GSCS24;
   extern const GasScatteringCrossSection GSCS25;
   extern const GasScatteringCrossSection GSCS26;
   extern const GasScatteringCrossSection GSCS27;
   extern const GasScatteringCrossSection GSCS28;
   extern const GasScatteringCrossSection GSCS29;
   extern const GasScatteringCrossSection GSCS30;
   extern const GasScatteringCrossSection GSCS31;
   extern const GasScatteringCrossSection GSCS32;
   extern const GasScatteringCrossSection GSCS33;
   extern const GasScatteringCrossSection GSCS34;
   extern const GasScatteringCrossSection GSCS35;
   extern const GasScatteringCrossSection GSCS36;
   extern const GasScatteringCrossSection GSCS37;
   extern const GasScatteringCrossSection GSCS38;
   extern const GasScatteringCrossSection GSCS39;
   extern const GasScatteringCrossSection GSCS40;
   extern const GasScatteringCrossSection GSCS41;
   extern const GasScatteringCrossSection GSCS42;
   extern const GasScatteringCrossSection GSCS43;
   extern const GasScatteringCrossSection GSCS44;
   extern const GasScatteringCrossSection GSCS45;
   extern const GasScatteringCrossSection GSCS46;
   extern const GasScatteringCrossSection GSCS47;
   extern const GasScatteringCrossSection GSCS48;
   extern const GasScatteringCrossSection GSCS49;
   extern const GasScatteringCrossSection GSCS50;
   extern const GasScatteringCrossSection GSCS51;
   extern const GasScatteringCrossSection GSCS52;
   extern const GasScatteringCrossSection GSCS53;
   extern const GasScatteringCrossSection GSCS54;
   extern const GasScatteringCrossSection GSCS55;
   extern const GasScatteringCrossSection GSCS56;
   extern const GasScatteringCrossSection GSCS57;
   extern const GasScatteringCrossSection GSCS58;
   extern const GasScatteringCrossSection GSCS59;
   extern const GasScatteringCrossSection GSCS60;
   extern const GasScatteringCrossSection GSCS61;
   extern const GasScatteringCrossSection GSCS62;
   extern const GasScatteringCrossSection GSCS63;
   extern const GasScatteringCrossSection GSCS64;
   extern const GasScatteringCrossSection GSCS65;
   extern const GasScatteringCrossSection GSCS66;
   extern const GasScatteringCrossSection GSCS67;
   extern const GasScatteringCrossSection GSCS68;
   extern const GasScatteringCrossSection GSCS69;
   extern const GasScatteringCrossSection GSCS70;
   extern const GasScatteringCrossSection GSCS71;
   extern const GasScatteringCrossSection GSCS72;
   extern const GasScatteringCrossSection GSCS73;
   extern const GasScatteringCrossSection GSCS74;
   extern const GasScatteringCrossSection GSCS75;
   extern const GasScatteringCrossSection GSCS76;
   extern const GasScatteringCrossSection GSCS77;
   extern const GasScatteringCrossSection GSCS78;
   extern const GasScatteringCrossSection GSCS79;
   extern const GasScatteringCrossSection GSCS80;
   extern const GasScatteringCrossSection GSCS81;
   extern const GasScatteringCrossSection GSCS82;
   extern const GasScatteringCrossSection GSCS83;
   extern const GasScatteringCrossSection GSCS84;
   extern const GasScatteringCrossSection GSCS85;
   extern const GasScatteringCrossSection GSCS86;
   extern const GasScatteringCrossSection GSCS87;
   extern const GasScatteringCrossSection GSCS88;
   extern const GasScatteringCrossSection GSCS89;
   extern const GasScatteringCrossSection GSCS90;
   extern const GasScatteringCrossSection GSCS91;
   extern const GasScatteringCrossSection GSCS92;
   extern const GasScatteringCrossSection GSCS93;
   extern const GasScatteringCrossSection GSCS94;
   extern const GasScatteringCrossSection GSCS95;
   extern const GasScatteringCrossSection GSCS96;

   extern GasScatteringCrossSection const * mScatter[113];

   class GasScatteringRandomizedScatterFactory : public RandomizedScatterFactoryT
   {
   protected:
      void initializeDefaultStrategy();

   public:
      GasScatteringRandomizedScatterFactory();

      const RandomizedScatterT& get(const ElementT& elm) const override;
   };

   extern const RandomizedScatterFactoryT& Factory;
}

#endif