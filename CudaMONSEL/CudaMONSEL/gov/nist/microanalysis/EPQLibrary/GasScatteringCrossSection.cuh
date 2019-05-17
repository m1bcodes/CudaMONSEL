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
      //GasScatteringCrossSection(const GasScatteringCrossSection& gscs);

      const RandomizedScatterT& getElasticModel();
      double ratioInelasticOverElastic() const;

      const ElementT& getElement() const override;
      double totalCrossSection(double energy) const override;
      double randomScatteringAngle(double energy) const override;

   private:
      const ElementT& mElement;
      const ScreenedRutherfordScatteringAngleT& mElastic;
   };
   
   class GasScatteringRandomizedScatterFactory : public RandomizedScatterFactoryT
   {
   protected:
      void initializeDefaultStrategy();

   public:
      GasScatteringRandomizedScatterFactory();

      const RandomizedScatterT& get(const ElementT& elm) const override;
   };

   extern const RandomizedScatterFactoryT& Factory;

   extern void init();
}

#endif