#ifndef _NIST_MOTT_RS_CUH_
#define _NIST_MOTT_RS_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatter.cuh"
#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatterFactory.cuh"

namespace NISTMottRS
{
   class NISTMottRS : public RandomizedScatterT
   {
   public:
      NISTMottRS(const ElementT& elm, int method);

      double randomScatteringAngle(double energy) const override;
      double totalCrossSection(double energy) const override;
      const ElementT& getElement() const override;

      StringT toString() const;
      int getMethod() const;
      //void setMethod(int method);

   private:
      //void loadData(int an); // same as NISTMottScatteringAngle

      int method;
      double extrapolateBelowEnergy;
      double MottXSatMinEnergy;

      const ElementT& mElement;
      const VectorXd& mSpwem;
      const MatrixXd& mX1;

      const ScreenedRutherfordScatteringAngleT& mRutherford;
      const BrowningEmpiricalCrossSectionT& mBrowning;
      double sfBrowning;
   };

   class NISTMottRSFactory : public RandomizedScatterFactoryT
   {
   public:
      NISTMottRSFactory(int method);

      const RandomizedScatterT& get(const ElementT& elm) const override;
      void initializeDefaultStrategy() override;

   private:
      int method;
   };

   extern const double MAX_NISTMOTT;
   extern const double MIN_NISTMOTT;

   extern const RandomizedScatterFactoryT& Factory;
   extern const RandomizedScatterFactoryT& Factory100;
   extern const RandomizedScatterFactoryT& Factory100Lin;

   extern const NISTMottRS& getNMRS1(int an);
   extern const NISTMottRS& getNMRS2(int an);
   extern const NISTMottRS& getNMRS3(int an);
}

#endif