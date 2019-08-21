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
      __host__ __device__ explicit NISTMottRS(const ElementT& elm, int method);

      __host__ __device__ double randomScatteringAngle(const double energy) const override;
      __host__ __device__ double totalCrossSection(const double energy) const override;
      __host__ __device__ const ElementT& getElement() const override;

      const char* toString();
      int getMethod() const;
      //void setMethod(int method);

   private:
      //void loadData(int an); // same as NISTMottScatteringAngle

      const ElementT& mElement;
      const VectorXf& mSpwem;
      const MatrixXf& mX1;

      const ScreenedRutherfordScatteringAngleT& mRutherford;
      const BrowningEmpiricalCrossSectionT& mBrowning;

      int method;
      double extrapolateBelowEnergy;

      double MottXSatMinEnergy;
      double sfBrowning;

      StringT name;
   };

   class NISTMottRSFactory : public RandomizedScatterFactoryT
   {
   public:
      __host__ __device__ NISTMottRSFactory(int method);

      __host__ __device__ const RandomizedScatterT& get(const ElementT& elm) const override;

   private:
      int method;
   };

   //extern const double MAX_NISTMOTT;
   //extern const double MIN_NISTMOTT;

   extern const RandomizedScatterFactoryT& Factory;
   extern const RandomizedScatterFactoryT& Factory100;
   extern const RandomizedScatterFactoryT& Factory100Lin;

   __device__ extern const RandomizedScatterFactoryT* dFactory;
   __device__ extern const RandomizedScatterFactoryT* dFactory100;
   __device__ extern const RandomizedScatterFactoryT* dFactory100Lin;

   extern void init();
   extern __global__ void initCuda();
   extern __global__ void initFactory();
}

#endif