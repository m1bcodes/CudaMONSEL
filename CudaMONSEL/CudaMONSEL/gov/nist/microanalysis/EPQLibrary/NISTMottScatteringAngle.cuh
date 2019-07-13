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
      __host__ __device__ explicit NISTMottScatteringAngle(const ElementT& elm);

      StringT toString() const;

      __host__ __device__ const ElementT& getElement() const override;
      __host__ __device__ double totalCrossSection(const double) const override;
      __host__ __device__ double randomScatteringAngle(const double) const override;

      __host__ __device__ const VectorXf& getSpwem() const;
      __host__ __device__ const MatrixXf& getX1() const;

      template<typename T>
      __device__ void copySpwem(const T *dSpwem, const unsigned int size)
      {
         mSpwem.assign(dSpwem, dSpwem + size);
      }

      template<typename T>
      __device__ void copyX1Row(const unsigned int r, const T * dSX1r, const unsigned int size)
      {
         mX1[r].assign(dSX1r, dSX1r + size);
      }

   private:
      void loadData(int an);

      const ElementT& mElement;
      VectorXf mSpwem;
      MatrixXf mX1;
      const ScreenedRutherfordScatteringAngleT& mRutherford;
   };

   //extern const double MAX_NISTMOTT;

   class NISTMottRandomizedScatterFactory : public RandomizedScatterFactoryT
   {
   public:
      __host__ __device__ NISTMottRandomizedScatterFactory();

      __host__ __device__ const RandomizedScatterT& get(const ElementT& elm) const override;
   };

   const NISTMottScatteringAngle* mScatter[];

   extern const RandomizedScatterFactoryT& Factory;
   extern __device__ const RandomizedScatterFactoryT* d_Factory;

   __host__ __device__ extern const NISTMottScatteringAngle& getNISTMSA(int an);

   extern void init();
   extern __global__ void initCuda();
   extern void copyDataToCuda();
   extern __global__ void initFactory();
}

#endif