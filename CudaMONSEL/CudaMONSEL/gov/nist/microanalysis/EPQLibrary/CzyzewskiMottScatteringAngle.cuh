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
      __host__ __device__ CzyzewskiMottScatteringAngle(const ElementT& el);

      void init(int an);

      StringT toString() const;

      __host__ __device__ double randomScatteringAngle(const double energy, const double rand) const;
      __host__ __device__ double randomScatteringAngle(const double energy) const override;
      __host__ __device__ double totalCrossSection(const double energy) const override;
      __host__ __device__ const ElementT& getElement() const override;

      double meanFreePath(double energy) const;

      template<typename T>
      __device__ void assignMeanFreePath(T* data, const unsigned int len)
      {
         mMeanFreePath.set_data(data, len);
      }

      template<typename T>
      __device__ void assignTotalCrossSection(T* data, const unsigned int len)
      {
         mTotalCrossSection.set_data(data, len);
      }

      template<typename T>
      __device__ void assignCummulativeDFRow(const unsigned int r, T* data, const unsigned int len)
      {
         mCummulativeDF[r].set_data(data, len);
      }

      __host__ __device__ const VectorXf& getMeanFreePath() const;
      __host__ __device__ const VectorXf& getTotalCrossSection() const;
      __host__ __device__ const MatrixXf& getCummulativeDF() const;

   private:
      __host__ __device__ double scatteringAngleForSpecialEnergy(int ei, double rand) const;

      const ElementT& mElement;
      VectorXf mMeanFreePath;
      VectorXf mTotalCrossSection;
      MatrixXf mCummulativeDF;
      const ScreenedRutherfordScatteringAngleT& mRutherford;
   };
   
   class CzyzewskiMottRandomizedScatterFactory : public RandomizedScatterFactoryT
   {
   public:
      __host__ __device__ CzyzewskiMottRandomizedScatterFactory();
      __host__ __device__ const RandomizedScatterT& get(const ElementT& elm) const override;
   };

   extern const RandomizedScatterFactoryT& Factory;

   extern void init();
   extern __global__ void initCuda();
   extern void transferDataToCuda();
}

#endif