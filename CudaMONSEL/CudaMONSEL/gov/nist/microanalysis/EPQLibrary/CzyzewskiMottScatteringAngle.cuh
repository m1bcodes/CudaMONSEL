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
      __device__ void copyMeanFreePath(const T* data, const unsigned int len)
      {
         mMeanFreePath.assign(data, data + len);
      }

      template<typename T>
      __device__ void copyTotalCrossSection(const T* data, const unsigned int len)
      {
         mTotalCrossSection.assign(data, data + len);
      }

      template<typename T>
      __device__ void copyCummulativeDFRow(const unsigned int r, const T* data, const unsigned int len)
      {
         mCummulativeDF[r].assign(data, data + len);
      }

      __host__ __device__ const VectorXd& getMeanFreePath() const;
      __host__ __device__ const VectorXd& getTotalCrossSection() const;
      __host__ __device__ const MatrixXd& getCummulativeDF() const;

   private:
      __host__ __device__ double scatteringAngleForSpecialEnergy(int ei, double rand) const;

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
      __host__ __device__ CzyzewskiMottRandomizedScatterFactory();
      __host__ __device__ const RandomizedScatterT& get(const ElementT& elm) const override;
   };

   extern const RandomizedScatterFactoryT& Factory;

   extern void init();
   extern __global__ void initCuda();
   extern void copyDataToCuda();
}

#endif