// file: gov\nist\microanalysis\EPQLibrary\CzyzewskiMottCrossSection.cuh

#ifndef _CZYZEWSKI_MOTT_CROSS_SECTION_CUH_
#define _CZYZEWSKI_MOTT_CROSS_SECTION_CUH_

#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

namespace CzyzewskiMottCrossSection
{
   class CzyzewskiMottCrossSection
   {
   public:
      CzyzewskiMottCrossSection(const ElementT& el);
      CzyzewskiMottCrossSection(int an);

      StringT toString() const;
      __host__ __device__ const ElementT& getElement() const;

      double totalCrossSection(const double energy) const;
      double partialCrossSection(double elevation, double azimuth, double energy) const;
      double partialCrossSection(double elevation, double energy) const;
      double meanFreePath(double energy) const;

   private:
      void loadTables(int atomicNo);

      const ElementT& mElement;
      MatrixXd mValues;
   };

   __host__ __device__ double getSpecialEnergy(const int index);
   __host__ __device__ double getSpecialAngle(const int index);
   __host__ __device__ int getEnergyIndex(const double energy);

   extern const int SpecialEnergyCount;
   extern const int SpecialAngleCount;
   extern const double MaxEnergy;
}

#endif