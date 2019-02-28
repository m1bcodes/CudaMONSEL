#ifndef MATERIAL_CUH
#define MATERIAL_CUH

#include "Composition.cuh"

namespace Material
{
   class Material : public Composition::Composition
   {
   protected:
      __device__ Material(Element::Element elms[], int elmsLen, double massFracs[], int massFracsLen, double density, char* name);

   public:
      __device__ Material(double density);
      __device__ Material(Composition comp, double density);
      __device__ Material(Element::Element elm, double density);

      __device__ void setDensity(double den);
      __device__ double getDensity();

      template<typename T>
      __device__ void defineByWeightFraction(LinkedListKV::Node<T, double>* map, double den)
      {
         mDensity = den;
         Composition::defineByWeightFraction(map);
      }

      __device__ void clear();
      __device__ double atomsPerCubicMeter(Element::Element elm);
      __device__ void defineByMaterialFraction(Material mats[], int matsLen, double matFracs[], int matFracsLen);
      __device__ void defineByMaterialFraction(Composition compositions[], int compositionsLen, double matFracs[], int matFracsLen, double density);
      //__device__ char*  descriptiveString(bool normalize);
      __device__ int compareTo(Composition obj);
      __device__ int compareTo(Material obj);
      __device__ void replicate(Material mat);
      __device__ Material clone();
      //__device__ int hashCode();
      __device__ bool equals(Material obj);
      __device__ bool almostEquals(Material other, double tol);

   private:
      double mDensity;
   };
}

#endif
