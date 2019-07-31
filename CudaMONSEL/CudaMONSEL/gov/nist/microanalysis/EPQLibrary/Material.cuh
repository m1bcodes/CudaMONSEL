#ifndef _MATERIAL_CUH_
#define _MATERIAL_CUH_

#include "Composition.cuh"

namespace Material
{
   typedef ::Composition::data_type data_type;

   class Material : public Composition::Composition
   {
   public:
      __host__ __device__ Material(data_type density);
      __host__ __device__ Material(const Material& comp);
      __host__ __device__ Material(const Composition& comp, data_type density);
      Material(const Element::Element* elm[], const data_type density[]);
      __host__ __device__ Material(const Element::Element* elms[], int elmsLen, const data_type massFracs[], int massFracsLen, data_type density, const char* name);

      __host__ __device__ Material& operator=(const Material&);
      bool operator==(const Material&) const;

      void setDensity(data_type den);
      __host__ __device__ data_type getDensity() const;

      //template<typename T>
      //void defineByWeightFraction(LinkedListKV::Node<T, data_type>* map, data_type den)
      //{
      //   mDensity = den;
      //   Composition::defineByWeightFraction(map);
      //}

      void clear();
      data_type atomsPerCubicMeter(const Element::Element& elm) const;
      void defineByMaterialFraction(const Material* mats[], int matsLen, data_type matFracs[], int matFracsLen);
      void defineByMaterialFraction(const Composition* compositions[], int compositionsLen, data_type matFracs[], int matFracsLen, data_type density);
      //char* descriptiveString(bool normalize);
      int compareTo(const Composition& obj) const;
      int compareTo(const Material& obj) const;
      void replicate(const Material& mat);
      Material clone() const;
      unsigned int hashCode() const;
      bool equals(const Material& obj) const;
      bool almostEquals(const Material& other, data_type tol) const;

      __host__ __device__ virtual bool isSEmaterial() const;

   private:
      data_type mDensity;
   };

   extern const Material Default;
   extern __device__ const Material *d_Default;
   extern __global__ void initCuda();
}
#endif
