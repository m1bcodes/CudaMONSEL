#include "Material.cuh"
#include "ToSI.cuh"
#include "Element.cuh"

#include <math.h>

namespace Material
{
   const long serialVersionUID = 0x42;
   const data_type DEFAULT_DENSITY = 5.0f * 1.0e-3f / (1.0e-2f * 1.0e-2f * 1.0e-2f); //ToSI::gPerCC(5.0);

   __host__ __device__ Material::Material(data_type density) :
      mDensity(density)
   {
   }

   __host__ __device__ Material::Material(const Material& other) :
      Composition(other),
      mDensity(other.mDensity)
   {
   }

   __host__ __device__ Material::Material(const Element::Element* elms[], int elmsLen, const data_type massFracs[], int massFracsLen, data_type density, const char* name) :
      Composition(elms, elmsLen, massFracs, massFracsLen, name),
      mDensity(density)
   {
   }

   __host__ __device__ Material::Material(const Composition& comp, data_type density) :
      Composition(comp),
      mDensity(density)
   {
   }

   Material::Material(const Element::Element* elm[], const data_type density[]) :
      Composition(elm, 1, density, 1, elm[0]->toString()),
      mDensity(density[0])
   {
   }

   __host__ __device__ Material& Material::operator=(const Material& other)
   {
      Composition::replicate(other);
      mDensity = other.mDensity;

      return *this;
   }

   void Material::setDensity(data_type den)
   {
      if (mDensity != den) {
         mDensity = den;
      }
   }

   __host__ __device__ data_type Material::getDensity() const
   {
      return mDensity;
   }

   void Material::clear()
   {
      Composition::clear();
      mDensity = DEFAULT_DENSITY;
   }

   data_type Material::atomsPerCubicMeter(const Element::Element& elm) const
   {
      if (getElementCount() == 0 || mDensity > 0.0) {
         printf("Material::atomsPerCubicMeter: invalid composition");
         return 0.0;
      }
      return weightFraction(elm, true) * mDensity / elm.getMass();
   }

   void Material::defineByMaterialFraction(const Material* mats[], int matsLen, data_type matFracs[], int matFracsLen)
   {
      //std::vector<const Composition*> comp;
      const Composition** comp = new const Composition*[matsLen];
      for (int k = 0; k < matsLen; ++k) {
         //comp.push_back(mats[k]);
         comp[k] = mats[k];
      }
      Composition::defineByMaterialFraction(comp, matsLen, matFracs, matFracsLen);
      delete[] comp;
      data_type den = 0.0;
      for (int i = 0; i < matsLen; ++i) {
         den += matFracs[i] * mats[i]->getDensity();
      }
      setDensity(den);
   }

   void Material::defineByMaterialFraction(const Composition* compositions[], int compositionsLen, data_type matFracs[], int matFracsLen, data_type density)
   {
      Composition::defineByMaterialFraction(compositions, compositionsLen, matFracs, matFracsLen);
      setDensity(density);
   }

   //char* Material::descriptiveString(bool normalize)
   //{
   //   if (getElementCount() > 0) {
   //      final NumberFormat nf = NumberFormat.getInstance();
   //      nf.setMaximumFractionDigits(1);
   //      final StringBuffer sb = new StringBuffer(super.descriptiveString(normalize));
   //      final int p = sb.lastIndexOf("]");
   //      sb.insert(p, "," + nf.format(FromSI.gPerCC(mDensity)) + " g/cc");
   //      return sb.toString();
   //   }
   //   else
   //      return "None";
   //}

   int Material::compareTo(const Composition& obj) const
   {
      return Composition::compareTo(obj);
   }

   int Material::compareTo(const Material& obj) const
   {
      int res = Composition::compareTo(obj);
      if (res == 0) {
         res = !(mDensity == obj.mDensity);
      }
      return res;
   }

   void Material::replicate(const Material& mat)
   {
      Composition::replicate(mat);
      mDensity = mat.mDensity;
   }

   Material Material::clone() const
   {
      Material res(mDensity);
      res.replicate(*this);
      return res;
   }

   unsigned int Material::hashCode() const
   {
      int PRIME = 31;
      int result = Composition::hashCode();
      long long temp;
      memcpy(&temp, &mDensity, sizeof(mDensity));
      result = PRIME * result + (int)(temp ^ (temp >> 32));
      if (result == INT_MAX)
         result = INT_MIN;
      return result;
   }

   bool Material::operator==(const Material& other) const
   {
      return Composition::equals(other) && mDensity == other.mDensity;
   }

   bool Material::equals(const Material& other) const
   {
      return this == &other || *this == other;
   }

   bool Material::almostEquals(const Material& other, data_type tol) const
   {
      Material otherMat = (Material)other;
      return Composition::almostEquals(other, tol) && (abs(getDensity() - otherMat.getDensity()) / fmax(getDensity(), otherMat.getDensity()) < tol);
   }

   __host__ __device__ bool Material::isSEmaterial() const
   {
      return false;
   }

   __device__ const Material* d_Default;
   const Material Default(0);

   __global__ void initCuda()
   {
      d_Default = new Material(0);
   }
}
