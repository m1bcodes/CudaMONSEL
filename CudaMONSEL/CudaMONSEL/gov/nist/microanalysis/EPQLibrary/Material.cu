#include "Material.cuh"
#include "ToSI.cuh"
#include "Element.cuh"

#include <math.h>

namespace Material
{
   __device__ const long serialVersionUID = 0x42;
   __device__ double DEFAULT_DENSITY = 5.0 * 1.0e-3 / (1.0e-2 * 1.0e-2 * 1.0e-2); //ToSI::gPerCC(5.0);

   __device__ Material::Material(double density)
   {
      mDensity = density;
   }

   __device__ Material::Material(Element::Element elms[], int elmsLen, double massFracs[], int massFracsLen, double density, char* name) : Composition(elms, elmsLen, massFracs, massFracsLen, name)
   {
      mDensity = density;
      renormalize();
   }

   __device__ Material::Material(Composition comp, double density)
   {
      mDensity = density;
      Composition::replicate(comp);
   }

   __device__ Material::Material(Element::Element elm[], double density[]) : Composition(elm, 1, density, 1, (char *)elm[0].toString())
   {
      mDensity = density[0];
      renormalize();
   }

   __device__ void Material::setDensity(double den)
   {
      if (mDensity != den) {
         mDensity = den;
      }
   }

   __device__ double Material::getDensity()
   {
      return mDensity;
   }

   __device__ void Material::clear()
   {
      Composition::clear();
      mDensity = DEFAULT_DENSITY;
   }

   __device__ double Material::atomsPerCubicMeter(Element::Element elm)
   {
      if ((getElementCount() == 0) || (mDensity > 0.0)) {
         printf("Material::atomsPerCubicMeter: invalid composition");
         return 0.0;
      }
      return weightFraction(elm, true) * mDensity / elm.getMass();
   }

   __device__ void Material::defineByMaterialFraction(Material mats[], int matsLen, double matFracs[], int matFracsLen)
   {
      Composition::defineByMaterialFraction(mats, matsLen, matFracs, matFracsLen);
      double den = 0.0;
      for (int i = 0; i < matsLen; ++i) {
         den += matFracs[i] * mats[i].getDensity();
      }
      setDensity(den);
   }

   __device__ void Material::defineByMaterialFraction(Composition compositions[], int compositionsLen, double matFracs[], int matFracsLen, double density)
   {
      Composition::defineByMaterialFraction(compositions, compositionsLen, matFracs, matFracsLen);
      setDensity(density);
   }

   //__device__ char* Material::descriptiveString(bool normalize)
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

   __device__ int Material::compareTo(Composition obj)
   {
      return Composition::compareTo(obj);
   }

   __device__ int Material::compareTo(Material obj)
   {
      int res = Composition::compareTo(obj);
      if ((res == 0)) {
         res = !(mDensity == obj.mDensity);
      }
      return res;
   }

   __device__ void Material::replicate(Material mat)
   {
      Composition::replicate(mat);
      mDensity = mat.mDensity;
   }

   __device__ Material Material::clone()
   {
      Material res(mDensity);
      res.replicate(*this);
      return res;
   }

   //__device__ int Material::hashCode()
   //{
   //   if (mHashCode == String::MAX_SIGNED_INTEGER) {
   //      int PRIME = 31;
   //      int result = Composition::hashCode();
   //      long temp;
   //      temp = Double.doubleToLongBits(mDensity);
   //      result = PRIME * result + (int)(temp ^ (temp >> > 32));
   //      if (result == Integer.MAX_VALUE)
   //         result = Integer.MIN_VALUE;
   //      mHashCode = result;
   //   }
   //   return mHashCode;
   //}

   __device__ bool Material::equals(Material other)
   {
      if (!Composition::equals(other)) {
         return false;
      }
      if (mDensity != other.mDensity) {
         return false;
      }
      return true;
   }

   __device__ bool Material::almostEquals(Material other, double tol)
   {
      Material otherMat = (Material)other;
      return Composition::almostEquals(other, tol) && (abs(getDensity() - otherMat.getDensity()) / fmax(getDensity(), otherMat.getDensity()) < tol);
   }
}
