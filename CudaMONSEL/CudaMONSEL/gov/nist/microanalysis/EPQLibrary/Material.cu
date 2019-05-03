#include "Material.cuh"
#include "ToSI.cuh"
#include "Element.cuh"

#include <math.h>

namespace Material
{
   const long serialVersionUID = 0x42;
   const double DEFAULT_DENSITY = 5.0 * 1.0e-3 / (1.0e-2 * 1.0e-2 * 1.0e-2); //ToSI::gPerCC(5.0);

   Material::Material(double density)
   {
      mDensity = density;
   }

   Material::Material(const Material& other)
   {
      Composition::replicate(other);
      mDensity = other.mDensity;
   }

   Material::Material(const Element::Element* elms[], int elmsLen, double massFracs[], int massFracsLen, double density, char* name) : Composition(elms, elmsLen, massFracs, massFracsLen, name)
   {
      mDensity = density;
      renormalize();
   }

   Material::Material(const Composition& comp, double density)
   {
      Composition::replicate(comp);
      mDensity = density;
   }

   Material::Material(const Element::Element* elm[], double density[]) : Composition(elm, 1, density, 1, elm[0]->toString())
   {
      mDensity = density[0];
      renormalize();
   }

   Material& Material::operator=(const Material& other)
   {
      Composition::replicate(other);
      mDensity = other.mDensity;

      return *this;
   }

   void Material::setDensity(double den)
   {
      if (mDensity != den) {
         mDensity = den;
      }
   }

   double Material::getDensity() const
   {
      return mDensity;
   }

   void Material::clear()
   {
      Composition::clear();
      mDensity = DEFAULT_DENSITY;
   }

   double Material::atomsPerCubicMeter(const Element::Element& elm) const
   {
      if (getElementCount() == 0 || mDensity > 0.0) {
         printf("Material::atomsPerCubicMeter: invalid composition");
         return 0.0;
      }
      return weightFraction(elm, true) * mDensity / elm.getMass();
   }

   void Material::defineByMaterialFraction(const Material* mats[], int matsLen, double matFracs[], int matFracsLen)
   {
      std::vector<const Composition*> comp;
      for (int k = 0; k < matsLen; ++k) {
         comp.push_back(mats[k]);
      }
      Composition::defineByMaterialFraction(comp.data(), matsLen, matFracs, matFracsLen);
      double den = 0.0;
      for (int i = 0; i < matsLen; ++i) {
         den += matFracs[i] * mats[i]->getDensity();
      }
      setDensity(den);
   }

   void Material::defineByMaterialFraction(const Composition* compositions[], int compositionsLen, double matFracs[], int matFracsLen, double density)
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

   bool Material::almostEquals(const Material& other, double tol) const
   {
      Material otherMat = (Material)other;
      return Composition::almostEquals(other, tol) && (abs(getDensity() - otherMat.getDensity()) / fmax(getDensity(), otherMat.getDensity()) < tol);
   }
}
