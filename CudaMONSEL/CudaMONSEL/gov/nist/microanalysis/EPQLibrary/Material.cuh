#ifndef _MATERIAL_CUH_
#define _MATERIAL_CUH_

#include "Composition.cuh"

namespace Material
{
   class Material : public Composition::Composition
   {
   protected:

   public:
      Material(double density);
      Material(const Material& comp);
      Material(const Composition& comp, double density);
      Material(const Element::Element* elm[], double density[]);
      Material(const Element::Element* elms[], int elmsLen, double massFracs[], int massFracsLen, double density, const char* name);

      Material& operator=(const Material&);
      bool operator==(const Material&) const;

      void setDensity(double den);
      double getDensity() const;

      //template<typename T>
      //void defineByWeightFraction(LinkedListKV::Node<T, double>* map, double den)
      //{
      //   mDensity = den;
      //   Composition::defineByWeightFraction(map);
      //}

      void clear();
      double atomsPerCubicMeter(const Element::Element& elm) const;
      void defineByMaterialFraction(const Material* mats[], int matsLen, double matFracs[], int matFracsLen);
      void defineByMaterialFraction(const Composition* compositions[], int compositionsLen, double matFracs[], int matFracsLen, double density);
      //char* descriptiveString(bool normalize);
      int compareTo(const Composition& obj) const;
      int compareTo(const Material& obj) const;
      void replicate(const Material& mat);
      Material clone() const;
      unsigned int hashCode() const;
      bool equals(const Material& obj) const;
      bool almostEquals(const Material& other, double tol) const;

      virtual bool isSEmaterial() const;

   private:
      double mDensity;
   };

   static Material Default(0);
}
#endif
