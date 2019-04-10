#ifndef _MATERIAL_CUH_
#define _MATERIAL_CUH_

#include "Composition.cuh"

namespace Material
{
   class Material : public Composition::Composition
   {
   protected:
      Material(Element::Element elms[], int elmsLen, double massFracs[], int massFracsLen, double density, char* name);

   public:
      Material(double density);
      Material(const Material& comp);
      Material(const Composition& comp, double density);
      Material(Element::Element elm[], double density[]);

      Material& operator=(const Material&);
      bool operator==(const Material&) const;

      void setDensity(double den);
      double getDensity();

      //template<typename T>
      //void defineByWeightFraction(LinkedListKV::Node<T, double>* map, double den)
      //{
      //   mDensity = den;
      //   Composition::defineByWeightFraction(map);
      //}

      void clear();
      double atomsPerCubicMeter(Element::Element elm);
      void defineByMaterialFraction(Material mats[], int matsLen, double matFracs[], int matFracsLen);
      void defineByMaterialFraction(Composition compositions[], int compositionsLen, double matFracs[], int matFracsLen, double density);
      //char*  descriptiveString(bool normalize);
      int compareTo(Composition obj);
      int compareTo(Material obj);
      void replicate(Material mat);
      Material clone();
      unsigned int hashCode();
      bool equals(const Material& obj) const;
      bool almostEquals(Material other, double tol);

   private:
      double mDensity;
   };
}

#endif
