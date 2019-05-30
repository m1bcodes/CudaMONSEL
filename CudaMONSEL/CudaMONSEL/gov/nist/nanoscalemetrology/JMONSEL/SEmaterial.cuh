#ifndef _SE_MATERIAL_CUH_
#define _SE_MATERIAL_CUH_

#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace SEmaterial
{
   class SEmaterial : public Material::Material
   {
   protected:
      void replicate(const SEmaterial& mat);

   public:
      SEmaterial();
      SEmaterial(const SEmaterial&);
      SEmaterial(const Composition& comp, double density);
      SEmaterial(const Element::Element* elms[], int elemLen, const double weightFracs[], int wfLen, double density, char* name);
      SEmaterial(const Material& mat);

      SEmaterial& operator=(const SEmaterial&);
      bool operator==(const SEmaterial&);

      //void addBindingEnergy(double bindingEnergy, double density);
      //void addBindingEnergy(double bindingEnergy, double kineticEnergy, double density);
      //void addBindingEnergy(const VectorXd& bindingEnergy, const VectorXd& density);
      //void addBindingEnergy(const VectorXd& bindingEnergy, const VectorXd& kineticEnergy, const VectorXd& density);

      void addCoreEnergy(double coreEnergy);
      void addCoreEnergy(const Setd& coreEnergy);

      //SEmaterial clone();
      long get_version() const;
      VectorXd getBindingEnergyArray() const;
      Setd getCoreEnergyArray() const;
      double getEFermi() const;
      VectorXd getElectronDensityArray() const;
      double getEnergyCBbottom() const;
      double getEplasmon() const;
      VectorXd getKineticEnergyArray() const;
      double getWorkfunction() const;
      double getBandgap() const;
      double getEpsr() const;
      double getDielectricBreakdownField() const;

      void removeBindingEnergy(int index);
      void removeCoreEnergy(double index);

      void setBindingEnergy(int index, double energy);
      void setCoreEnergy();
      void setCoreEnergy(const double coreEnergy[], int len);
      void setElectronDensity(int index, double density);
      void setEnergyCBbottom(double energyCBbottom);
      void setBandgap(double bandgap);
      void setEplasmon(double eplasmon);
      void setEstimatedCoreEnergy();
      void setEstimatedCoreEnergy(double cutoff);
      void setKEtoDefault();
      void setKineticEnergy(int index, double energy);
      void setWorkfunction(double workfunction);
      void setEpsr(double epsr);
      void setDielectricBreakdownField(double breakdownField);

      bool isSEmaterial() const override;

   private:
      double workfunction; // work function
      double dielectricBreakdownField = INFINITY;
      double epsr = 1.; // relative dielectric function
      double eplasmon; // plasmon resonance energy
      double energyCBbottom; // energy of conduction band bottom
      double bandgap = 0.; // width of the bandgap
      long version = 0L; // Updates each time the data change

      VectorXd bindingEnergy;
      VectorXd electronDensity;
      VectorXd kineticEnergy;
      bool userSetKE = false; // Flag = true when user sets a kinetic

      Setd coreEnergy;
   };
}

#endif