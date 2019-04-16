//#ifndef _SE_MATERIAL_CUH_
//#define _SE_MATERIAL_CUH_
//
//#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
//
//namespace SEmaterial
//{
//   class SEmaterial : Material::Material
//   {
//   protected:
//      void replicate(const SEmaterial& mat);
//
//   public:
//      typedef std::vector<double> SEmaterialArrayT;
//      typedef std::set<double> SEmaterialSetT;
//
//   public:
//      SEmaterial();
//      SEmaterial(const SEmaterial&);
//      SEmaterial(const Composition& comp, double density);
//      SEmaterial(Element::Element elms[], int elemLen, double weightFracs[], int wfLen, double density, char* name);
//      SEmaterial(const Material& mat);
//
//      SEmaterial& operator=(const SEmaterial&);
//      bool operator==(const SEmaterial&);
//
//      void addBindingEnergy(double bindingEnergy, double density);
//      void addBindingEnergy(double bindingEnergy, double kineticEnergy, double density);
//      void addBindingEnergy(const SEmaterialArrayT& bindingEnergy, const SEmaterialArrayT& density);
//      void addBindingEnergy(const SEmaterialArrayT& bindingEnergy, const SEmaterialArrayT& kineticEnergy, const SEmaterialArrayT& density);
//
//      void addCoreEnergy(double coreEnergy);
//      void addCoreEnergy(const SEmaterialSetT& coreEnergy);
//
//      //SEmaterial clone();
//      long get_version() const;
//      SEmaterialArrayT getBindingEnergyArray() const;
//      SEmaterialSetT getCoreEnergyArray() const;
//      double getEFermi() const;
//      SEmaterialArrayT getElectronDensityArray() const;
//      double getEnergyCBbottom() const;
//      double getEplasmon() const;
//      SEmaterialArrayT getKineticEnergyArray() const;
//      double getWorkfunction() const;
//      double getBandgap() const;
//      double getEpsr() const;
//      double getDielectricBreakdownField() const;
//
//      void removeBindingEnergy(int index);
//      void removeCoreEnergy(double index);
//
//      void setBindingEnergy(int index, double energy);
//      void setCoreEnergy();
//      void setCoreEnergy(const SEmaterialSetT& coreEnergy);
//      void setElectronDensity(int index, double density);
//      void setEnergyCBbottom(double energyCBbottom);
//      void setBandgap(double bandgap);
//      void setEplasmon(double eplasmon);
//      //void setEstimatedCoreEnergy();
//      //void setEstimatedCoreEnergy(double cutoff);
//      void setKEtoDefault();
//      void setKineticEnergy(int index, double energy);
//      void setWorkfunction(double workfunction);
//      void setEpsr(double epsr);
//      void setDielectricBreakdownField(double breakdownField);
//
//   private:
//      double workfunction; // work function
//      double dielectricBreakdownField = INFINITY;
//      double epsr = 1.; // relative dielectric function
//      double eplasmon; // plasmon resonance energy
//      double energyCBbottom; // energy of conduction band bottom
//      double bandgap = 0.; // width of the bandgap
//      long version = 0L; // Updates each time the data change
//
//      SEmaterialArrayT bindingEnergy;
//      SEmaterialArrayT electronDensity;
//      SEmaterialArrayT kineticEnergy;
//      bool userSetKE = false; // Flag = true when user sets a kinetic
//
//      SEmaterialSetT coreEnergy;
//   };
//
//}
//
//#endif