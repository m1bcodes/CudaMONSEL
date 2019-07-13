#ifndef _SE_MATERIAL_CUH_
#define _SE_MATERIAL_CUH_

#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

#include "Amphibian\unordered_set.cuh"

namespace SEmaterial
{
   typedef amp::unordered_set<double, Hasher::DoubleHashFcn, Comparator::DoubleCompareFcn> Setd;

   class SEmaterial : public Material::Material
   {
   protected:
      //void replicate(const SEmaterial& mat);

   public:
      __host__ __device__ SEmaterial();
      __host__ __device__ SEmaterial(const SEmaterial&);
      __host__ __device__ SEmaterial(const Composition& comp, double density);
      __host__ __device__ SEmaterial(const Element::Element* elms[], int elemLen, const double weightFracs[], int wfLen, double density, char* name);
      __host__ __device__ SEmaterial(const Material& mat);

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
      //VectorXd getBindingEnergyArray() const;
      __host__ __device__ Setd getCoreEnergyArray() const;
      __host__ __device__ double getEFermi() const;
      //VectorXd getElectronDensityArray() const;
      __host__ __device__ double getEnergyCBbottom() const;
      double getEplasmon() const;
      //VectorXd getKineticEnergyArray() const;
      __host__ __device__ double getWorkfunction() const;
      __host__ __device__ double getBandgap() const;
      double getEpsr() const;
      double getDielectricBreakdownField() const;

      //void removeBindingEnergy(int index);
      //void removeCoreEnergy(double index);

      //void setBindingEnergy(int index, double energy);
      __host__ __device__ void setCoreEnergy();
      __host__ __device__ void setCoreEnergy(const double coreEnergy[], int len);
      //void setElectronDensity(int index, double density);
      __host__ __device__ void setEnergyCBbottom(double energyCBbottom);
      __host__ __device__ void setBandgap(double bandgap);
      void setEplasmon(double eplasmon);
      void setEstimatedCoreEnergy();
      void setEstimatedCoreEnergy(double cutoff);
      __host__ __device__ void setKEtoDefault();
      //void setKineticEnergy(int index, double energy);
      __host__ __device__ void setWorkfunction(double workfunction);
      void setEpsr(double epsr);
      void setDielectricBreakdownField(double breakdownField);

      __host__ __device__ bool isSEmaterial() const override;

   private:
      double workfunction; // work function
      double dielectricBreakdownField = INFINITY;
      double epsr = 1.; // relative dielectric function
      double eplasmon; // plasmon resonance energy
      double energyCBbottom; // energy of conduction band bottom
      double bandgap = 0.; // width of the bandgap
      long version = 0L; // Updates each time the data change

      //VectorXd bindingEnergy;
      //VectorXd electronDensity;
      //VectorXd kineticEnergy;
      bool userSetKE = false; // Flag = true when user sets a kinetic

      Setd coreEnergy;
   };
}

#endif