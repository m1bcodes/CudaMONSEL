#ifndef _SE_MATERIAL_CUH_
#define _SE_MATERIAL_CUH_

#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

#include "Amphibian\unordered_set.cuh"

namespace SEmaterial
{
   typedef ::Material::data_type data_type;
   typedef amp::unordered_set<data_type, Hasher::DoubleHashFcn, Comparator::DoubleCompareFcn> Setd;

   class SEmaterial : public Material::Material
   {
   protected:
      //void replicate(const SEmaterial& mat);

   public:
      __host__ __device__ SEmaterial();
      __host__ __device__ SEmaterial(const SEmaterial&);
      __host__ __device__ SEmaterial(const Composition& comp, data_type density);
      __host__ __device__ SEmaterial(const Element::Element* elms[], int elemLen, const data_type weightFracs[], int wfLen, data_type density, char* name);
      __host__ __device__ SEmaterial(const Material& mat);

      SEmaterial& operator=(const SEmaterial&);
      bool operator==(const SEmaterial&);

      //void addBindingEnergy(data_type bindingEnergy, data_type density);
      //void addBindingEnergy(data_type bindingEnergy, data_type kineticEnergy, data_type density);
      //void addBindingEnergy(const VectorXd& bindingEnergy, const VectorXd& density);
      //void addBindingEnergy(const VectorXd& bindingEnergy, const VectorXd& kineticEnergy, const VectorXd& density);

      void addCoreEnergy(data_type coreEnergy);
      void addCoreEnergy(const Setd& coreEnergy);

      //SEmaterial clone();
      long get_version() const;
      //VectorXd getBindingEnergyArray() const;
      __host__ __device__ Setd getCoreEnergyArray() const;
      __host__ __device__ data_type getEFermi() const;
      //VectorXd getElectronDensityArray() const;
      __host__ __device__ data_type getEnergyCBbottom() const;
      data_type getEplasmon() const;
      //VectorXd getKineticEnergyArray() const;
      __host__ __device__ data_type getWorkfunction() const;
      __host__ __device__ data_type getBandgap() const;
      data_type getEpsr() const;
      data_type getDielectricBreakdownField() const;

      //void removeBindingEnergy(int index);
      //void removeCoreEnergy(data_type index);

      //void setBindingEnergy(int index, data_type energy);
      __host__ __device__ void setCoreEnergy();
      __host__ __device__ void setCoreEnergy(const data_type coreEnergy[], int len);
      //void setElectronDensity(int index, data_type density);
      __host__ __device__ void setEnergyCBbottom(data_type energyCBbottom);
      __host__ __device__ void setBandgap(data_type bandgap);
      void setEplasmon(data_type eplasmon);
      void setEstimatedCoreEnergy();
      void setEstimatedCoreEnergy(data_type cutoff);
      __host__ __device__ void setKEtoDefault();
      //void setKineticEnergy(int index, data_type energy);
      __host__ __device__ void setWorkfunction(data_type workfunction);
      void setEpsr(data_type epsr);
      void setDielectricBreakdownField(data_type breakdownField);

      __host__ __device__ bool isSEmaterial() const override;

   private:
      data_type workfunction; // work function
      data_type dielectricBreakdownField = INFINITY;
      data_type epsr = 1.f; // relative dielectric function
      data_type eplasmon; // plasmon resonance energy
      data_type energyCBbottom; // energy of conduction band bottom
      data_type bandgap = 0.; // width of the bandgap
      long version = 0L; // Updates each time the data change

      //VectorXd bindingEnergy;
      //VectorXd electronDensity;
      //VectorXd kineticEnergy;
      bool userSetKE = false; // Flag = true when user sets a kinetic

      Setd coreEnergy;
   };
}

#endif