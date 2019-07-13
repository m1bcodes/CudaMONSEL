#include "SEmaterial.cuh"

#include "gov\nist\microanalysis\EPQLIbrary\Composition.cuh"
#include "gov\nist\microanalysis\EPQLIbrary\AtomicSHell.cuh"

static const long serialVersionUID = 0x42;

namespace SEmaterial
{
   // Material properties appropriate to vacuum
   __host__ __device__ SEmaterial::SEmaterial() :
      Material(0),
      workfunction(0.),
      energyCBbottom(0.),
      eplasmon(0)
   {
   }

   __host__ __device__ SEmaterial::SEmaterial(const SEmaterial& other) :
      Material(other),
      workfunction(other.workfunction),
      energyCBbottom(other.energyCBbottom),
      eplasmon(other.eplasmon)
   {
   }

   __host__ __device__ SEmaterial::SEmaterial(const Composition& comp, double density) :
      Material(comp, density),
      workfunction(0.),
      energyCBbottom(0.),
      eplasmon(0)
   {
   }

   __host__ __device__ SEmaterial::SEmaterial(const Element::Element* elms[], int elemLen, const double weightFracs[], int wfLen, double density, char* name) :
      Material(elms, elemLen, weightFracs, wfLen, density, name),
      workfunction(0.),
      energyCBbottom(0.),
      eplasmon(0)
   {
   }

   __host__ __device__ SEmaterial::SEmaterial(const Material& mat) :
      Material(mat),
      workfunction(0.),
      energyCBbottom(0.),
      eplasmon(0)
   {
   }

   SEmaterial& SEmaterial::operator=(const SEmaterial& other)
   {
      if (this == &other) return *this;
      Material::operator=(other);
      workfunction = other.workfunction;
      energyCBbottom = other.energyCBbottom;
      eplasmon = other.eplasmon;

      return *this;
   }

   bool SEmaterial::operator==(const SEmaterial& other)
   {
      return this == &other;
   }

   //void SEmaterial::addBindingEnergy(double bindingEnergy, double density)
   //{
   //   if (bindingEnergy < 0.0) {
   //      printf("Binding energies must be positive.");
   //      return;
   //   }

   //   if (density < 0.0) {
   //      printf("Electron density must be positive.");
   //      return;
   //   }
   //   this->bindingEnergy.push_back(bindingEnergy);
   //   if (-bindingEnergy > energyCBbottom) {
   //      this->kineticEnergy.push_back(-bindingEnergy - energyCBbottom);
   //   }
   //   else {
   //      this->kineticEnergy.push_back(bindingEnergy + energyCBbottom);
   //   }
   //   electronDensity.push_back(density);
   //   version = (version == _UI32_MAX) ? 0L : version + 1L;
   //}

   //void SEmaterial::addBindingEnergy(double bindingEnergy, double kineticEnergy, double density)
   //{
   //   if (bindingEnergy < 0.0) {
   //      printf("Binding energies must be positive.");
   //      return;
   //   }
   //   if (kineticEnergy < 0.0) {
   //      printf("Kinetic energies must be positive.");
   //      return;
   //   }
   //   if (density < 0.0) {
   //      printf("Electron density must be positive.");
   //      return;
   //   }
   //   this->bindingEnergy.push_back(bindingEnergy);
   //   this->kineticEnergy.push_back(kineticEnergy);
   //   userSetKE = true;
   //   electronDensity.push_back(density);
   //   version = (version == _UI32_MAX) ? 0L : version + 1L;
   //}

   //void SEmaterial::addBindingEnergy(const VectorXd& bindingEnergy, const VectorXd& density)
   //{
   //   if (bindingEnergy.size() != density.size()) {
   //      printf("Unequal # of binding energies and densities");
   //   }
   //   for (double b : bindingEnergy) {
   //      if (b < 0.0) {
   //         printf("Binding energies must be positive.");
   //         return;
   //      }
   //   }
   //   for (double d : density) {
   //      if (d < 0.0) {
   //         printf("Electron density must be positive.");
   //         return;
   //      }
   //   }
   //   this->bindingEnergy.insert(this->bindingEnergy.end(), bindingEnergy.begin(), bindingEnergy.end());
   //   // Use default kinetic energy
   //   for (double b : bindingEnergy) {
   //      if (-b > energyCBbottom) {
   //         this->kineticEnergy.push_back(-b - energyCBbottom);
   //      }
   //      else {
   //         this->kineticEnergy.push_back(b + energyCBbottom);
   //      }
   //   }
   //   electronDensity.insert(electronDensity.end(), density.begin(), density.end());
   //   version = (version == _UI32_MAX) ? 0L : version + 1L;
   //}

   //void SEmaterial::addBindingEnergy(const VectorXd& bindingEnergy, const VectorXd& kineticEnergy, const VectorXd& density)
   //{
   //   // Error checking
   //   if ((bindingEnergy.size() != density.size()) || (kineticEnergy.size() != density.size())) {
   //      printf("Lists of energies and densities must be equal length");
   //      return;
   //   }
   //   for (double b : bindingEnergy) {
   //      if (b < 0.0) {
   //         printf("Binding energies must be positive.");
   //         return;
   //      }
   //   }
   //   for (double b : kineticEnergy) {
   //      if (b < 0.0) {
   //         printf("Kinetic energies must be positive.");
   //         return;
   //      }
   //   }
   //   for (double d : density) {
   //      if (d < 0.0) {
   //         printf("Electron density must be positive.");
   //         return;
   //      }
   //   }
   //   this->bindingEnergy.insert(this->bindingEnergy.end(), bindingEnergy.begin(), bindingEnergy.end());
   //   this->kineticEnergy.insert(this->kineticEnergy.end(), kineticEnergy.begin(), kineticEnergy.end());
   //   userSetKE = true;
   //   this->electronDensity.insert(this->electronDensity.end(), density.begin(), density.end());
   //   version = (version == _UI32_MAX) ? 0L : version + 1L;
   //}

   void SEmaterial::addCoreEnergy(double coreEnergy)
   {
      if (coreEnergy < 0.0) {
         printf("Core energies must be positive.");
         return;
      }
      this->coreEnergy.insert(coreEnergy);

      version = (version == _UI32_MAX) ? 0L : version + 1L;
   }

   void SEmaterial::addCoreEnergy(const Setd& coreEnergy)
   {
      // Error checking
      for (double cE : coreEnergy) {
         if (cE < 0.0) {
            printf("Core energies must be positive.");
         }
      }
      this->coreEnergy.insert(coreEnergy.begin(), coreEnergy.end());
      version = (version == _UI32_MAX) ? 0L : version + 1L;
   }

   //SEmaterial clone()
   //{
   //   SEmaterial res;
   //   res.replicate(*this);
   //   return res;
   //}

   long SEmaterial::get_version() const
   {
      return version;
   }

   //VectorXd SEmaterial::getBindingEnergyArray() const
   //{
   //   return bindingEnergy;
   //}

   __host__ __device__ Setd SEmaterial::getCoreEnergyArray() const
   {
      return coreEnergy;
   }

   __host__ __device__ double SEmaterial::getEFermi() const
   {
      return -energyCBbottom - workfunction;
   }

   //VectorXd SEmaterial::getElectronDensityArray() const
   //{
   //   return electronDensity;
   //}

   __host__ __device__ double SEmaterial::getEnergyCBbottom() const
   {
      return energyCBbottom;
   }

   double SEmaterial::getEplasmon() const
   {
      return eplasmon;
   }

   //VectorXd SEmaterial::getKineticEnergyArray() const
   //{
   //   return kineticEnergy;
   //}

   __host__ __device__ double SEmaterial::getWorkfunction() const
   {
      return workfunction;
   }

   __host__ __device__ double SEmaterial::getBandgap() const
   {
      return bandgap;
   }

   double SEmaterial::getEpsr() const
   {
      return epsr;
   }

   double SEmaterial::getDielectricBreakdownField() const
   {
      return dielectricBreakdownField;
   }

   //void SEmaterial::removeBindingEnergy(int index)
   //{
   //   bindingEnergy.erase(bindingEnergy.begin() + index);
   //   kineticEnergy.erase(kineticEnergy.begin() + index);
   //   electronDensity.erase(kineticEnergy.begin() + index);
   //   version = (version == _UI32_MAX) ? 0L : version + 1L;
   //}

   //void SEmaterial::removeCoreEnergy(double index)
   //{
   //   coreEnergy.erase(index);
   //   version = (version == _UI32_MAX) ? 0L : version + 1L;
   //}

   //void SEmaterial::replicate(const SEmaterial& mat)
   //{
   //   Material::replicate(mat);
   //   workfunction = mat.getWorkfunction();
   //   energyCBbottom = mat.getEnergyCBbottom();
   //   eplasmon = mat.getEplasmon();
   //   bindingEnergy.insert(bindingEnergy.end(), mat.bindingEnergy.begin(), mat.bindingEnergy.end());
   //   electronDensity.insert(electronDensity.end(), mat.electronDensity.begin(), mat.electronDensity.end());
   //   kineticEnergy.insert(kineticEnergy.end(), mat.kineticEnergy.begin(), mat.kineticEnergy.end());
   //   userSetKE = mat.userSetKE;
   //   coreEnergy.insert(mat.coreEnergy.begin(), mat.coreEnergy.end());
   //   bandgap = mat.getBandgap();
   //   epsr = mat.getEpsr();
   //   dielectricBreakdownField = mat.getDielectricBreakdownField();
   //}

   //void SEmaterial::setBindingEnergy(int index, double energy)
   //{
   //   if (energy < 0.) {
   //      printf("Binding energies must be positive.");
   //      return;
   //   }
   //   bindingEnergy[index] = energy;
   //   if (!userSetKE) {
   //      if (-energy > energyCBbottom) {
   //         this->kineticEnergy[index] = -energy - energyCBbottom;
   //      }
   //      else {
   //         this->kineticEnergy[index] = energy + energyCBbottom;
   //      }
   //   }
   //   version = (version == _UI32_MAX) ? 0L : version + 1L;
   //}

   __host__ __device__ void SEmaterial::setCoreEnergy()
   {
      this->coreEnergy.clear();
   }

   __host__ __device__ void SEmaterial::setCoreEnergy(const double coreEnergy[], int len)
   {
      setCoreEnergy();
      for (int i = 0; i < len; ++i) {
         if (coreEnergy[i] < 0.0) {
            printf("Core energies must be positive %.10e\n", coreEnergy[i]);
         }
         this->coreEnergy.insert(coreEnergy[i]);
      }

      //this->coreEnergy.insert(coreEnergy, coreEnergy + len);
      //if (*(this->coreEnergy.begin()) < 0) printf("Core energies must be positive.");

      version = (version == _UI32_MAX) ? 0L : version + 1L;
   }

   //void SEmaterial::setElectronDensity(int index, double density)
   //{
   //   if (density < 0.) {
   //      printf("Electron density must be positive.");
   //   }
   //   electronDensity[index] = density;
   //   version = (version == _UI32_MAX) ? 0L : version + 1L;
   //}

   __host__ __device__ void SEmaterial::setEnergyCBbottom(double energyCBbottom)
   {
      this->energyCBbottom = energyCBbottom;
      if (!userSetKE) {
         setKEtoDefault();
      }
      version = (version == _UI32_MAX) ? 0L : version + 1L;
   }

   __host__ __device__ void SEmaterial::setBandgap(double bandgap)
   {
      this->bandgap = bandgap;
      version = (version == _UI32_MAX) ? 0L : version + 1L;
   }

   void SEmaterial::setEplasmon(double eplasmon)
   {
      this->eplasmon = eplasmon;
      version = (version == _UI32_MAX) ? 0L : version + 1L;
   }

   static const double coreEnergyCutoff = ToSI::eV(20.);

   void SEmaterial::setEstimatedCoreEnergy()
   {
      setEstimatedCoreEnergy(coreEnergyCutoff);
   }

   void SEmaterial::setEstimatedCoreEnergy(double cutoff)
   {
      setCoreEnergy(); // Clear any existing ones.
      auto constituentElements = this->getElementSet();
      for (auto el : constituentElements) {
         int i = 0;
         while (true) {
            AtomicShellT as(*el, i);
            double shellenergy = as.getGroundStateOccupancy() > 0 ? as.getEdgeEnergy() : NAN;
            if (shellenergy <= cutoff) break;
            addCoreEnergy(shellenergy);
            i++;
         }
      }
   }

   __host__ __device__ void SEmaterial::setKEtoDefault()
   {
      //for (int i = 0; i < kineticEnergy.size(); i++) {
      //   if (-bindingEnergy[i] > energyCBbottom) {
      //      kineticEnergy[i] = -bindingEnergy[i] - energyCBbottom;
      //   }
      //   else {
      //      kineticEnergy[i] = bindingEnergy[i] + energyCBbottom;
      //   }
      //}
      userSetKE = false;
      version = (version == _UI32_MAX) ? 0L : version + 1L;
   }

   //void SEmaterial::setKineticEnergy(int index, double energy)
   //{
   //   if (energy < 0.) {
   //      printf("Kinetic energies must be positive.");
   //      return;
   //   }
   //   kineticEnergy[index] = energy;
   //   version = (version == _UI32_MAX) ? 0L : version + 1L;
   //}

   __host__ __device__ void SEmaterial::setWorkfunction(double workfunction)
   {
      this->workfunction = workfunction;
      if (energyCBbottom > -workfunction) {
         energyCBbottom = -workfunction;
      }
      version = (version == _UI32_MAX) ? 0L : version + 1L;
   }

   void SEmaterial::setEpsr(double epsr)
   {
      this->epsr = epsr;
      version = (version == _UI32_MAX) ? 0L : version + 1L;
   }

   void SEmaterial::setDielectricBreakdownField(double breakdownField)
   {
      dielectricBreakdownField = breakdownField;
      version = (version == _UI32_MAX) ? 0L : version + 1L;
   }

   __host__ __device__ bool SEmaterial::isSEmaterial() const
   {
      return true;
   }
}
