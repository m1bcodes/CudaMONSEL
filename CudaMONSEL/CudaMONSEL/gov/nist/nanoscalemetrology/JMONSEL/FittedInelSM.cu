#include "gov\nist\nanoscalemetrology\JMONSEL\FittedInelSM.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\SEmaterial.cuh"

#include "Amphibian\random.cuh"

namespace FittedInelSM
{
   __host__ __device__ FittedInelSM::FittedInelSM(const SEmaterialT& mat, float energySEgen, const SlowingDownAlgT& sdAlg) : sdAlg(sdAlg), energySEgen(energySEgen)
   {
      setMaterial(&mat);
   }

   __host__ __device__ ElectronT* FittedInelSM::scatter(ElectronT* pe)
   {
      const double phi = 2.f * Math2::PI * Random::random();
      const double theta = ::acosf(1.f - (2.f * Random::random()));

      ElectronT* newElectron = new ElectronT(*pe, theta, phi, energySEgen + eFermi);
      if (!newElectron) printf("FittedInelSM::FittedInelSM: failed creating electron.\n");
      return newElectron;
   }

   __host__ __device__::ScatterMechanism::data_type FittedInelSM::scatterRate(const ElectronT* pe)
   {
      if (pe->getEnergy() <= (energySEgen + eFermi))
         return 0.f;
      return (-sdAlg.compute(1.e-10f, pe) * 1.e10f) / energySEgen;
   }

   __host__ __device__ void FittedInelSM::setMaterial(const MaterialT* mat)
   {
      if (!mat->isSEmaterial()) {
         printf("FittedInelSM::setMaterial: not SEmaterial\n");
      }
      eFermi = ((SEmaterialT*)mat)->getEFermi();
   }

   StringT FittedInelSM::toString() const
   {
      return "FittedInelSM(" + amp::to_string(eFermi) + "," + amp::to_string(energySEgen) + "," + sdAlg.toString() + ")";
   }
}
