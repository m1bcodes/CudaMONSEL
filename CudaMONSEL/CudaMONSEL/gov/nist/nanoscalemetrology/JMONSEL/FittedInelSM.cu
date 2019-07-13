#include "gov\nist\nanoscalemetrology\JMONSEL\FittedInelSM.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\SEmaterial.cuh"

#include "Amphibian\random.cuh"

namespace FittedInelSM
{
   __host__ __device__ FittedInelSM::FittedInelSM(const SEmaterialT& mat, double energySEgen, const SlowingDownAlgT& sdAlg) : sdAlg(sdAlg), energySEgen(energySEgen)
   {
      setMaterial(&mat);
   }

   __host__ __device__ ElectronT* FittedInelSM::scatter(ElectronT* pe)
   {
      const double phi = 2 * Math2::PI * Random::random();
      const double theta = ::acos(1. - (2. * Random::random()));
      return new ElectronT(*pe, theta, phi, energySEgen + eFermi);
   }

   __host__ __device__ double FittedInelSM::scatterRate(const ElectronT* pe)
   {
      if (pe->getEnergy() <= (energySEgen + eFermi))
         return 0.;
      return (-sdAlg.compute(1.e-10, pe) * 1.e10) / energySEgen;
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
