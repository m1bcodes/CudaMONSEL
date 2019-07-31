#include "gov\nist\nanoscalemetrology\JMONSEL\GanachaudMokraniPolaronTrapSM.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\SEmaterial.cuh"

namespace GanachaudMokraniPolaronTrapSM
{
   __host__ __device__ GanachaudMokraniPolaronTrapSM::GanachaudMokraniPolaronTrapSM(data_type prefactor, data_type extinction) :
      prefactor(prefactor),
      extinction(extinction),
      CUTOFF(-::logf(10.f / prefactor) / extinction)
   {
      if (prefactor < 0.f) printf("ERROR: Nonpositive prefactor in GanachaudMokraniPolaronTrapSM constructor.");
      if (extinction < 0.f) printf("ERROR: Nonpositive extinction in GanachaudMokraniPolaronTrapSM constructor.");
   }

   __host__ __device__ ElectronT* GanachaudMokraniPolaronTrapSM::scatter(ElectronT* pe)
   {
      pe->setEnergy(0.f); // So listeners, if any, will record the energy change.
      pe->setTrajectoryComplete(true); // It's trapped
      return NULL;
   }

   __host__ __device__ data_type GanachaudMokraniPolaronTrapSM::scatterRate(const ElectronT* pe)
   {
      const data_type kE0 = pe->getEnergy();
      if (kE0 > CUTOFF) return 0.f;
      const data_type result = prefactor * ::expf(-extinction * kE0);
      return result;
   }

   __host__ __device__ void GanachaudMokraniPolaronTrapSM::setMaterial(const MaterialT* mat)
   {
   }

   StringT GanachaudMokraniPolaronTrapSM::toString() const
   {
      return "FittedInelSM(" + amp::to_string(prefactor) + "," + amp::to_string(extinction) + "," + amp::to_string(CUTOFF) + ")";
   }
}
