#include "gov\nist\nanoscalemetrology\JMONSEL\GanachaudMokraniPhononInelasticSM.cuh"
#include "gov\nist\microanalysis\EPQLibrary\PhysicalConstants.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

#include "Amphibian\random.cuh"

namespace GanachaudMokraniPhononInelasticSM
{
   __host__ __device__ GanachaudMokraniPhononInelasticSM::GanachaudMokraniPhononInelasticSM(data_type ratemultiplier, data_type phononE, data_type temperature, data_type eps0, data_type epsInfinity) : 
      ratemultiplier(ratemultiplier),
      phononE(phononE),
      temperature(temperature),
      eps0(eps0),
      epsInfinity(epsInfinity),
      occupationFactor(0.5f * (1.f + (1.f / (::expf(phononE / (PhysicalConstants::BoltzmannConstant * temperature)) - 1.f)))),
      epsRatio((eps0 - epsInfinity) / eps0 / epsInfinity),
      prefactor((ratemultiplier * occupationFactor * epsRatio) / PhysicalConstants::BohrRadius)
   {
      if (ratemultiplier <= 0.) printf("Nonpositive ratemultiplier in GanachaudMokraniPhononInelasticSM constructor.");
      if (phononE <= 0.) printf("Nonpositive phononE in GanachaudMokraniPhononInelasticSM constructor.");
      if (epsRatio <= 0.) printf("(eps0-epsInfinity)/eps0/epsInfinity < 0 in GanachaudMokraniPhononInelasticSM constructor.");
   }

   __host__ __device__ ElectronT* GanachaudMokraniPhononInelasticSM::scatter(ElectronT* pe)
   {
      const data_type kE0 = pe->getEnergy();
      if (kE0 < phononE)
         return nullptr;

      const data_type x = phononE / kE0; // Energy ratio

      const data_type randoms[] = { Random::random(), Random::random() };

      data_type costheta; // scattering angle
      if (x < 0.1f)
         costheta = 1.f + (((x * x) - (::powf(16.f, randoms[0]) * ::powf(x, 2.f - (2.f * randoms[0])))) / 8.f);
      else { // Using general formula
         const data_type root = ::sqrtf(1. - x);
         const data_type temp = ::powf(((-2. * (1 + root)) + x) / ((-2. * (1 - root)) + x), randoms[0]);
         costheta = temp + (((x - 2.) * (temp - 1)) / 2. / root);
      }
      const data_type phi = 2.f * Math2::PI * randoms[1];

      pe->updateDirection(::acosf(costheta), phi);
      pe->setEnergy(kE0 - phononE);
      pe->setPreviousEnergy(kE0);

      return nullptr;
   }

   __host__ __device__ data_type GanachaudMokraniPhononInelasticSM::scatterRate(const ElectronT* pe)
   {
      const data_type kE = pe->getEnergy();
      if (kE < phononE)
         return 0.f;
      const data_type x = phononE / kE; // Energy ratio
      /*
      * In the usual case the PE has energy of a few eV. Phonons typically have
      * energies ~0.1 eV or lower, so the above ratio is usually ~1/50 or so.
      * For such small x we calculate using a series expansion form. This is
      * simpler than the general expression, it is faster, and it avoids
      * numerical round-off problems. For larger x, which we expect to
      * encounter only rarely, we use the exact expression.
      */
      data_type result;
      if (x < 0.1f) {
         result = prefactor * x * ::logf(4.f / x);
         return result;
      }
      else if (x >= 1.f) // phonon energy >= PE energy: no scattering possible
         return 0.;
      else {
         const data_type temp = ::sqrt(1.f - x);
         result = prefactor * x * ::logf((1.f + temp) / (1.f - temp));
         return result;
      }
   }

   __host__ __device__ void GanachaudMokraniPhononInelasticSM::setMaterial(const MaterialT* mat)
   {
      /*
      * There's nothing to do here. This is a required method, but the phonon
      * model doesn't require any parameters not already passed in the
      * constructor.
      */
   }

   const char* GanachaudMokraniPhononInelasticSM::toString()
   {
      sprintf(buff, "GanachaudMokraniPhononInelasticSM(%.10e, %.10e, %.10e, %.10e, %.10e)", ratemultiplier, phononE, temperature, eps0, epsInfinity);
      return buff;
   }
}