#include "gov\nist\nanoscalemetrology\JMONSEL\GanachaudMokraniPhononInelasticSM.cuh"
#include "gov\nist\microanalysis\EPQLibrary\PhysicalConstants.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

namespace GanachaudMokraniPhononInelasticSM
{
   GanachaudMokraniPhononInelasticSM::GanachaudMokraniPhononInelasticSM(double ratemultiplier, double phononE, double temperature, double eps0, double epsInfinity) : 
      ratemultiplier(ratemultiplier),
      phononE(phononE),
      temperature(temperature),
      eps0(eps0),
      epsInfinity(epsInfinity),
      occupationFactor(0.5 * (1. + (1. / (::exp(phononE / (PhysicalConstants::BoltzmannConstant * temperature)) - 1.)))),
      epsRatio((eps0 - epsInfinity) / eps0 / epsInfinity),
      prefactor((ratemultiplier * occupationFactor * epsRatio) / PhysicalConstants::BohrRadius)
   {
      if (ratemultiplier <= 0.) printf("Nonpositive ratemultiplier in GanachaudMokraniPhononInelasticSM constructor.");
      if (phononE <= 0.) printf("Nonpositive phononE in GanachaudMokraniPhononInelasticSM constructor.");
      if (epsRatio <= 0.) printf("(eps0-epsInfinity)/eps0/epsInfinity < 0 in GanachaudMokraniPhononInelasticSM constructor.");
   }

   ElectronT* GanachaudMokraniPhononInelasticSM::scatter(ElectronT* pe)
   {
      const double kE0 = pe->getEnergy();
      if (kE0 < phononE)
         return NULL;

      const double x = phononE / kE0; // Energy ratio

      const double randoms[] = { Math2::random(), Math2::random() };

      double costheta; // scattering angle
      if (x < 0.1)
         costheta = 1 + (((x * x) - (::pow(16., randoms[0]) * ::pow(x, 2. - (2. * randoms[0])))) / 8.);
      else { // Using general formula
         const double root = ::sqrt(1. - x);
         const double temp = ::pow(((-2. * (1 + root)) + x) / ((-2. * (1 - root)) + x), randoms[0]);
         costheta = temp + (((x - 2.) * (temp - 1)) / 2. / root);
      }
      const double phi = 2. * Math2::PI * randoms[1];

      pe->updateDirection(::acos(costheta), phi);
      pe->setEnergy(kE0 - phononE);
      pe->setPreviousEnergy(kE0);

      return NULL;
   }

   double GanachaudMokraniPhononInelasticSM::scatterRate(const ElectronT* pe)
   {
      const double kE = pe->getEnergy();
      if (kE < phononE)
         return 0.;
      const double x = phononE / kE; // Energy ratio
      /*
      * In the usual case the PE has energy of a few eV. Phonons typically have
      * energies ~0.1 eV or lower, so the above ratio is usually ~1/50 or so.
      * For such small x we calculate using a series expansion form. This is
      * simpler than the general expression, it is faster, and it avoids
      * numerical round-off problems. For larger x, which we expect to
      * encounter only rarely, we use the exact expression.
      */
      double result;
      if (x < 0.1) {
         result = prefactor * x * ::logf(4. / x);
         return result;
      }
      else if (x >= 1.) // phonon energy >= PE energy: no scattering possible
         return 0.;
      else {
         const double temp = ::sqrt(1. - x);
         result = prefactor * x * ::logf((1. + temp) / (1. - temp));
         return result;
      }
   }

   void GanachaudMokraniPhononInelasticSM::setMaterial(const MaterialT* mat)
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