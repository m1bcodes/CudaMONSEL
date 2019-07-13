#include "gov\nist\nanoscalemetrology\JMONSEL\JoyLuoNieminenCSD.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\EPQLibrary\MeanIonizationPotential.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\SEmaterial.cuh"

namespace JoyLuoNieminenCSD
{
   __host__ __device__ JoyLuoNieminenCSD::JoyLuoNieminenCSD(SEmaterialT& mat, double breakE) : mat(mat)
   {
      if (breakE > (mat.getWorkfunction() + ToSI::eV(1.)))
         this->breakE = breakE;
      else
         printf("JoyLuoNieminenCSD::JoyLuoNieminenCSD: Supplied breakpoint energy is too small.");
      setMaterial(&mat);
   }

   JoyLuoNieminenCSD::JoyLuoNieminenCSD(MaterialT& mat, double bh, double breakE) : mat(mat)
   {
      if (breakE > (bh + ToSI::eV(1.)))
         this->breakE = breakE;
      else
         printf("JoyLuoNieminenCSD::JoyLuoNieminenCSD: Supplied breakpoint energy is too small.");
      setMaterial(mat, bh);
   }

   void JoyLuoNieminenCSD::setMaterial(const MaterialT& mat, double bh)
   {
      this->mat = mat;
      this->bh = bh;
      init();
   }

   __host__ __device__ void JoyLuoNieminenCSD::setMaterial(const SEmaterialT* mat)
   {
      this->mat = *mat;
      bh = -mat->getEnergyCBbottom();
      init();
   }

   __host__ __device__ void JoyLuoNieminenCSD::init()
   {
      nce = mat.getElementCount();
      if (nce == 0)
         return;
      recipJ.resize(nce);
      coef.resize(nce);
      beta.resize(nce);
      bhplus1eV = bh + ToSI::eV(1.);
      if (breakE < bhplus1eV)
         breakE = bhplus1eV;

      int i = 0;
      for (const auto &el : mat.getElementSet()) {
         recipJ[i] = 1.166 / (el->getAtomicNumber() < 13 ? MeanIonizationPotential::computeWilson41(*el) : MeanIonizationPotential::computeSternheimer64(*el));
         beta[i] = 1. - (recipJ[i] * bhplus1eV);
         coef[i] = 2.01507E-28 * mat.getDensity() * mat.weightFraction(*el, true) * el->getAtomicNumber() / el->getAtomicWeight();
         ++i;
      }
      minEforTracking = bh;

      // Determine the proportionality constant
      gamma = 0.;
      for (int i = 0; i < nce; i++)
         gamma += coef[i] * ::log((recipJ[i] * breakE) + beta[i]);
      gamma /= ::pow(breakE, 3.5); // note: breaks if powf is used instead of pow (gamma becomes #1.ind due to division by 0)
#ifdef _DEBUG
      if (gamma != gamma) printf("JoyLuoNieminenCSD::init(): gamma is NaN");
#endif
   }

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ static double MaxLossFraction = 0.1;
#else
   static double MaxLossFraction = 0.1;
#endif

   __host__ __device__ double JoyLuoNieminenCSD::compute(double len, const ElectronT *pe) const
   {
      return compute(len, pe->getEnergy());
   }

   __host__ __device__ double JoyLuoNieminenCSD::compute(const double len, const double kE) const
   {
      if ((nce == 0) || (kE < minEforTracking) || (kE <= 0.))
         return 0.;
      if (kE <= breakE)
         return (kE / ::pow(1. + (1.5 * gamma * len * kE * ::sqrt(kE)), 2. / 3.)) - kE;

      // Otherwise, Joy/Luo formula
      double loss = 0.;
      for (int i = 0; i < nce; i++)
         loss += coef[i] * ::log((recipJ[i] * kE) + beta[i]);
      loss *= len / kE;

      if (loss <= (MaxLossFraction * kE))
         return -loss;
      loss = compute(len / 2., kE); // loss from 1st half of step
#ifdef _DEBUG
      if (loss != loss) printf("JoyLuoNieminenCSD::compute(): loss is NaN");
#endif
      return loss + compute(len / 2., kE + loss); // add loss from second half
   }

   double JoyLuoNieminenCSD::getBreakE() const
   {
      return breakE;
   }

   void JoyLuoNieminenCSD::setBreakE(double breakE)
   {
      this->breakE = breakE;
      setMaterial(mat, bh);
   }

   __host__ __device__ StringT JoyLuoNieminenCSD::toString() const
   {
      return "JoyLuoNieminenCSD(" + StringT(mat.toString()) + "," + amp::to_string(bh) + "," + amp::to_string(breakE) + ")";
   }
}
