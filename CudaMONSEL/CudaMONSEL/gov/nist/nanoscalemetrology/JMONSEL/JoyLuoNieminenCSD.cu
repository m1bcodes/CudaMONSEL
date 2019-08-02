#include "gov\nist\nanoscalemetrology\JMONSEL\JoyLuoNieminenCSD.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\EPQLibrary\MeanIonizationPotential.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\SEmaterial.cuh"

namespace JoyLuoNieminenCSD
{
   __host__ __device__ JoyLuoNieminenCSD::JoyLuoNieminenCSD(SEmaterialT& mat, float breakE) : mat(mat), breakE(breakE)
   {
      if (!(breakE > (mat.getWorkfunction() + ToSI::eV(1.f))))
         printf("JoyLuoNieminenCSD::JoyLuoNieminenCSD: Supplied breakpoint energy is too small.");
      setMaterial(&mat);
   }

   JoyLuoNieminenCSD::JoyLuoNieminenCSD(MaterialT& mat, float bh, float breakE) : mat(mat)
   {
      if (breakE > (bh + ToSI::eV(1.f)))
         this->breakE = breakE;
      else
         printf("JoyLuoNieminenCSD::JoyLuoNieminenCSD: Supplied breakpoint energy is too small.");
      setMaterial(mat, bh);
   }

   void JoyLuoNieminenCSD::setMaterial(const MaterialT& mat, float bh)
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
      bhplus1eV = bh + ToSI::eV(1.f);
      if (breakE < bhplus1eV)
         breakE = bhplus1eV;

      unsigned int i = 0;
      for (const auto &el : mat.getElementSet()) {
         recipJ[i] = 1.166f / (el->getAtomicNumber() < 13 ? MeanIonizationPotential::computeWilson41(*el) : MeanIonizationPotential::computeSternheimer64(*el));
         beta[i] = 1.f - (recipJ[i] * bhplus1eV);
         coef[i] = 2.01507e-28f * mat.getDensity() * mat.weightFraction(*el, true) * el->getAtomicNumber() / el->getAtomicWeight();
         ++i;
      }
      minEforTracking = bh;

      // Determine the proportionality constant
      gamma = 0.f;
      for (int i = 0; i < nce; i++)
         gamma += coef[i] * ::logf((recipJ[i] * breakE) + beta[i]);
      //gamma /= ::pow(breakE, 3.5); // note: breaks if powf is used instead of pow (gamma becomes #1.ind due to division by 0)
      gamma /= breakE;
      gamma /= breakE;
      gamma /= breakE;
      gamma /= ::sqrtf(breakE);
      if (gamma != gamma) printf("JoyLuoNieminenCSD::init(): gamma is NaN");
   }

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ static float MaxLossFraction = 0.1;
#else
   static float MaxLossFraction = 0.1;
#endif

   __host__ __device__ SlowingDownAlg::data_type JoyLuoNieminenCSD::compute(::SlowingDownAlg::data_type len, const ElectronT *pe) const
   {
      return compute(len, pe->getEnergy());
   }

   __host__ __device__ SlowingDownAlg::data_type JoyLuoNieminenCSD::compute(const ::SlowingDownAlg::data_type len, const ::SlowingDownAlg::data_type kE) const
   {
      if ((nce == 0) || (kE < minEforTracking) || (kE <= 0.f))
         return 0.f;
      if (kE <= breakE)
         return (kE / ::powf(1.f + (1.5f * gamma * len * kE * ::sqrtf(kE)), 2.f / 3.f)) - kE;

      // Otherwise, Joy/Luo formula
      float loss = 0.f;
      for (int i = 0; i < nce; i++)
         loss += coef[i] * ::logf((recipJ[i] * kE) + beta[i]);
      loss *= len / kE;

      if (loss <= (MaxLossFraction * kE))
         return -loss;
      loss = compute(len / 2.f, kE); // loss from 1st half of step
      if (loss != loss) printf("JoyLuoNieminenCSD::compute(): loss is NaN");
      return loss + compute(len / 2.f, kE + loss); // add loss from second half
   }

   float JoyLuoNieminenCSD::getBreakE() const
   {
      return breakE;
   }

   void JoyLuoNieminenCSD::setBreakE(float breakE)
   {
      this->breakE = breakE;
      setMaterial(mat, bh);
   }

   __host__ __device__ StringT JoyLuoNieminenCSD::toString() const
   {
      return "JoyLuoNieminenCSD(" + StringT(mat.toString()) + "," + amp::to_string(bh) + "," + amp::to_string(breakE) + ")";
   }
}
