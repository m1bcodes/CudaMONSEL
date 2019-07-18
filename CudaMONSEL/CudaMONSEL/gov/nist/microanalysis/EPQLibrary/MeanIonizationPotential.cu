#include "gov\nist\microanalysis\EPQLibrary\MeanIonizationPotential.cuh"
#include "gov\nist\microanalysis\EPQLibrary\CaveatBase.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Composition.cuh"
#include "gov\nist\microanalysis\EPQLibrary\FromSI.cuh"

#include "CudaUtil.h"

namespace MeanIonizationPotential
{
   __host__ __device__ MeanIonizationPotential::MeanIonizationPotential(StringT name, const ReferenceT &reference) : AlgorithmClass("Mean Ionization Potential", name, reference)
   {
   }

   __host__ __device__ void MeanIonizationPotential::initializeDefaultStrategy()
   {
   }

   StringT caveat(const ElementT &el)
   {
      return CaveatBase::None;
   }

   StringT caveat(const CompositionT &comp)
   {
      StringT res(CaveatBase::None);
      for (auto el : comp.getElementSet())
         res = CaveatBase::append(res, caveat(*el));
      return res;
   }

   double MeanIonizationPotential::computeLn(const CompositionT &comp) const
   {
      double m = 0.0;
      double lnJ = 0.0;
      for (auto &el : comp.getElementSet()) {
         double cz_a = comp.weightFraction(*el, true) * el->getAtomicNumber() / el->getAtomicWeight();
         m += cz_a;
         lnJ += cz_a * ::log(FromSI::keV(compute(*el)));
      }
      return ToSI::keV(::exp(lnJ / m));
   }

   Reference::CrudeReference SternheimerCR("Sternheimer quoted in Berger MJ, Seltzer S. NASA Technical Publication SP-4012 (1964)");
   __host__ __device__ Sternheimer64MeanIonizationPotential::Sternheimer64MeanIonizationPotential() :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      MeanIonizationPotential("Sternheimer 1964", *Reference::d_NullReference)
#else
      MeanIonizationPotential("Sternheimer 1964", SternheimerCR)
#endif
   {
   }

   __host__ __device__ double Sternheimer64MeanIonizationPotential::compute(const ElementT &el) const
   {
      const double z = el.getAtomicNumber();
      return ToSI::eV(9.76 * z + 58.8 * ::powf(z, -0.19));
   }

   const Sternheimer64MeanIonizationPotential Sternheimer64Ref;
   const MeanIonizationPotential &Sternheimer64 = Sternheimer64Ref;
   __device__ const MeanIonizationPotential *dSternheimer64;

   __host__ __device__ double computeSternheimer64(const ElementT &el)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      if (dSternheimer64 == nullptr) {
         printf("computeSternheimer64: not initialized on device");
         return NAN;
      }
      return dSternheimer64->compute(el);
#else
      return Sternheimer64.compute(el);
#endif
   }

   Reference::CrudeReference BergerSeltzerCR("Berger and Seltzer as implemented by CITZAF 3.06");
   __host__ __device__ BergerAndSeltzerCITZAFMeanIonizationPotential::BergerAndSeltzerCITZAFMeanIonizationPotential() :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      MeanIonizationPotential("Berger & Seltzer as per JTA", *Reference::d_NullReference)
#else
      MeanIonizationPotential("Berger & Seltzer as per JTA", BergerSeltzerCR)
#endif
   {
   }

   __host__ __device__ double BergerAndSeltzerCITZAFMeanIonizationPotential::compute(const ElementT &el) const
   {
      const double z = el.getAtomicNumber();
      return ToSI::eV(9.76 * z + 58.5 * ::powf(z, -0.19));
   }

   const BergerAndSeltzerCITZAFMeanIonizationPotential BergerAndSeltzerCITZAFRef;
   const MeanIonizationPotential &BergerAndSeltzerCITZAF = BergerAndSeltzerCITZAFRef;
   __device__ const MeanIonizationPotential *dBergerAndSeltzerCITZAF;

   Reference::CrudeReference Bloch33CR("Bloch F, F. Z. Phys. 81, 363 (1933)");
   __host__ __device__ Bloch33MeanIonizationPotential::Bloch33MeanIonizationPotential() :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      MeanIonizationPotential("Bloch 1933", *Reference::d_NullReference)
#else
      MeanIonizationPotential("Bloch 1933", Bloch33CR)
#endif
   {
   }

   __host__ __device__ double Bloch33MeanIonizationPotential::compute(const ElementT &el) const
   {
      const double z = el.getAtomicNumber();
      return ToSI::eV(13.5 * z);
   }

   const Bloch33MeanIonizationPotential Bloch33Ref;
   const MeanIonizationPotential &Bloch33 = Bloch33Ref;
   __device__ const MeanIonizationPotential *dBloch33;

   Reference::CrudeReference Wilson41CR("Wilson RR. Phys Rev. 60. 749 (1941)");
   __host__ __device__ Wilson41MeanIonizationPotential::Wilson41MeanIonizationPotential() :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      MeanIonizationPotential("Wilson 1941", *Reference::d_NullReference)
#else
      MeanIonizationPotential("Wilson 1941", Wilson41CR)
#endif
   {
   }

   __host__ __device__ double Wilson41MeanIonizationPotential::compute(const ElementT &el) const
   {
      const double z = el.getAtomicNumber();
      return ToSI::eV(11.5 * z);
   }

   const Wilson41MeanIonizationPotential Wilson41Ref;
   const MeanIonizationPotential &Wilson41 = Wilson41Ref;
   __device__ const MeanIonizationPotential *dWilson41;

   __host__ __device__ double computeWilson41(const ElementT &el)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      if (dWilson41 == nullptr) {
         printf("computeWilson41: not initialized on device");
         return NAN;
      }
      return dWilson41->compute(el);
#else
      return Wilson41.compute(el);
#endif
   }

   Reference::CrudeReference Springer67CR("Springer G. Meues Jahrbuch Fuer Mineralogie, Monatshefte (1967) 9/10, p. 304");
   __host__ __device__ Springer67MeanIonizationPotential::Springer67MeanIonizationPotential() :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      MeanIonizationPotential("Springer 1967", *Reference::d_NullReference)
#else
      MeanIonizationPotential("Springer 1967", Springer67CR)
#endif
   {
   }

   __host__ __device__ double Springer67MeanIonizationPotential::compute(const ElementT &el) const
   {
      const double z = el.getAtomicNumber();
      return ToSI::eV(z * (9.0 * (1.0 + ::powf(z, -0.67)) + 0.03 * z));
   }

   const Springer67MeanIonizationPotential Springer67Ref;
   const MeanIonizationPotential &Springer67 = Springer67Ref;
   __device__ const MeanIonizationPotential *dSpringer67;

   Reference::CrudeReference Heinrich70CR("Heinrich KFJ, Yakowitz H. Mikrochim Acta (1970) p 123");
   __host__ __device__ Heinrich70MeanIonizationPotential::Heinrich70MeanIonizationPotential() :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      MeanIonizationPotential("Heinrich & Yakowitz 1970", *Reference::d_NullReference)
#else
      MeanIonizationPotential("Heinrich & Yakowitz 1970", Heinrich70CR)
#endif
   {
   }

   __host__ __device__ double Heinrich70MeanIonizationPotential::compute(const ElementT &el) const
   {
      const double z = el.getAtomicNumber();
      return ToSI::eV(z * (12.4 + 0.027 * z));
   }

   const Heinrich70MeanIonizationPotential Heinrich70Ref;
   const MeanIonizationPotential &Heinrich70 = Heinrich70Ref;
   __device__ const MeanIonizationPotential *dHeinrich70;

   Reference::CrudeReference Duncumb69CR("Duncumb P, Shields-Mason PK, DeCasa C. Proc. 5th Int. Congr. on X-ray Optics and Microanalysis, Springer, Berlin, 1969 p. 146");
   __host__ __device__ Duncumb69MeanIonizationPotential::Duncumb69MeanIonizationPotential() :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      MeanIonizationPotential("Duncumb & DeCasa 1969", *Reference::d_NullReference)
#else
      MeanIonizationPotential("Duncumb & DeCasa 1969", Duncumb69CR)
#endif
   {
   }

   __host__ __device__ double Duncumb69MeanIonizationPotential::compute(const ElementT &el) const
   {
      const double z = el.getAtomicNumber();
      return ToSI::eV((14.0 * (1.0 - ::expf(-0.1 * z)) + 75.5 / ::powf(z, z / 7.5) - z / (100 + z)) * z);
   }

   const Duncumb69MeanIonizationPotential Duncumb69Ref;
   const MeanIonizationPotential &Duncumb69 = Duncumb69Ref;
   __device__ const MeanIonizationPotential *dDuncumb69;

   Reference::CrudeReference Zeller75CR("Zeller C in Ruste J, Gantois M, J. Phys. D. Appl. Phys 8, 872 (1975)");
   __host__ __device__ Zeller75MeanIonizationPotential::Zeller75MeanIonizationPotential() :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      MeanIonizationPotential("Zeller 1975", *Reference::d_NullReference)
#else
      MeanIonizationPotential("Zeller 1975", Zeller75CR)
#endif
   {
   }

   __host__ __device__ double Zeller75MeanIonizationPotential::compute(const ElementT &el) const
   {
      const double z = el.getAtomicNumber();
      return ToSI::eV((10.04 + 8.25 * ::expf(-z / 11.22)) * z);
   }

   const Zeller75MeanIonizationPotential Zeller75Ref;
   const MeanIonizationPotential &Zeller75 = Zeller75Ref;
   __device__ const MeanIonizationPotential *dZeller75;

   // https://www.oreilly.com/library/view/c-cookbook/0596007612/ch03s06.html
   static double sciToDub(const std::string& str)
   {
      std::stringstream ss(str);
      double d = 0;
      ss >> d;

      if (ss.fail()) {
         std::string s = "Unable to format ";
         s += str;
         s += " as a number!";
         throw (s);
      }

      return d;
   }

   Reference::CrudeReference Berger64CR("Berger MJ, Seltzer S. NASA Technical Publication SP-4012 (1964)");
   __host__ __device__ Berger64MeanIonizationPotential::Berger64MeanIonizationPotential() :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      MeanIonizationPotential("Berger & Seltzer 1964", *Reference::d_NullReference)
#else
      MeanIonizationPotential("Berger & Seltzer 1964", Berger64CR)
#endif
   {
   }

   void Berger64MeanIonizationPotential::readTabulatedValues()
   {
      if (!mMeasured.empty()) return;

      std::string name(".\\gov\\nist\\microanalysis\\EPQLibrary\\BergerSeltzer64.csv");
      printf("Reading: %s\n", name.c_str());
      try {
         std::ifstream file(name);
         if (!file.good()) throw 0;
         mMeasured.reserve(92);
         for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
            if ((*loop)[0][0] != '/')
               mMeasured.push_back(ToSI::eV(sciToDub((*loop)[0])));
         }
         file.close();
      }
      catch (std::exception&) {
         printf("Fatal error while attempting to load the mean ionization potential data file.");
      }
   }

   __host__ __device__ const VectorXd &Berger64MeanIonizationPotential::getData() const
   {
      return mMeasured;
   }

   //__device__ void Berger64MeanIonizationPotential::copyData(const double *data, const unsigned int len)
   //{
   //   mMeasured.resize(len);
   //   memcpy(mMeasured.data(), data, len * sizeof(double));
   //}

   __host__ __device__ double Berger64MeanIonizationPotential::compute(const ElementT &el) const
   {
      return mMeasured[el.getAtomicNumber() - 1];
   }

   Berger64MeanIonizationPotential Berger64Ref;
   Berger64MeanIonizationPotential &Berger64 = Berger64Ref;
   __device__ Berger64MeanIonizationPotential *dBerger64;

   Reference::CrudeReference Berger83CR("Berger MJ, Seltzer S. NBSIR 82-2550-A - US Dept of Commerce, Washington DC (1983)");
   __host__ __device__ Berger83MeanIonizationPotential::Berger83MeanIonizationPotential() :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      MeanIonizationPotential("Berger & Seltzer 1983", *Reference::d_NullReference)
#else
      MeanIonizationPotential("Berger & Seltzer 1983", Berger83CR)
#endif
   {
   }

   void Berger83MeanIonizationPotential::readTabulatedValues()
   {
      if (!mMeasured.empty()) return;

      std::string name(".\\gov\\nist\\microanalysis\\EPQLibrary\\BergerSeltzer83.csv");
      printf("Reading: %s\n", name.c_str());
      try {
         std::ifstream file(name);
         if (!file.good()) throw 0;
         mMeasured.reserve(100);
         for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
            if ((*loop)[0][0] != '/')
               mMeasured.push_back(ToSI::eV(sciToDub((*loop)[0])));
         }
         file.close();
      }
      catch (std::exception&) {
         printf("Fatal error while attempting to load the mean ionization potential data file.");
      }
   }

   __host__ __device__ const VectorXd &Berger83MeanIonizationPotential::getData() const
   {
      return mMeasured;
   }

   //__device__ void Berger83MeanIonizationPotential::copyData(const double *data, const unsigned int len)
   //{
   //   mMeasured.resize(len);
   //   memcpy(mMeasured.data(), data, len * sizeof(double));
   //}

   __host__ __device__ double Berger83MeanIonizationPotential::compute(const ElementT &el) const
   {
      return mMeasured[el.getAtomicNumber() - 1];
   }

   Berger83MeanIonizationPotential Berger83Ref;
   Berger83MeanIonizationPotential &Berger83 = Berger83Ref;
   __device__ Berger83MeanIonizationPotential *dBerger83;

   const AlgorithmClassT * mAllImplementations[] = {
      &Berger64,
      &Berger83,
      &Bloch33,
      &Duncumb69,
      &BergerAndSeltzerCITZAF,
      &Heinrich70,
      &Springer67,
      &Sternheimer64,
      &Wilson41,
      &Zeller75
   };

   //__device__ const AlgorithmClassT * dAllImplementations[] = {
   //   dBerger64,
   //   dBerger83,
   //   dBloch33,
   //   dDuncumb69,
   //   dBergerAndSeltzerCITZAF,
   //   dHeinrich70,
   //   dSpringer67,
   //   dSternheimer64,
   //   dWilson41,
   //   dZeller75
   //};

   __global__ void initCuda()
   {
      dBerger64 = new Berger64MeanIonizationPotential();
      dBerger83 = new Berger83MeanIonizationPotential();
      dBloch33 = new Bloch33MeanIonizationPotential();
      dDuncumb69 = new Duncumb69MeanIonizationPotential();
      dBergerAndSeltzerCITZAF = new BergerAndSeltzerCITZAFMeanIonizationPotential();
      dHeinrich70 = new Heinrich70MeanIonizationPotential();
      dSpringer67 = new Springer67MeanIonizationPotential();
      dSternheimer64 = new Sternheimer64MeanIonizationPotential();
      dWilson41 = new Wilson41MeanIonizationPotential();
      dZeller75 = new Zeller75MeanIonizationPotential();
   }

   __global__ void copyDataToBerger64(double *data, int len)
   {
      dBerger64->copyData<double>(data, len);
   }

   __global__ void copyDataToBerger83(double *data, int len)
   {
      dBerger83->copyData<double>(data, len);
   }

   void copyDataToCuda()
   {
      double *Berger64data = nullptr;
      checkCudaErrors(cudaMalloc((void **)&Berger64data, Berger64.getData().size() * sizeof(double)));
      checkCudaErrors(cudaMemcpy(Berger64data, Berger64.getData().data(), Berger64.getData().size() * sizeof(double), cudaMemcpyHostToDevice));
      copyDataToBerger64 << <1, 1 >> >(Berger64data, Berger64.getData().size());
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      checkCudaErrors(cudaFree(Berger64data));

      double *Berger83data = nullptr;
      checkCudaErrors(cudaMalloc((void **)&Berger83data, Berger83.getData().size() * sizeof(double)));
      checkCudaErrors(cudaMemcpy(Berger83data, Berger83.getData().data(), Berger83.getData().size() * sizeof(double), cudaMemcpyHostToDevice));
      copyDataToBerger83 << <1, 1 >> >(Berger83data, Berger83.getData().size());
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      checkCudaErrors(cudaFree(Berger83data));
   }

   AlgorithmClassT const * const * MeanIonizationPotential::getAllImplementations() const
   {
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      return dAllImplementations;
//#else
//      return mAllImplementations;
//#endif
      return mAllImplementations;
   }
}
