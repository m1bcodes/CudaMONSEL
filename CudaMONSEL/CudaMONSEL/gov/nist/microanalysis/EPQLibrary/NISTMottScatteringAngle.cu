#include "gov\nist\microanalysis\EPQLibrary\NISTMottScatteringAngle.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Reference.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\PhysicalConstants.cuh"
#include "gov\nist\microanalysis\EPQLibrary\FromSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ScreenedRutherfordScatteringAngle.cuh"

#include "CudaUtil.h"

#include "Amphibian\random.cuh"

namespace NISTMottScatteringAngle
{
   static const Reference::Author* auRef[] = { &Reference::CPowell, &Reference::FSalvat, &Reference::AJablonski };
   static const Reference::WebSite REFERENCE("http://www.nist.gov/srd/nist64.htm", "NIST Electron Elastic-Scattering Cross-Section Database version 3.1", "2007 AUGUST 24", auRef, 3);


#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const int SPWEM_LEN = 61;
   __constant__ const int X1_LEN = 201;
   __constant__ const float DL50 = 3.91202300543f;
   __constant__ const float PARAM = 0.09985774245f;

   __constant__ static const float MAX_NISTMOTT = 3.2043531e-15;
#else
   const int SPWEM_LEN = 61;
   const int X1_LEN = 201;
   const float DL50 = ::log(50.0);
   const float PARAM = (::log(2.0e4) - DL50) / 60.0;

   static const float MAX_NISTMOTT = ToSI::keV(20.0);
#endif

   __host__ __device__ static double value(double a, double b, double c, double y0, double y1, double y2, double x)
   {
      return (x - b) * (x - c) * y0 / ((a - b) * (a - c)) + (x - a) * (x - c) * y1 / ((b - a) * (b - c)) + (x - a) * (x - b) * y2 / ((c - a) * (c - b));
   }

   // https://www.oreilly.com/library/view/c-cookbook/0596007612/ch03s06.html
   static double sciToDub(const std::string& str)
   {
      std::string tmp = str.substr(str.find_first_not_of(" "));
      std::stringstream ss(tmp);
      double d = 0;
      ss >> d;

      if (ss.fail()) {
         std::string s = "Unable to format ";
         s += tmp;
         s += " as a number!";
         throw (s);
      }

      return (d);
   }

   void NISTMottScatteringAngle::loadData(int an)
   {
      std::string name(an < 10 ? ".\\gov\\nist\\microanalysis\\EPQLibrary\\NistXSec/E0" + std::to_string(an) + ".D64" : ".\\gov\\nist\\microanalysis\\EPQLibrary\\NistXSec/E" + std::to_string(an) + ".D64");
      printf("Reading: %s\n", name.c_str());
      try {
         std::ifstream t(name);
         if (!t.good()) throw 0;
         std::string line;
         std::getline(t, line);
         for (int j = 0; j < SPWEM_LEN; ++j) {
            std::getline(t, line);
            mSpwem[j] = sciToDub(line);
            for (int i = 0; i < X1_LEN; ++i) {
               std::getline(t, line);
               mX1[j][i] = sciToDub(line);
            }
         }
         t.close();
      }
      catch (std::exception&) {
         printf("Unable to construct NISTMottScatteringAngle: %s\n", name.c_str());
      }
   }

   __host__ __device__ NISTMottScatteringAngle::NISTMottScatteringAngle(const ElementT& elm) :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      RandomizedScatterT("NIST Elastic cross-section", *Reference::d_NullReference),
#else
      RandomizedScatterT("NIST Elastic cross-section", REFERENCE),
#endif
      mElement(elm),
      mSpwem(SPWEM_LEN, 0),
      mX1(SPWEM_LEN, VectorXf(X1_LEN, 0)),
      mRutherford(ScreenedRutherfordScatteringAngle::getSRSA(elm.getAtomicNumber()))
   {
#if (!(defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0)))
      loadData(elm.getAtomicNumber());
#endif
   }

   StringT NISTMottScatteringAngle::toString() const
   {
      return "CrossSection[NIST-Mott," + StringT(mElement.toAbbrev()) + "]";
   }

   __host__ __device__ const ElementT& NISTMottScatteringAngle::getElement() const
   {
      return mElement;
   }
   
   __host__ __device__ double NISTMottScatteringAngle::totalCrossSection(const double energy) const
   {
      if (energy < MAX_NISTMOTT) {
         const double scale = PhysicalConstants::BohrRadius * PhysicalConstants::BohrRadius;
         const double logE = ::logf(FromSI::eV(energy));
         const int j = 1 + (int)((logE - DL50) / PARAM);
         if (j == 1)
            return value(DL50, DL50 + PARAM, DL50 + 2.0 * PARAM, mSpwem[0], mSpwem[1], mSpwem[2], logE) * scale;
         else if (j == SPWEM_LEN)
            return value(DL50 + 58.0 * PARAM, DL50 + 59.0 * PARAM, DL50 + 60.0 * PARAM, mSpwem[SPWEM_LEN - 3], mSpwem[SPWEM_LEN - 2], mSpwem[SPWEM_LEN - 1], logE) * scale;
         else {
            double e0 = DL50 + (j - 2) * PARAM;
            return value(e0, e0 + PARAM, e0 + 2.0 * PARAM, mSpwem[j - 2], mSpwem[j - 1], mSpwem[j], logE) * scale;
         }
      }
      else {
         return mRutherford.totalCrossSection(energy);
      }
   }

   __host__ __device__ double NISTMottScatteringAngle::randomScatteringAngle(const double energy) const
   {
      if (energy < MAX_NISTMOTT) {
         const double logE = ::log(FromSI::eV(energy));
         const int j = (int)((logE - DL50) / PARAM); // offset to zero-based
         const double e2 = DL50 + (j + 1) * PARAM;
         const double e1 = e2 - PARAM;
         const int i = (logE - e1 < e2 - logE ? j : j + 1); // offset to zero-based
         if (!((i >= 0) && (i < SPWEM_LEN))) printf("%d %s %lf %s %lf %s %lf\n", i, "\t", FromSI::eV(energy), "\t", e1, "\t", e2);
         // via j
         const int k = Random::randomInt(200); // offset to zero-based
         const double x = (mX1[i][k + 1] - mX1[i][k]) * Random::random();
         const double q = mX1[i][k] + x;
         const double com = 1.0 - 2.0 * q * q;
         return com > -1.0 ? (com < 1.0 ? ::acos(com) : 0.0) : PhysicalConstants::PI;
      }
      else {
         return mRutherford.randomScatteringAngle(energy);
      }
   }

   __host__ __device__ const VectorXf& NISTMottScatteringAngle::getSpwem() const
   {
      return mSpwem;
   }

   __host__ __device__ const MatrixXf& NISTMottScatteringAngle::getX1() const
   {
      return mX1;
   }

   const NISTMottScatteringAngle* mScatter[113];

   void init()
   {
      mScatter[1] = new NISTMottScatteringAngle(Element::H);
      mScatter[2] = new NISTMottScatteringAngle(Element::He);
      mScatter[3] = new NISTMottScatteringAngle(Element::Li);
      mScatter[4] = new NISTMottScatteringAngle(Element::Be);
      mScatter[5] = new NISTMottScatteringAngle(Element::B);
      mScatter[6] = new NISTMottScatteringAngle(Element::C);
      mScatter[7] = new NISTMottScatteringAngle(Element::N);
      mScatter[8] = new NISTMottScatteringAngle(Element::O);
      mScatter[9] = new NISTMottScatteringAngle(Element::F);
      mScatter[10] = new NISTMottScatteringAngle(Element::Ne);
      mScatter[11] = new NISTMottScatteringAngle(Element::Na);
      mScatter[12] = new NISTMottScatteringAngle(Element::Mg);
      mScatter[13] = new NISTMottScatteringAngle(Element::Al);
      mScatter[14] = new NISTMottScatteringAngle(Element::Si);
      mScatter[15] = new NISTMottScatteringAngle(Element::P);
      mScatter[16] = new NISTMottScatteringAngle(Element::S);
      mScatter[17] = new NISTMottScatteringAngle(Element::Cl);
      mScatter[18] = new NISTMottScatteringAngle(Element::Ar);
      mScatter[19] = new NISTMottScatteringAngle(Element::K);
      mScatter[20] = new NISTMottScatteringAngle(Element::Ca);
      mScatter[21] = new NISTMottScatteringAngle(Element::Sc);
      mScatter[22] = new NISTMottScatteringAngle(Element::Ti);
      mScatter[23] = new NISTMottScatteringAngle(Element::V);
      mScatter[24] = new NISTMottScatteringAngle(Element::Cr);
      mScatter[25] = new NISTMottScatteringAngle(Element::Mn);
      mScatter[26] = new NISTMottScatteringAngle(Element::Fe);
      mScatter[27] = new NISTMottScatteringAngle(Element::Co);
      mScatter[28] = new NISTMottScatteringAngle(Element::Ni);
      mScatter[29] = new NISTMottScatteringAngle(Element::Cu);
      mScatter[30] = new NISTMottScatteringAngle(Element::Zn);
      mScatter[31] = new NISTMottScatteringAngle(Element::Ga);
      mScatter[32] = new NISTMottScatteringAngle(Element::Ge);
      mScatter[33] = new NISTMottScatteringAngle(Element::As);
      mScatter[34] = new NISTMottScatteringAngle(Element::Se);
      mScatter[35] = new NISTMottScatteringAngle(Element::Br);
      mScatter[36] = new NISTMottScatteringAngle(Element::Kr);
      mScatter[37] = new NISTMottScatteringAngle(Element::Rb);
      mScatter[38] = new NISTMottScatteringAngle(Element::Sr);
      mScatter[39] = new NISTMottScatteringAngle(Element::Y);
      mScatter[40] = new NISTMottScatteringAngle(Element::Zr);
      mScatter[41] = new NISTMottScatteringAngle(Element::Nb);
      mScatter[42] = new NISTMottScatteringAngle(Element::Mo);
      mScatter[43] = new NISTMottScatteringAngle(Element::Tc);
      mScatter[44] = new NISTMottScatteringAngle(Element::Ru);
      mScatter[45] = new NISTMottScatteringAngle(Element::Rh);
      mScatter[46] = new NISTMottScatteringAngle(Element::Pd);
      mScatter[47] = new NISTMottScatteringAngle(Element::Ag);
      mScatter[48] = new NISTMottScatteringAngle(Element::Cd);
      mScatter[49] = new NISTMottScatteringAngle(Element::In);
      mScatter[50] = new NISTMottScatteringAngle(Element::Sn);
      mScatter[51] = new NISTMottScatteringAngle(Element::Sb);
      mScatter[52] = new NISTMottScatteringAngle(Element::Te);
      mScatter[53] = new NISTMottScatteringAngle(Element::I);
      mScatter[54] = new NISTMottScatteringAngle(Element::Xe);
      mScatter[55] = new NISTMottScatteringAngle(Element::Cs);
      mScatter[56] = new NISTMottScatteringAngle(Element::Ba);
      mScatter[57] = new NISTMottScatteringAngle(Element::La);
      mScatter[58] = new NISTMottScatteringAngle(Element::Ce);
      mScatter[59] = new NISTMottScatteringAngle(Element::Pr);
      mScatter[60] = new NISTMottScatteringAngle(Element::Nd);
      mScatter[61] = new NISTMottScatteringAngle(Element::Pm);
      mScatter[62] = new NISTMottScatteringAngle(Element::Sm);
      mScatter[63] = new NISTMottScatteringAngle(Element::Eu);
      mScatter[64] = new NISTMottScatteringAngle(Element::Gd);
      mScatter[65] = new NISTMottScatteringAngle(Element::Tb);
      mScatter[66] = new NISTMottScatteringAngle(Element::Dy);
      mScatter[67] = new NISTMottScatteringAngle(Element::Ho);
      mScatter[68] = new NISTMottScatteringAngle(Element::Er);
      mScatter[69] = new NISTMottScatteringAngle(Element::Tm);
      mScatter[70] = new NISTMottScatteringAngle(Element::Yb);
      mScatter[71] = new NISTMottScatteringAngle(Element::Lu);
      mScatter[72] = new NISTMottScatteringAngle(Element::Hf);
      mScatter[73] = new NISTMottScatteringAngle(Element::Ta);
      mScatter[74] = new NISTMottScatteringAngle(Element::W);
      mScatter[75] = new NISTMottScatteringAngle(Element::Re);
      mScatter[76] = new NISTMottScatteringAngle(Element::Os);
      mScatter[77] = new NISTMottScatteringAngle(Element::Ir);
      mScatter[78] = new NISTMottScatteringAngle(Element::Pt);
      mScatter[79] = new NISTMottScatteringAngle(Element::Au);
      mScatter[80] = new NISTMottScatteringAngle(Element::Hg);
      mScatter[81] = new NISTMottScatteringAngle(Element::Tl);
      mScatter[82] = new NISTMottScatteringAngle(Element::Pb);
      mScatter[83] = new NISTMottScatteringAngle(Element::Bi);
      mScatter[84] = new NISTMottScatteringAngle(Element::Po);
      mScatter[85] = new NISTMottScatteringAngle(Element::At);
      mScatter[86] = new NISTMottScatteringAngle(Element::Rn);
      mScatter[87] = new NISTMottScatteringAngle(Element::Fr);
      mScatter[88] = new NISTMottScatteringAngle(Element::Ra);
      mScatter[89] = new NISTMottScatteringAngle(Element::Ac);
      mScatter[90] = new NISTMottScatteringAngle(Element::Th);
      mScatter[91] = new NISTMottScatteringAngle(Element::Pa);
      mScatter[92] = new NISTMottScatteringAngle(Element::U);
      mScatter[93] = new NISTMottScatteringAngle(Element::Np);
      mScatter[94] = new NISTMottScatteringAngle(Element::Pu);
      mScatter[95] = new NISTMottScatteringAngle(Element::Am);
      mScatter[96] = new NISTMottScatteringAngle(Element::Cm);
   }

   __device__ NISTMottScatteringAngle *d_mScatter[113];

   __global__ void initCuda()
   {
      d_mScatter[1] = new NISTMottScatteringAngle(*Element::dH);
      d_mScatter[2] = new NISTMottScatteringAngle(*Element::dHe);
      d_mScatter[3] = new NISTMottScatteringAngle(*Element::dLi);
      d_mScatter[4] = new NISTMottScatteringAngle(*Element::dBe);
      d_mScatter[5] = new NISTMottScatteringAngle(*Element::dB);
      d_mScatter[6] = new NISTMottScatteringAngle(*Element::dC);
      d_mScatter[7] = new NISTMottScatteringAngle(*Element::dN);
      d_mScatter[8] = new NISTMottScatteringAngle(*Element::dO);
      d_mScatter[9] = new NISTMottScatteringAngle(*Element::dF);
      d_mScatter[10] = new NISTMottScatteringAngle(*Element::dNe);
      d_mScatter[11] = new NISTMottScatteringAngle(*Element::dNa);
      d_mScatter[12] = new NISTMottScatteringAngle(*Element::dMg);
      d_mScatter[13] = new NISTMottScatteringAngle(*Element::dAl);
      d_mScatter[14] = new NISTMottScatteringAngle(*Element::dSi);
      d_mScatter[15] = new NISTMottScatteringAngle(*Element::dP);
      d_mScatter[16] = new NISTMottScatteringAngle(*Element::dS);
      d_mScatter[17] = new NISTMottScatteringAngle(*Element::dCl);
      d_mScatter[18] = new NISTMottScatteringAngle(*Element::dAr);
      d_mScatter[19] = new NISTMottScatteringAngle(*Element::dK);
      d_mScatter[20] = new NISTMottScatteringAngle(*Element::dCa);
      d_mScatter[21] = new NISTMottScatteringAngle(*Element::dSc);
      d_mScatter[22] = new NISTMottScatteringAngle(*Element::dTi);
      d_mScatter[23] = new NISTMottScatteringAngle(*Element::dV);
      d_mScatter[24] = new NISTMottScatteringAngle(*Element::dCr);
      d_mScatter[25] = new NISTMottScatteringAngle(*Element::dMn);
      d_mScatter[26] = new NISTMottScatteringAngle(*Element::dFe);
      d_mScatter[27] = new NISTMottScatteringAngle(*Element::dCo);
      d_mScatter[28] = new NISTMottScatteringAngle(*Element::dNi);
      d_mScatter[29] = new NISTMottScatteringAngle(*Element::dCu);
      d_mScatter[30] = new NISTMottScatteringAngle(*Element::dZn);
      d_mScatter[31] = new NISTMottScatteringAngle(*Element::dGa);
      d_mScatter[32] = new NISTMottScatteringAngle(*Element::dGe);
      d_mScatter[33] = new NISTMottScatteringAngle(*Element::dAs);
      d_mScatter[34] = new NISTMottScatteringAngle(*Element::dSe);
      d_mScatter[35] = new NISTMottScatteringAngle(*Element::dBr);
      d_mScatter[36] = new NISTMottScatteringAngle(*Element::dKr);
      d_mScatter[37] = new NISTMottScatteringAngle(*Element::dRb);
      d_mScatter[38] = new NISTMottScatteringAngle(*Element::dSr);
      d_mScatter[39] = new NISTMottScatteringAngle(*Element::dY);
      d_mScatter[40] = new NISTMottScatteringAngle(*Element::dZr);
      d_mScatter[41] = new NISTMottScatteringAngle(*Element::dNb);
      d_mScatter[42] = new NISTMottScatteringAngle(*Element::dMo);
      d_mScatter[43] = new NISTMottScatteringAngle(*Element::dTc);
      d_mScatter[44] = new NISTMottScatteringAngle(*Element::dRu);
      d_mScatter[45] = new NISTMottScatteringAngle(*Element::dRh);
      d_mScatter[46] = new NISTMottScatteringAngle(*Element::dPd);
      d_mScatter[47] = new NISTMottScatteringAngle(*Element::dAg);
      d_mScatter[48] = new NISTMottScatteringAngle(*Element::dCd);
      d_mScatter[49] = new NISTMottScatteringAngle(*Element::dIn);
      d_mScatter[50] = new NISTMottScatteringAngle(*Element::dSn);
      d_mScatter[51] = new NISTMottScatteringAngle(*Element::dSb);
      d_mScatter[52] = new NISTMottScatteringAngle(*Element::dTe);
      d_mScatter[53] = new NISTMottScatteringAngle(*Element::dI);
      d_mScatter[54] = new NISTMottScatteringAngle(*Element::dXe);
      d_mScatter[55] = new NISTMottScatteringAngle(*Element::dCs);
      d_mScatter[56] = new NISTMottScatteringAngle(*Element::dBa);
      d_mScatter[57] = new NISTMottScatteringAngle(*Element::dLa);
      d_mScatter[58] = new NISTMottScatteringAngle(*Element::dCe);
      d_mScatter[59] = new NISTMottScatteringAngle(*Element::dPr);
      d_mScatter[60] = new NISTMottScatteringAngle(*Element::dNd);
      d_mScatter[61] = new NISTMottScatteringAngle(*Element::dPm);
      d_mScatter[62] = new NISTMottScatteringAngle(*Element::dSm);
      d_mScatter[63] = new NISTMottScatteringAngle(*Element::dEu);
      d_mScatter[64] = new NISTMottScatteringAngle(*Element::dGd);
      d_mScatter[65] = new NISTMottScatteringAngle(*Element::dTb);
      d_mScatter[66] = new NISTMottScatteringAngle(*Element::dDy);
      d_mScatter[67] = new NISTMottScatteringAngle(*Element::dHo);
      d_mScatter[68] = new NISTMottScatteringAngle(*Element::dEr);
      d_mScatter[69] = new NISTMottScatteringAngle(*Element::dTm);
      d_mScatter[70] = new NISTMottScatteringAngle(*Element::dYb);
      d_mScatter[71] = new NISTMottScatteringAngle(*Element::dLu);
      d_mScatter[72] = new NISTMottScatteringAngle(*Element::dHf);
      d_mScatter[73] = new NISTMottScatteringAngle(*Element::dTa);
      d_mScatter[74] = new NISTMottScatteringAngle(*Element::dW);
      d_mScatter[75] = new NISTMottScatteringAngle(*Element::dRe);
      d_mScatter[76] = new NISTMottScatteringAngle(*Element::dOs);
      d_mScatter[77] = new NISTMottScatteringAngle(*Element::dIr);
      d_mScatter[78] = new NISTMottScatteringAngle(*Element::dPt);
      d_mScatter[79] = new NISTMottScatteringAngle(*Element::dAu);
      d_mScatter[80] = new NISTMottScatteringAngle(*Element::dHg);
      d_mScatter[81] = new NISTMottScatteringAngle(*Element::dTl);
      d_mScatter[82] = new NISTMottScatteringAngle(*Element::dPb);
      d_mScatter[83] = new NISTMottScatteringAngle(*Element::dBi);
      d_mScatter[84] = new NISTMottScatteringAngle(*Element::dPo);
      d_mScatter[85] = new NISTMottScatteringAngle(*Element::dAt);
      d_mScatter[86] = new NISTMottScatteringAngle(*Element::dRn);
      d_mScatter[87] = new NISTMottScatteringAngle(*Element::dFr);
      d_mScatter[88] = new NISTMottScatteringAngle(*Element::dRa);
      d_mScatter[89] = new NISTMottScatteringAngle(*Element::dAc);
      d_mScatter[90] = new NISTMottScatteringAngle(*Element::dTh);
      d_mScatter[91] = new NISTMottScatteringAngle(*Element::dPa);
      d_mScatter[92] = new NISTMottScatteringAngle(*Element::dU);
      d_mScatter[93] = new NISTMottScatteringAngle(*Element::dNp);
      d_mScatter[94] = new NISTMottScatteringAngle(*Element::dPu);
      d_mScatter[95] = new NISTMottScatteringAngle(*Element::dAm);
      d_mScatter[96] = new NISTMottScatteringAngle(*Element::dCm);
   }

   //__device__ void NISTMottScatteringAngle::copySpwem(float *dSpwem, unsigned int size)
   //{
   //   mSpwem.assign(dSpwem, dSpwem + size);
   //}

   //__device__ void NISTMottScatteringAngle::copyX1j(unsigned int r, float *dSpwem, unsigned int size)
   //{
   //   mX1[r].assign(dSpwem, dSpwem + size);
   //}

   __global__ void assignSpwem(unsigned int i, float *dSpwem, unsigned int size)
   {
      d_mScatter[i]->assignSpwem<float>(dSpwem, size);
   }

   __global__ void assignX1Row(unsigned int i, unsigned int r, float *dSX1r, unsigned int size)
   {
      d_mScatter[i]->assignX1Row<float>(r, dSX1r, size);
   }

   void transferDataToCuda()
   {
      float **dSpwem = new float*[97];
      float ***dX1 = new float**[97];
      for (int i = 1; i <= 96; ++i) {
         const VectorXf& spwem = mScatter[i]->getSpwem();
         dSpwem[i] = new float[spwem.size()];
         checkCudaErrors(cudaMalloc((void **)&dSpwem[i], sizeof(float) * spwem.size()));
         checkCudaErrors(cudaMemcpy(dSpwem[i], spwem.data(), sizeof(float) * spwem.size(), cudaMemcpyHostToDevice));
         assignSpwem << <1, 1 >> >(i, dSpwem[i], spwem.size());

         const MatrixXf& x1 = mScatter[i]->getX1();
         dX1[i] = new float*[x1.size()];
         for (int r = 0; r < x1.size(); ++r) {
            dX1[i][r] = new float[x1[r].size()];
            checkCudaErrors(cudaMalloc((void **)&dX1[i][r], sizeof(float) * x1[r].size()));
            checkCudaErrors(cudaMemcpy(dX1[i][r], x1[r].data(), sizeof(float) * x1[r].size(), cudaMemcpyHostToDevice));
            assignX1Row << <1, 1 >> >(i, r, dX1[i][r], x1[r].size());
         }
      }
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      delete[] dSpwem;
      for (int i = 1; i <= 96; ++i) {
         delete[] dX1[i];
      }
      delete[] dX1;
   }

   __host__ __device__ const NISTMottScatteringAngle& getNISTMSA(int an)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      return *d_mScatter[an];
#else
      return *mScatter[an];
#endif
   }

   __host__ __device__ NISTMottRandomizedScatterFactory::NISTMottRandomizedScatterFactory() :
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      RandomizedScatterFactoryT("NIST Mott Inelastic Cross-Section", *Reference::d_NullReference)
#else
      RandomizedScatterFactoryT("NIST Mott Inelastic Cross-Section", REFERENCE)
#endif
   {
   }

   __host__ __device__ const RandomizedScatterT& NISTMottRandomizedScatterFactory::get(const ElementT& elm) const
   {
      return getNISTMSA(elm.getAtomicNumber());
   }

   const NISTMottRandomizedScatterFactory NISTMottRandomizedFactory;
   const RandomizedScatterFactoryT& Factory = NISTMottRandomizedFactory;
   __device__ const RandomizedScatterFactoryT* d_Factory = nullptr;

   __global__ void initFactory()
   {
      d_Factory = new NISTMottRandomizedScatterFactory();
   }
}
