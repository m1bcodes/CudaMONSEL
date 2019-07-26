#include "gov\nist\nanoscalemetrology\JMONSELTests\LinesOnLayers0.cuh"

#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\NISTMonte\NullMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"
#include "gov\nist\microanalysis\NISTMonte\GaussianBeam.cuh"
#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
#include "gov\nist\microanalysis\NISTMonte\Sphere.cuh"
#include "gov\nist\microanalysis\NISTMonte\BackscatterStats.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\Utility\Histogram.cuh"

#include "gov\nist\nanoscalemetrology\JMONSEL\SEmaterial.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ExpQMBarrierSM.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\MONSEL_MaterialScatterModel.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\SelectableElasticSM.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NISTMottRS.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\JoyLuoNieminenCSD.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\FittedInelSM.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\GanachaudMokraniPolaronTrapSM.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\TabulatedInelasticSM.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\GanachaudMokraniPhononInelasticSM.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalMultiPlaneShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalIntersectionShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NShapes.cuh"

#include "gov\nist\nanoscalemetrology\JMONSEL\NUTableInterpolation.cuh"

#include "Amphibian\random.cuh"

#include <fstream>

#include <chrono>

#include "CudaUtil.h"

namespace LinesOnLayers
{
   __global__ void verifyNUTable1d(const char* fn)
   {
      const NUTableInterpolationT* table = NUTableInterpolation::getInstance(fn);
      const VectorXd& data = table->gettable1d();
      printf("GPU %s\n", fn);
      for (auto v : data) {
         printf("%.5e ", v);
      }
      printf("\n");
   }

   __global__ void verifyNUTable2d(const char* fn, const int r)
   {
      const NUTableInterpolationT* table = NUTableInterpolation::getInstance(fn);
      const MatrixXd& data = table->gettable2d();
      printf("GPU %s: row %d\n", fn, r);
      for (auto v : data[r]) {
         printf("%.5e ", v);
      }
      printf("\n");
   }

   __global__ void verifyNUTable3d(const char* fn, const int r, const int c)
   {
      const NUTableInterpolationT* table = NUTableInterpolation::getInstance(fn);
      const Matrix3DXf& data = table->gettable3d();
      printf("GPU %s: row %d, col %d\n", fn, r, c);
      for (auto v : data[r][c]) {
         printf("%.5e ", v);
      }
      printf("\n");
   }

   void transferDataToCuda()
   {
      NUTableInterpolation::initFactory << <1, 1 >> >();
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());

      StringT tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\";
      const StringT IIMFPPennInterpglassy = tablePath + "IIMFPPennInterpglassyCSI.csv";
      const StringT SimReducedDeltaEglassy = tablePath + "interpNUSimReducedDeltaEglassyCSI.csv";
      const StringT simTableThetaNUglassy = tablePath + "interpsimTableThetaNUglassyCSI.csv";
      const StringT SimESE0NUglassy = tablePath + "interpSimESE0NUglassyCSI.csv";
      NUTableInterpolation::transferDataToCuda(IIMFPPennInterpglassy.c_str());
      NUTableInterpolation::transferDataToCuda(SimReducedDeltaEglassy.c_str());
      NUTableInterpolation::transferDataToCuda(simTableThetaNUglassy.c_str());
      NUTableInterpolation::transferDataToCuda(SimESE0NUglassy.c_str());

      //tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\";
      //const StringT IIMFPFullPennInterpSiSI = tablePath + "IIMFPFullPennInterpSiSI.csv";
      //const StringT interpNUSimReducedDeltaEFullPennSiSI = tablePath + "interpNUSimReducedDeltaEFullPennSiSI.csv";
      //const StringT interpNUThetaFullPennSiBGSI = tablePath + "interpNUThetaFullPennSiBGSI.csv";
      //const StringT interpSimESE0NUSiBGSI = tablePath + "interpSimESE0NUSiBGSI.csv";
      //NUTableInterpolation::copyDataToCuda(IIMFPFullPennInterpSiSI.c_str());
      //NUTableInterpolation::copyDataToCuda(interpNUSimReducedDeltaEFullPennSiSI.c_str());
      //NUTableInterpolation::copyDataToCuda(interpNUThetaFullPennSiBGSI.c_str());
      //NUTableInterpolation::copyDataToCuda(interpSimESE0NUSiBGSI.c_str());

      const char* fn = IIMFPPennInterpglassy.c_str();
      char* d_fn = nullptr;

      checkCudaErrors(cudaMalloc((void **)&d_fn, (IIMFPPennInterpglassy.size() + 1) * sizeof(char)));
      checkCudaErrors(cudaMemset(d_fn, 0, (IIMFPPennInterpglassy.size() + 1) * sizeof(char)));
      checkCudaErrors(cudaMemcpy(d_fn, fn, (IIMFPPennInterpglassy.size() + 1) * sizeof(char), cudaMemcpyHostToDevice));
      verifyNUTable1d << <1, 1 >> >(d_fn);
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      checkCudaErrors(cudaFree(d_fn));
      const NUTableInterpolationT* table0 = NUTableInterpolation::getInstance(fn);
      const VectorXd& data0 = table0->gettable1d();
      printf("CPU %s\n", fn);
      for (auto v : data0) {
         printf("%.5e ", v);
      }
      printf("\n");

      fn = SimReducedDeltaEglassy.c_str();
      checkCudaErrors(cudaMalloc((void **)&d_fn, (SimReducedDeltaEglassy.size() + 1) * sizeof(char)));
      checkCudaErrors(cudaMemset(d_fn, 0, (SimReducedDeltaEglassy.size() + 1) * sizeof(char)));
      checkCudaErrors(cudaMemcpy(d_fn, fn, (SimReducedDeltaEglassy.size() + 1) * sizeof(char), cudaMemcpyHostToDevice));
      verifyNUTable2d << <1, 1 >> >(d_fn, 0);
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      checkCudaErrors(cudaFree(d_fn));
      const NUTableInterpolationT* table1 = NUTableInterpolation::getInstance(fn);
      const MatrixXd& data1 = table1->gettable2d();
      printf("CPU %s: row %d\n", fn, 0);
      for (auto v : data1[0]) {
         printf("%.5e ", v);
      }
      printf("\n");

      fn = simTableThetaNUglassy.c_str();
      checkCudaErrors(cudaMalloc((void **)&d_fn, (simTableThetaNUglassy.size() + 1) * sizeof(char)));
      checkCudaErrors(cudaMemset(d_fn, 0, (simTableThetaNUglassy.size() + 1) * sizeof(char)));
      checkCudaErrors(cudaMemcpy(d_fn, fn, (simTableThetaNUglassy.size() + 1) * sizeof(char), cudaMemcpyHostToDevice));
      verifyNUTable3d << <1, 1 >> >(d_fn, 50, 50);
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      checkCudaErrors(cudaFree(d_fn));
      const NUTableInterpolationT* table2 = NUTableInterpolation::getInstance(fn);
      const Matrix3DXf& data2 = table2->gettable3d();
      printf("CPU %s: row %d, col %d\n", fn, 50, 50);
      for (auto v : data2[50][50]) {
         printf("%.5e ", v);
      }
      printf("\n");

      tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\";
      const StringT IIMFPFullPennInterpSiSI = tablePath + "IIMFPFullPennInterpSiSI.csv";
      const StringT interpNUSimReducedDeltaEFullPennSiSI = tablePath + "interpNUSimReducedDeltaEFullPennSiSI.csv";
      const StringT interpNUThetaFullPennSiBGSI = tablePath + "interpNUThetaFullPennSiBGSI.csv";
      const StringT interpSimESE0NUSiBGSI = tablePath + "interpSimESE0NUSiBGSI.csv";
      NUTableInterpolation::transferDataToCuda(IIMFPFullPennInterpSiSI.c_str());
      NUTableInterpolation::transferDataToCuda(interpNUSimReducedDeltaEFullPennSiSI.c_str());
      NUTableInterpolation::transferDataToCuda(interpNUThetaFullPennSiBGSI.c_str());
      NUTableInterpolation::transferDataToCuda(interpSimESE0NUSiBGSI.c_str());
   }

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const int nTrajectories = 100;

   __constant__ const double pitchnm = 180;
   __constant__ const int nlines = 3;
   __constant__ const double hnm = 120;
   __constant__ const double wnm = 80;
   __constant__ const double linelengthnm = 1000;
   __constant__ const double thetardeg = 3;
   __constant__ const double thetaldeg = 3;
   __constant__ const double radrnm = 20;
   __constant__ const double radlnm = 20;
   __constant__ const double layer1thicknessnm = 80;
   __constant__ const double layer2thicknessnm = 200;

   __constant__ const double beamEeVvals[] = { 500. };
   __constant__ const int beamEeVvalsLen = 1;
   __constant__ const double beamsizenm = 0.5;
   __constant__ const double deepnm = 15;

   __constant__ const bool trajImg = true;
   __constant__ const int trajImgMaxTraj = 50;
   __constant__ const double trajImgSize = 200e-9;

   __constant__ const bool VRML = false;
   __constant__ const int VRMLImgMaxTraj = 0;

   //__device__ SEmaterialT* vacuum = nullptr;
   //__device__ ExpQMBarrierSMT* vacuumBarrier = nullptr;
   //__device__ ZeroCSDT* sZeroCSD = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* vacuumMSM = nullptr;

   __constant__ const double PMMAbreakE = 1.60217653e-19 * 45.;
   __constant__ const double PMMAdensity = 1190.;
   __constant__ const double PMMAworkfun = 5.5;
   __constant__ const double PMMAbandgap = 5.;
   __constant__ const double PMMAEFermi = -5.;//-PMMAbandgap;
   __constant__ const double PMMApotU = -5.5 - (-5.);

   //__device__ SEmaterialT* PMMA = nullptr;

   //__device__ SelectableElasticSMT* PMMANISTMott = nullptr;

   //__device__ JoyLuoNieminenCSDT* PMMACSD = nullptr;
   //__device__ FittedInelSMT* PMMAfittedInel = nullptr;
   //__device__ GanachaudMokraniPolaronTrapSMT* PMMApolaron = nullptr;

   //__device__ ExpQMBarrierSMT* pmmaeqmbsm = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* PMMAMSM = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* PMMAMSMDeep = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* ARCMSM = nullptr;

   __constant__ const double glCdensity = 1800.;
   __constant__ const double glCworkfun = 5.0;
   __constant__ const double glCbandgap = 0.;
   __constant__ const double glCEFermi = 20.4;
   __constant__ const double glCpotU = -5. - 20.4;

   //__device__ SEmaterialT* glC = nullptr;

   //__device__ SelectableElasticSMT* glCNISTMott = nullptr;

   //__device__ TabulatedInelasticSMT* glCDS = nullptr;

   //__device__ ExpQMBarrierSMT* glceqmbsm = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* glCMSM = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* glCMSMDeep = nullptr;

   __constant__ const double phononE = 0.063;
   __constant__ const double phononStrength = 3.;

   __constant__ const double Sidensity = 2330.;
   __constant__ const double Siworkfun = 4.85;
   __constant__ const double Sibandgap = 1.1;
   __constant__ const double SiEFermi = -1.1;//-Sibandgap;
   __constant__ const double SipotU = -44.85 - (-1.1);//-Siworkfun - SiEFermi;

   //__device__ SEmaterialT* Si = nullptr;

   //__device__ SelectableElasticSMT* SiNISTMott = nullptr;

   //__device__ TabulatedInelasticSMT* SiDS = nullptr;

   //__device__ GanachaudMokraniPhononInelasticSMT* Siphonon = nullptr;

   //__device__ ExpQMBarrierSMT* sieqmbsm = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* SiMSM = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* SiMSMDeep = nullptr;

   //__device__ SphereT* sphere = nullptr;

   __constant__ const double pitch = 180 * 1.e-9;
   __constant__ const double h = 120 * 1.e-9;
   __constant__ const double w = 80 * 1.e-9;
   __constant__ const double linelength = 1000 * 1.e-9;

   __constant__ const double radperdeg = 3.14159265358979323846 / 180.;
   __constant__ const double thetar = 3 * 3.14159265358979323846 / 180.;
   __constant__ const double thetal = 3 * 3.14159265358979323846 / 180.;
   __constant__ const double radr = 20 * 1.e-9;
   __constant__ const double radl = 20 * 1.e-9;
   __constant__ const double layer1thickness = 80 * 1.e-9;
   __constant__ const double layer2thickness = 200 * 1.e-9;
   __constant__ const double beamsize = 0.5 * 1.e-9;
   __constant__ const double deep = 15 * 1.e-9;

   __constant__ const double center[] = {
      0.0,
      0.0,
      0.0
   };

   __constant__ const double beamEeV = 500.;
   __constant__ const double beamE = 1.60217653e-19 * 500;
   __constant__ const double binSizeEV = 10.;

   //__device__ NullMaterialScatterModelT* NULL_MSM = nullptr;

   //__device__ RegionT* chamber = nullptr;

   __constant__ const double normalvector[] = { 0., 0., -1. };
   __constant__ const double layer1Pos[] = { 0., 0., 0. };

   //__device__ NormalMultiPlaneShapeT* layer1 = nullptr;
   //__device__ PlaneT* pl1 = nullptr;
   //__device__ RegionT* layer1Region = nullptr;

   __constant__ const double layer2Pos[] = { 0., 0., 80 * 1.e-9 };
   //__device__ NormalMultiPlaneShapeT* layer2 = nullptr;
   //__device__ PlaneT* pl2 = nullptr;
   //__device__ RegionT* layer2Region = nullptr;

   __constant__ const double layer3Pos[] = { 0., 0., 80 * 1.e-9 + 200 * 1.e-9 };
   //__device__ NormalMultiPlaneShapeT* layer3 = nullptr;
   //__device__ PlaneT* pl3 = nullptr;
   //__device__ RegionT* layer3Region = nullptr;

   __constant__ const double layer4Pos[] = { 0., 0., 80 * 1.e-9 + 200 * 1.e-9 + 15 * 1.e-9 };
   //__device__ NormalMultiPlaneShapeT* layer4 = nullptr;
   //__device__ PlaneT* pl4 = nullptr;
   //__device__ RegionT* layer4Region = nullptr;

   //__device__ RegionT* deepRegion = nullptr;

   __constant__ const double leftmostLineCenterx = -180. * 1.e-9 * (3. / 2.);
   __constant__ const double xcenter = -180. * 1.e-9 * (3. / 2.) + 0 * 180 * 1.e-9;

   //__device__ NormalIntersectionShapeT* line = nullptr;
   //__device__ RegionT* lineRegion = nullptr;

   __device__ double* yvals = nullptr;
   __device__ double* xvals = nullptr;
   __device__ unsigned int yvalsSize = 0;
   __device__ unsigned int xvalsSize = 0;
#else
   const int nTrajectories = 100;

   const double pitchnm = 180;
   const int nlines = 3;
   const double hnm = 120;
   const double wnm = 80;
   const double linelengthnm = 1000;
   const double thetardeg = 3;
   const double thetaldeg = 3;
   const double radrnm = 20;
   const double radlnm = 20;
   const double layer1thicknessnm = 80;
   const double layer2thicknessnm = 200;

   const double beamEeVvals[] = { 500. };
   const int beamEeVvalsLen = 1;
   const double beamsizenm = 0.5;
   const double deepnm = 15;

   const bool trajImg = true;
   const int trajImgMaxTraj = 50;
   const double trajImgSize = 200e-9;

   const bool VRML = false;
   const int VRMLImgMaxTraj = 0;

   //SEmaterialT* vacuum = nullptr;
   //ExpQMBarrierSMT* vacuumBarrier = nullptr;
   //ZeroCSDT* sZeroCSD = nullptr;

   //MONSEL_MaterialScatterModelT* vacuumMSM = nullptr;

   const double PMMAbreakE = 1.60217653e-19 * 45.;
   const double PMMAdensity = 1190.;
   const double PMMAworkfun = 5.5;
   const double PMMAbandgap = 5.;
   const double PMMAEFermi = -5.;//-PMMAbandgap;
   const double PMMApotU = -5.5 - (-5.);

   //SEmaterialT* PMMA = nullptr;

   //SelectableElasticSMT* PMMANISTMott = nullptr;

   //JoyLuoNieminenCSDT* PMMACSD = nullptr;
   //FittedInelSMT* PMMAfittedInel = nullptr;
   //GanachaudMokraniPolaronTrapSMT* PMMApolaron = nullptr;

   //ExpQMBarrierSMT* pmmaeqmbsm = nullptr;

   //MONSEL_MaterialScatterModelT* PMMAMSM = nullptr;

   //MONSEL_MaterialScatterModelT* PMMAMSMDeep = nullptr;

   //MONSEL_MaterialScatterModelT* ARCMSM = nullptr;

   const double glCdensity = 1800.;
   const double glCworkfun = 5.0;
   const double glCbandgap = 0.;
   const double glCEFermi = 20.4;
   const double glCpotU = -5. - 20.4;

   //SEmaterialT* glC = nullptr;

   //SelectableElasticSMT* glCNISTMott = nullptr;

   //TabulatedInelasticSMT* glCDS = nullptr;

   //ExpQMBarrierSMT* glceqmbsm = nullptr;

   //MONSEL_MaterialScatterModelT* glCMSM = nullptr;

   //MONSEL_MaterialScatterModelT* glCMSMDeep = nullptr;

   const double phononE = 0.063;
   const double phononStrength = 3.;

   const double Sidensity = 2330.;
   const double Siworkfun = 4.85;
   const double Sibandgap = 1.1;
   const double SiEFermi = -1.1;//-Sibandgap;
   const double SipotU = -44.85 - (-1.1);//-Siworkfun - SiEFermi;

   //SEmaterialT* Si = nullptr;

   //SelectableElasticSMT* SiNISTMott = nullptr;

   //TabulatedInelasticSMT* SiDS = nullptr;

   //GanachaudMokraniPhononInelasticSMT* Siphonon = nullptr;

   //ExpQMBarrierSMT* sieqmbsm = nullptr;

   //MONSEL_MaterialScatterModelT* SiMSM = nullptr;

   //MONSEL_MaterialScatterModelT* SiMSMDeep = nullptr;

   //SphereT* sphere = nullptr;

   const double pitch = 180 * 1.e-9;
   const double h = 120 * 1.e-9;
   const double w = 80 * 1.e-9;
   const double linelength = 1000 * 1.e-9;

   const double radperdeg = 3.14159265358979323846 / 180.;
   const double thetar = 3 * 3.14159265358979323846 / 180.;
   const double thetal = 3 * 3.14159265358979323846 / 180.;
   const double radr = 20 * 1.e-9;
   const double radl = 20 * 1.e-9;
   const double layer1thickness = 80 * 1.e-9;
   const double layer2thickness = 200 * 1.e-9;
   const double beamsize = 0.5 * 1.e-9;
   const double deep = 15 * 1.e-9;

   const double center[] = {
      0.0,
      0.0,
      0.0
   };

   const double beamEeV = 500.;
   const double beamE = 1.60217653e-19 * 500;
   const double binSizeEV = 10.;

   //NullMaterialScatterModelT* NULL_MSM = nullptr;

   //RegionT* chamber = nullptr;

   const double normalvector[] = { 0., 0., -1. };
   const double layer1Pos[] = { 0., 0., 0. };

   //NormalMultiPlaneShapeT* layer1 = nullptr;
   //PlaneT* pl1 = nullptr;
   //RegionT* layer1Region = nullptr;

   const double layer2Pos[] = { 0., 0., 80 * 1.e-9 };
   //NormalMultiPlaneShapeT* layer2 = nullptr;
   //PlaneT* pl2 = nullptr;
   //RegionT* layer2Region = nullptr;

   const double layer3Pos[] = { 0., 0., 80 * 1.e-9 + 200 * 1.e-9 };
   //NormalMultiPlaneShapeT* layer3 = nullptr;
   //PlaneT* pl3 = nullptr;
   //RegionT* layer3Region = nullptr;

   const double layer4Pos[] = { 0., 0., 80 * 1.e-9 + 200 * 1.e-9 + 15 * 1.e-9 };
   //NormalMultiPlaneShapeT* layer4 = nullptr;
   //PlaneT* pl4 = nullptr;
   //RegionT* layer4Region = nullptr;

   //RegionT* deepRegion = nullptr;

   const double leftmostLineCenterx = -180. * 1.e-9 * (3. / 2.);
   const double xcenter = -180. * 1.e-9 * (3. / 2.) + 0 * 180 * 1.e-9;

   //NormalIntersectionShapeT* line = nullptr;
   //RegionT* lineRegion = nullptr;

   double* yvals = nullptr;
   double* xvals = nullptr;
   unsigned int yvalsSize = 0;
   unsigned int xvalsSize = 0;
#endif

   __global__ void initCuda()
   //void initCuda()
   {
      printf("LinesOnLayers: initCuda\n");
      for (int i = 0; i < 10; ++i) {
         printf("%.10e\n", Random::random());
      }

      //vacuum = new SEmaterialT();
      //vacuum->setName("SE vacuum");
      //vacuumBarrier = new ExpQMBarrierSMT(vacuum);
      //sZeroCSD = new ZeroCSDT();

      //vacuumMSM = new MONSEL_MaterialScatterModelT(vacuum, vacuumBarrier, sZeroCSD);

//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      const ElementT* componentsCOH[] = { Element::dC, Element::dO, Element::dH };
//#else
//      const ElementT* componentsCOH[] = { &Element::C, &Element::O, &Element::H };
//#endif
      //CompositionT* PMMAcomp = new CompositionT();
      //const double compositionCOH[] = { 5, 2, 8 };
      //PMMAcomp->defineByMoleFraction(componentsCOH, 3, compositionCOH, 3);
      //SEmaterialT* PMMA = new SEmaterialT(*PMMAcomp, PMMAdensity);
      //PMMA->setName("PMMA");
      //PMMA->setWorkfunction(ToSI::eV(PMMAworkfun));
      //PMMA->setBandgap(ToSI::eV(PMMAbandgap));
      //PMMA->setEnergyCBbottom(ToSI::eV(PMMApotU));

//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      PMMANISTMott = new SelectableElasticSMT(*PMMA, *NISTMottRS::d_Factory);
//#else
//      PMMANISTMott = new SelectableElasticSMT(*PMMA, NISTMottRS::Factory);
//#endif

      //PMMACSD = new JoyLuoNieminenCSDT(*PMMA, PMMAbreakE);
      //PMMAfittedInel = new FittedInelSMT(*PMMA, ToSI::eV(65.4), *PMMACSD);
      //PMMApolaron = new GanachaudMokraniPolaronTrapSMT(2.e7, 1. / ToSI::eV(4.));

      //pmmaeqmbsm = new ExpQMBarrierSMT(PMMA);

      //PMMAMSM = new MONSEL_MaterialScatterModelT(PMMA, pmmaeqmbsm, sZeroCSD);
      //PMMAMSM->addScatterMechanism(PMMANISTMott);
      //PMMAMSM->addScatterMechanism(PMMAfittedInel);
      //PMMAMSM->addScatterMechanism(PMMApolaron);

      //PMMAMSM->setCSD(PMMACSD);

      //PMMAMSMDeep = new MONSEL_MaterialScatterModelT(PMMA, pmmaeqmbsm, sZeroCSD);
      //PMMAMSMDeep->addScatterMechanism(PMMANISTMott);
      //PMMAMSMDeep->addScatterMechanism(PMMAfittedInel);
      //PMMAMSMDeep->addScatterMechanism(PMMApolaron);

      //PMMAMSMDeep->setCSD(PMMACSD);
      //PMMAMSMDeep->setMinEforTracking(ToSI::eV(50.));

      //ARCMSM = PMMAMSM;

//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      const ElementT* glCComponents[] = { Element::dC };
//#else
//      const ElementT* glCComponents[] = { &Element::C };
//#endif
//      const double glCComposition[] = { 1. };
//      glC = new SEmaterialT(glCComponents, 1, glCComposition, 1, glCdensity, "glassy Carbon");
//      glC->setWorkfunction(ToSI::eV(glCworkfun));
//      glC->setEnergyCBbottom(ToSI::eV(glCpotU));
//      glC->setBandgap(ToSI::eV(glCbandgap));
//      const double glCCoreEnergy[] = { ToSI::eV(284.2) };
//      glC->setCoreEnergy(glCCoreEnergy, 1);
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      glCNISTMott = new SelectableElasticSMT(*glC, *NISTMottRS::d_Factory);
//#else
//      glCNISTMott = new SelectableElasticSMT(*glC, NISTMottRS::Factory);
//#endif

      //StringT tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\";
      //const StringT IIMFPPennInterpglassy = tablePath + "IIMFPPennInterpglassyCSI.csv";
      //const StringT SimReducedDeltaEglassy = tablePath + "interpNUSimReducedDeltaEglassyCSI.csv";
      //const StringT simTableThetaNUglassy = tablePath + "interpsimTableThetaNUglassyCSI.csv";
      //const StringT SimESE0NUglassy = tablePath + "interpSimESE0NUglassyCSI.csv";
      //const char* glCTables[] = {
      //   IIMFPPennInterpglassy.c_str(),
      //   SimReducedDeltaEglassy.c_str(),
      //   simTableThetaNUglassy.c_str(),
      //   SimESE0NUglassy.c_str()
      //};

      //glCDS = new TabulatedInelasticSMT(*glC, 3, glCTables);

      //glceqmbsm = new ExpQMBarrierSMT(glC);

      //glCMSM = new MONSEL_MaterialScatterModelT(glC, glceqmbsm, sZeroCSD);
      //glCMSM->addScatterMechanism(glCNISTMott);
      //glCMSM->addScatterMechanism(glCDS);

      //glCMSMDeep = new MONSEL_MaterialScatterModelT(glC, glceqmbsm, sZeroCSD);
      //glCMSMDeep->addScatterMechanism(glCNISTMott);
      //glCMSMDeep->addScatterMechanism(glCDS);

      //glCMSMDeep->setMinEforTracking(ToSI::eV(50.));

//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      const ElementT* SiComponent[] = { Element::dSi };
//#else
//      const ElementT* SiComponent[] = { &Element::Si };
//#endif
      //const double SiComposition[] = { 1. };
      //Si = new SEmaterialT(SiComponent, 1, SiComposition, 1, Sidensity, "Silicon");
      //Si->setWorkfunction(ToSI::eV(Siworkfun));
      //Si->setEnergyCBbottom(ToSI::eV(SipotU));
      //Si->setBandgap(ToSI::eV(Sibandgap));
      //const double SiCoreEnergy[] = { ToSI::eV(99.2), ToSI::eV(99.8), ToSI::eV(149.7), ToSI::eV(1839.) };
      //Si->setCoreEnergy(SiCoreEnergy, 4);

//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      SiNISTMott = new SelectableElasticSMT(*Si, *NISTMottRS::d_Factory);
//#else
//      SiNISTMott = new SelectableElasticSMT(*Si, NISTMottRS::Factory);
//#endif

      //tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\";
      //const StringT IIMFPFullPennInterpSiSI = tablePath + "IIMFPFullPennInterpSiSI.csv";
      //const StringT interpNUSimReducedDeltaEFullPennSiSI = tablePath + "interpNUSimReducedDeltaEFullPennSiSI.csv";
      //const StringT interpNUThetaFullPennSiBGSI = tablePath + "interpNUThetaFullPennSiBGSI.csv";
      //const StringT interpSimESE0NUSiBGSI = tablePath + "interpSimESE0NUSiBGSI.csv";
      //const char* SiTables[] = {
      //   IIMFPFullPennInterpSiSI.c_str(),
      //   interpNUSimReducedDeltaEFullPennSiSI.c_str(),
      //   interpNUThetaFullPennSiBGSI.c_str(),
      //   interpSimESE0NUSiBGSI.c_str()
      //};

      //SiDS = new TabulatedInelasticSMT(*Si, 3, SiTables, ToSI::eV(13.54));

      //Siphonon = new GanachaudMokraniPhononInelasticSMT(phononStrength, ToSI::eV(phononE), 300., 11.7, 1.);

      //sieqmbsm = new ExpQMBarrierSMT(Si);

      //SiMSM = new MONSEL_MaterialScatterModelT(Si, sieqmbsm, sZeroCSD);
      //SiMSM->addScatterMechanism(SiNISTMott);
      //SiMSM->addScatterMechanism(SiDS);
      //SiMSM->addScatterMechanism(Siphonon);

      //SiMSMDeep = new MONSEL_MaterialScatterModelT(Si, sieqmbsm, sZeroCSD);
      //SiMSMDeep->addScatterMechanism(SiNISTMott);
      //SiMSMDeep->addScatterMechanism(SiDS);
      //SiMSMDeep->addScatterMechanism(Siphonon);

      //SiMSMDeep->setMinEforTracking(ToSI::eV(50.));

      //sphere = new SphereT(center, MonteCarloSS::ChamberRadius);

      //NULL_MSM = new NullMaterialScatterModelT();
      //chamber = new RegionT(nullptr, NULL_MSM, sphere);
      //chamber->updateMaterial(*(chamber->getScatterModel()), *vacuumMSM);

      //layer1 = new NormalMultiPlaneShapeT();
      //pl1 = new PlaneT(normalvector, layer1Pos);
      //layer1->addPlane(pl1);
      //layer1Region = new RegionT(chamber, ARCMSM, (NormalShapeT*)layer1);

      //layer2 = new NormalMultiPlaneShapeT();
      //pl2 = new PlaneT(normalvector, layer2Pos);
      //layer2->addPlane(pl2);
      //layer2Region = new RegionT(layer1Region, glCMSM, (NormalShapeT*)layer2);

      //layer3 = new NormalMultiPlaneShapeT();
      //pl3 = new PlaneT(normalvector, layer2Pos);
      //layer3->addPlane(pl3);
      //layer3Region = new RegionT(layer2Region, SiMSM, (NormalShapeT*)layer3);

      //layer4 = new NormalMultiPlaneShapeT();
      //pl4 = new PlaneT(normalvector, layer4Pos);
      //layer4Region = new RegionT(layer3Region, SiMSM, (NormalShapeT*)layer4);

      //deepRegion = new RegionT(layer3Region, SiMSMDeep, (NormalShapeT*)layer4);

      //for (int i = 0; i < nlines; ++i) {
      //   double xcenter = leftmostLineCenterx + i*pitch;
      //   NormalIntersectionShapeT* line = (NormalIntersectionShapeT*)NShapes::createLine(-h, w, linelength, thetal, thetar, radl, radr);
      //   const double newLinePos[] = { xcenter, 0., 0. };
      //   line->translate(newLinePos);
      //   RegionT lineRegion(&chamber, &PMMAMSM, line);
      //}

      //line = (NormalIntersectionShapeT*)NShapes::createLine(-h, w, linelength, thetal, thetar, radl, radr);
      //lineRegion = new RegionT(chamber, PMMAMSM, line);

      VectorXd yvalstmp(128);
      for (int i = -64; i < 64; i += 1) {
         yvalstmp.push_back(i);
      }

      const double xbottom = wnm / 2.;
      const double xtop = wnm / 2. - hnm * ::tan(thetar);
      const double xstart = xbottom - 100.5;
      const double xstop = xbottom + 100.5;
      const double xfinestart = xtop - 20.5;
      double xfinestop;
      if (thetar < 0.) xfinestop = xtop + 20.5;
      else xfinestop = wnm / 2. + 20.5;

      VectorXd xvalstmp(128);
      double deltax = 5.;
      double x = xstart;
      while (x < xfinestart) {
         xvalstmp.push_back(x);
         x += deltax;
      }
      x = xfinestart;
      deltax = 1;
      while (x < xfinestop) {
         xvalstmp.push_back(x);
         x += deltax;
      }
      x = xfinestop;
      deltax = 5.;
      while (x < xstop) {
         xvalstmp.push_back(x);
         x += deltax;
      }
      xvalstmp.push_back(xstop);

      yvalsSize = yvalstmp.size();
      xvalsSize = xvalstmp.size();

      yvals = new double[yvalsSize];
      xvals = new double[xvalsSize];

      memcpy(yvals, yvalstmp.data(), yvalsSize * sizeof(double));
      memcpy(xvals, xvalstmp.data(), xvalsSize * sizeof(double));
   }

   __global__ void runCudaSinglePixel()
   {
      int r = blockIdx.y*blockDim.y + threadIdx.y;
      int c = blockIdx.x*blockDim.x + threadIdx.x;

      int blockId = blockIdx.x + blockIdx.y * gridDim.x;
      int threadId = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;
      printf("%d, %d (%d) began\n", r, c, threadId);

      const double ynm = yvals[r];
      const double y = ynm * 1.e-9;
      const double xnm = xvals[c];
      const double x = xnm * 1.e-9;

      SEmaterialT vacuum;
      vacuum.setName("SE vacuum");
      ExpQMBarrierSMT vacuumBarrier(&vacuum);
      ZeroCSDT sZeroCSD;

      MONSEL_MaterialScatterModelT vacuumMSM(&vacuum, &vacuumBarrier, &sZeroCSD);

      const ElementT* componentsCOH[] = { Element::dC, Element::dO, Element::dH };

      CompositionT PMMAcomp;
      const double compositionCOH[] = { 5, 2, 8 };
      PMMAcomp.defineByMoleFraction(componentsCOH, 3, compositionCOH, 3);
      SEmaterialT PMMA(PMMAcomp, PMMAdensity);
      PMMA.setName("PMMA");
      PMMA.setWorkfunction(ToSI::eV(PMMAworkfun));
      PMMA.setBandgap(ToSI::eV(PMMAbandgap));
      PMMA.setEnergyCBbottom(ToSI::eV(PMMApotU));

      SelectableElasticSMT PMMANISTMott(PMMA, *NISTMottRS::d_Factory);

      JoyLuoNieminenCSDT PMMACSD(PMMA, PMMAbreakE);
      FittedInelSMT PMMAfittedInel(PMMA, ToSI::eV(65.4), PMMACSD);
      GanachaudMokraniPolaronTrapSMT PMMApolaron(2.e7, 1. / ToSI::eV(4.));

      ExpQMBarrierSMT pmmaeqmbsm(&PMMA);

      MONSEL_MaterialScatterModelT PMMAMSM(&PMMA, &pmmaeqmbsm, &sZeroCSD);
      PMMAMSM.addScatterMechanism(&PMMANISTMott);
      PMMAMSM.addScatterMechanism(&PMMAfittedInel);
      PMMAMSM.addScatterMechanism(&PMMApolaron);

      PMMAMSM.setCSD(&PMMACSD);

      MONSEL_MaterialScatterModelT PMMAMSMDeep(&PMMA, &pmmaeqmbsm, &sZeroCSD);
      PMMAMSMDeep.addScatterMechanism(&PMMANISTMott);
      PMMAMSMDeep.addScatterMechanism(&PMMAfittedInel);
      PMMAMSMDeep.addScatterMechanism(&PMMApolaron);

      PMMAMSMDeep.setCSD(&PMMACSD);
      PMMAMSMDeep.setMinEforTracking(ToSI::eV(50.));

      MONSEL_MaterialScatterModelT& ARCMSM = PMMAMSM;

      const ElementT* glCComponents[] = { Element::dC };
      const double glCComposition[] = { 1. };
      SEmaterialT glC(glCComponents, 1, glCComposition, 1, glCdensity, "glassy Carbon");
      glC.setWorkfunction(ToSI::eV(glCworkfun));
      glC.setEnergyCBbottom(ToSI::eV(glCpotU));
      glC.setBandgap(ToSI::eV(glCbandgap));
      const double glCCoreEnergy[] = { ToSI::eV(284.2) };
      glC.setCoreEnergy(glCCoreEnergy, 1);

      SelectableElasticSMT glCNISTMott(glC, *NISTMottRS::d_Factory);

      StringT tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\";
      const StringT IIMFPPennInterpglassy = tablePath + "IIMFPPennInterpglassyCSI.csv";
      const StringT SimReducedDeltaEglassy = tablePath + "interpNUSimReducedDeltaEglassyCSI.csv";
      const StringT simTableThetaNUglassy = tablePath + "interpsimTableThetaNUglassyCSI.csv";
      const StringT SimESE0NUglassy = tablePath + "interpSimESE0NUglassyCSI.csv";
      const char* glCTables[] = {
         IIMFPPennInterpglassy.c_str(),
         SimReducedDeltaEglassy.c_str(),
         simTableThetaNUglassy.c_str(),
         SimESE0NUglassy.c_str()
      };
      TabulatedInelasticSMT glCDS(glC, 3, glCTables);

      ExpQMBarrierSMT glceqmbsm(&glC);

      MONSEL_MaterialScatterModelT glCMSM(&glC, &glceqmbsm, &sZeroCSD);
      glCMSM.addScatterMechanism(&glCNISTMott);
      glCMSM.addScatterMechanism(&glCDS);

      MONSEL_MaterialScatterModelT glCMSMDeep(&glC, &glceqmbsm, &sZeroCSD);
      glCMSMDeep.addScatterMechanism(&glCNISTMott);
      glCMSMDeep.addScatterMechanism(&glCDS);

      glCMSMDeep.setMinEforTracking(ToSI::eV(50.));

      const ElementT* SiComponent[] = { Element::dSi };
      const double SiComposition[] = { 1. };
      SEmaterialT Si(SiComponent, 1, SiComposition, 1, Sidensity, "Silicon");
      Si.setWorkfunction(ToSI::eV(Siworkfun));
      Si.setEnergyCBbottom(ToSI::eV(SipotU));
      Si.setBandgap(ToSI::eV(Sibandgap));
      const double SiCoreEnergy[] = { ToSI::eV(99.2), ToSI::eV(99.8), ToSI::eV(149.7), ToSI::eV(1839.) };
      Si.setCoreEnergy(SiCoreEnergy, 4);

      SelectableElasticSMT SiNISTMott(Si, *NISTMottRS::d_Factory);

      tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\";
      const StringT IIMFPFullPennInterpSiSI = tablePath + "IIMFPFullPennInterpSiSI.csv";
      const StringT interpNUSimReducedDeltaEFullPennSiSI = tablePath + "interpNUSimReducedDeltaEFullPennSiSI.csv";
      const StringT interpNUThetaFullPennSiBGSI = tablePath + "interpNUThetaFullPennSiBGSI.csv";
      const StringT interpSimESE0NUSiBGSI = tablePath + "interpSimESE0NUSiBGSI.csv";
      const char* SiTables[] = {
         IIMFPFullPennInterpSiSI.c_str(),
         interpNUSimReducedDeltaEFullPennSiSI.c_str(),
         interpNUThetaFullPennSiBGSI.c_str(),
         interpSimESE0NUSiBGSI.c_str()
      };
      
      TabulatedInelasticSMT SiDS(Si, 3, SiTables, ToSI::eV(13.54));

      GanachaudMokraniPhononInelasticSMT Siphonon(phononStrength, ToSI::eV(phononE), 300., 11.7, 1.);

      ExpQMBarrierSMT sieqmbsm(&Si);

      MONSEL_MaterialScatterModelT SiMSM(&Si, &sieqmbsm, &sZeroCSD);
      SiMSM.addScatterMechanism(&SiNISTMott);
      SiMSM.addScatterMechanism(&SiDS);
      SiMSM.addScatterMechanism(&Siphonon);

      MONSEL_MaterialScatterModelT SiMSMDeep(&Si, &sieqmbsm, &sZeroCSD);
      SiMSMDeep.addScatterMechanism(&SiNISTMott);
      SiMSMDeep.addScatterMechanism(&SiDS);
      SiMSMDeep.addScatterMechanism(&Siphonon);

      SiMSMDeep.setMinEforTracking(ToSI::eV(50.));

      SphereT sphere(center, MonteCarloSS::ChamberRadius);

      NullMaterialScatterModelT NULL_MSM;

      RegionT chamber(nullptr, &NULL_MSM, &sphere);
      chamber.updateMaterial(*(chamber.getScatterModel()), vacuumMSM);

      NormalMultiPlaneShapeT layer1;
      PlaneT pl1(normalvector, layer1Pos);
      layer1.addPlane(&pl1);
      RegionT layer1Region(&chamber, &ARCMSM, (NormalShapeT*)&layer1);

      NormalMultiPlaneShapeT layer2;
      PlaneT pl2(normalvector, layer2Pos);
      layer2.addPlane(&pl2);
      RegionT layer2Region(&layer1Region, &glCMSM, (NormalShapeT*)&layer2);

      NormalMultiPlaneShapeT layer3;
      PlaneT pl3(normalvector, layer2Pos);
      layer3.addPlane(&pl3);
      RegionT layer3Region(&layer2Region, &SiMSM, (NormalShapeT*)&layer3);

      NormalMultiPlaneShapeT layer4;
      PlaneT pl4(normalvector, layer4Pos);
      RegionT layer4Region(&layer3Region, &SiMSM, (NormalShapeT*)&layer4);

      RegionT deepRegion(&layer3Region, &SiMSMDeep, (NormalShapeT*)&layer4);

      NormalIntersectionShapeT* line = (NormalIntersectionShapeT*)NShapes::createLine(-h, w, linelength, thetal, thetar, radl, radr);
      RegionT lineRegion(&chamber, &PMMAMSM, line);

      GaussianBeamT eg(beamsize, beamE, center);
      const double egCenter[] = { x, y, -h - 20. * 1.e-9 };
      eg.setCenter(egCenter);
      MonteCarloSST monte(&eg, &chamber, nullptr);

      const int nbins = (int)(beamEeV / binSizeEV);
      BackscatterStatsT back(monte, nbins); //printf("48\n");
      monte.addActionListener(back);

      monte.runMultipleTrajectories(nTrajectories);

      const HistogramT& hist = back.backscatterEnergyHistogram(); //printf("49\n");

      const double energyperbineV = beamEeV / hist.binCount();
      const double maxSEbin = 50. / energyperbineV;
      int totalSE = 0;
      for (int j = 0; j < (int)maxSEbin; ++j) {
         totalSE = totalSE + hist.counts(j);
      }

      const double SEf = (float)totalSE / nTrajectories;
      const double bsf = back.backscatterFraction() - SEf;
      printf("%lf %lf %lf %lf %lf\n", beamEeV, xnm, ynm, bsf, SEf);
      monte.removeActionListener(back);

      printf("%d, %d (%d) ended\n", r, c, threadId);
   }

   //__global__ void runCudaSinglePixel(int r, int c)
   //{
   //   int blockId = blockIdx.x + blockIdx.y * gridDim.x;
   //   int threadId = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;
   //   printf("%d, %d (%d) began\n", r, c, threadId);

   //   const double ynm = yvals[r];
   //   const double y = ynm * 1.e-9;
   //   const double xnm = xvals[c];
   //   const double x = xnm * 1.e-9;

   //   GaussianBeamT eg(beamsize, beamE, center);
   //   const double egCenter[] = { x, y, -h - 20. * 1.e-9 };
   //   eg.setCenter(egCenter);
   //   MonteCarloSST monte(&eg, chamber, nullptr);

   //   const int nbins = (int)(beamEeV / binSizeEV);
   //   BackscatterStatsT back(monte, nbins); //printf("48\n");
   //   monte.addActionListener(back);

   //   monte.runMultipleTrajectories(nTrajectories);

   //   const HistogramT& hist = back.backscatterEnergyHistogram(); //printf("49\n");

   //   const double energyperbineV = beamEeV / hist.binCount();
   //   const double maxSEbin = 50. / energyperbineV;
   //   int totalSE = 0;
   //   for (int j = 0; j < (int)maxSEbin; ++j) {
   //      totalSE = totalSE + hist.counts(j);
   //   }

   //   const double SEf = (float)totalSE / nTrajectories;
   //   const double bsf = back.backscatterFraction() - SEf;
   //   printf("%lf %lf %lf %lf %lf\n", beamEeV, xnm, ynm, bsf, SEf);
   //   monte.removeActionListener(back);

   //   printf("%d, %d (%d) ended\n", r, c, threadId);
   //}
//
//   __global__ void runCuda()
//   //void runCuda()
//   {
//      printf("LinesOnLayers: runCuda (%d, %d)\n", yvalsSize, xvalsSize);
//
//      //printf("\n# Trajectories at each landing position: %d", nTrajectories);
//      //printf("\n# Pitch of lines (nm): %.10e", pitchnm);
//      //printf("\n# lines: %d", nlines);
//      //printf("\nLine height (nm): %.10e", hnm);
//      //printf("\nLine bottom width (nm): %.10e", wnm);
//      //printf("\nLine length (nm): %.10e", linelengthnm);
//      //printf("\nLeft and right sidewall angles (deg): %.10e %.10e", thetaldeg, thetardeg);
//      //printf("\nLeft and right top corner radii (nm): %.10e %.10e", radlnm, radrnm);
//      //printf("\nThicknesses of 1st and second layers (nm): %.10e %.10e", layer1thicknessnm, layer2thicknessnm);
//      //printf("\nBeam landing energies (eV): ");
//
//      //for (int i = 0; i < beamEeVvalsLen; i++) {
//      //   printf("\n%.10e", beamEeVvals[i]);
//      //}
//      //printf("\nBeam size (standard deviation, in nm): %.10e", beamsizenm);
//
//      //printf("\n");
//      //printf("\nbeamE (eV)\t x(nm)\t y (nm)\t BSE yield\t SE yield");
//      GaussianBeamT* eg = new GaussianBeamT(beamsize, beamE, center);
//      MonteCarloSST* monte = new MonteCarloSST(eg, chamber, nullptr);
//
//#if (!(defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0)))
//      std::string output;
//      auto start = std::chrono::system_clock::now();
//#endif
//
//      for (int r = 0; r < yvalsSize; ++r) {
//         const double ynm = yvals[r];
//         const double y = ynm * 1.e-9;
//         for (int c = 0; c < xvalsSize; ++c) {
//            double xnm = xvals[c];
//            const double x = xnm * 1.e-9;
//            const double egCenter[] = { x, y, -h - 20. * 1.e-9 };
//            eg->setCenter(egCenter);
//
//            const int nbins = (int)(beamEeV / binSizeEV);
//            BackscatterStatsT* back = new BackscatterStatsT(*monte, nbins); //printf("48\n");
//            monte->addActionListener(*back);
//
//            monte->runMultipleTrajectories(nTrajectories);
//
//            const HistogramT& hist = back->backscatterEnergyHistogram(); //printf("49\n");
//
//            const double energyperbineV = beamEeV / hist.binCount();
//            const double maxSEbin = 50. / energyperbineV;
//            int totalSE = 0;
//            for (int j = 0; j < (int)maxSEbin; ++j) {
//               totalSE = totalSE + hist.counts(j);
//            }
//
//            const double SEf = (float)totalSE / nTrajectories;
//            const double bsf = back->backscatterFraction() - SEf;
//            printf("%lf %lf %lf %lf %lf\n", beamEeV, xnm, ynm, bsf, SEf);
//#if (!(defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0)))
//            output += std::to_string(beamEeV) + " " + std::to_string(xnm) + " " + std::to_string(ynm) + " " + std::to_string(bsf) + " " + std::to_string(SEf) + "\n";
//#endif
//            monte->removeActionListener(*back);
//            delete back;
//         }
//      }
//#if (!(defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0)))
//      auto end = std::chrono::system_clock::now();
//      std::chrono::duration<double> elapsed_seconds = end - start;
//      std::time_t end_time = std::chrono::system_clock::to_time_t(end);
//      std::cout << std::endl << "finished computation at " << std::ctime(&end_time) << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
//      output += std::to_string(elapsed_seconds.count());
//      
//      std::ofstream myfile;
//      myfile.open("output.txt");
//      myfile << output.c_str();
//      myfile.close();
//#endif
//   }
}