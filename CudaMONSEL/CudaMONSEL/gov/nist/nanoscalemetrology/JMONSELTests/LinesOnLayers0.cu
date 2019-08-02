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
#include <time.h>

#include "CudaUtil.h"

namespace LinesOnLayers
{
   __global__ void verifyNUTable1d(const char* fn)
   {
      const NUTableInterpolationT* table = NUTableInterpolation::getInstance(fn);
      const VectorXf& data = table->gettable1d();
      printf("GPU %s\n", fn);
      for (auto v : data) {
         printf("%.5e ", v);
      }
      printf("\n");
   }

   __global__ void verifyNUTable2d(const char* fn, const int r)
   {
      const NUTableInterpolationT* table = NUTableInterpolation::getInstance(fn);
      const MatrixXf& data = table->gettable2d();
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

   void loadNUTable()
   {
       StringT tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\";
       const StringT IIMFPPennInterpglassy = tablePath + "IIMFPPennInterpglassyCSI.csv";
       const StringT SimReducedDeltaEglassy = tablePath + "interpNUSimReducedDeltaEglassyCSI.csv";
       const StringT simTableThetaNUglassy = tablePath + "interpsimTableThetaNUglassyCSI.csv";
       const StringT SimESE0NUglassy = tablePath + "interpSimESE0NUglassyCSI.csv";
       NUTableInterpolation::getInstance(IIMFPPennInterpglassy.c_str());
       NUTableInterpolation::getInstance(SimReducedDeltaEglassy.c_str());
       NUTableInterpolation::getInstance(simTableThetaNUglassy.c_str());
       NUTableInterpolation::getInstance(SimESE0NUglassy.c_str());

       tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\";
       const StringT IIMFPFullPennInterpSiSI = tablePath + "IIMFPFullPennInterpSiSI.csv";
       const StringT interpNUSimReducedDeltaEFullPennSiSI = tablePath + "interpNUSimReducedDeltaEFullPennSiSI.csv";
       const StringT interpNUThetaFullPennSiBGSI = tablePath + "interpNUThetaFullPennSiBGSI.csv";
       const StringT interpSimESE0NUSiBGSI = tablePath + "interpSimESE0NUSiBGSI.csv";
       NUTableInterpolation::getInstance(IIMFPFullPennInterpSiSI.c_str());
       NUTableInterpolation::getInstance(interpNUSimReducedDeltaEFullPennSiSI.c_str());
       NUTableInterpolation::getInstance(interpNUThetaFullPennSiBGSI.c_str());
       NUTableInterpolation::getInstance(interpSimESE0NUSiBGSI.c_str());
   }

   void transferNUTableToCuda()
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
      const VectorXf& data0 = table0->gettable1d();
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
      const MatrixXf& data1 = table1->gettable2d();
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

   __constant__ const float pitchnm = 180.f;
   __constant__ const int nlines = 3;
   __constant__ const float hnm = 120.f;
   __constant__ const float wnm = 80.f;
   __constant__ const float linelengthnm = 1000.f;
   __constant__ const float thetardeg = 3.f;
   __constant__ const float thetaldeg = 3.f;
   __constant__ const float radrnm = 20.f;
   __constant__ const float radlnm = 20.f;
   __constant__ const float layer1thicknessnm = 80.f;
   __constant__ const float layer2thicknessnm = 200.f;

   __constant__ const float beamEeVvals[] = { 500.f };
   __constant__ const int beamEeVvalsLen = 1;
   __constant__ const float beamsizenm = 0.5f;
   __constant__ const float deepnm = 15.f;

   __constant__ const bool trajImg = true;
   __constant__ const int trajImgMaxTraj = 50;
   __constant__ const float trajImgSize = 200e-9f;

   __constant__ const bool VRML = false;
   __constant__ const int VRMLImgMaxTraj = 0;

   //__device__ SEmaterialT* vacuum = nullptr;
   //__device__ ExpQMBarrierSMT* vacuumBarrier = nullptr;
   //__device__ ZeroCSDT* sZeroCSD = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* vacuumMSM = nullptr;

   __constant__ const float PMMAbreakE = 1.60217653e-19 * 45.f;
   __constant__ const float PMMAdensity = 1190.f;
   __constant__ const float PMMAworkfun = 5.5f;
   __constant__ const float PMMAbandgap = 5.f;
   __constant__ const float PMMAEFermi = -5.f;//-PMMAbandgap;
   __constant__ const float PMMApotU = -5.5f - (-5.f);

   //__device__ SEmaterialT* PMMA = nullptr;

   //__device__ SelectableElasticSMT* PMMANISTMott = nullptr;

   //__device__ JoyLuoNieminenCSDT* PMMACSD = nullptr;
   //__device__ FittedInelSMT* PMMAfittedInel = nullptr;
   //__device__ GanachaudMokraniPolaronTrapSMT* PMMApolaron = nullptr;

   //__device__ ExpQMBarrierSMT* pmmaeqmbsm = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* PMMAMSM = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* PMMAMSMDeep = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* ARCMSM = nullptr;

   __constant__ const float glCdensity = 1800.f;
   __constant__ const float glCworkfun = 5.0f;
   __constant__ const float glCbandgap = 0.f;
   __constant__ const float glCEFermi = 20.4f;
   __constant__ const float glCpotU = -5.f - 20.4f;

   //__device__ SEmaterialT* glC = nullptr;

   //__device__ SelectableElasticSMT* glCNISTMott = nullptr;

   //__device__ TabulatedInelasticSMT* glCDS = nullptr;

   //__device__ ExpQMBarrierSMT* glceqmbsm = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* glCMSM = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* glCMSMDeep = nullptr;

   __constant__ const float phononE = 0.063f;
   __constant__ const float phononStrength = 3.f;

   __constant__ const float Sidensity = 2330.f;
   __constant__ const float Siworkfun = 4.85f;
   __constant__ const float Sibandgap = 1.1f;
   __constant__ const float SiEFermi = -1.1f;//-Sibandgap;
   __constant__ const float SipotU = -44.85f - (-1.1f);//-Siworkfun - SiEFermi;

   //__device__ SEmaterialT* Si = nullptr;

   //__device__ SelectableElasticSMT* SiNISTMott = nullptr;

   //__device__ TabulatedInelasticSMT* SiDS = nullptr;

   //__device__ GanachaudMokraniPhononInelasticSMT* Siphonon = nullptr;

   //__device__ ExpQMBarrierSMT* sieqmbsm = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* SiMSM = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* SiMSMDeep = nullptr;

   //__device__ SphereT* sphere = nullptr;

   __constant__ const float pitch = 180.f * 1.e-9f;
   __constant__ const float h = 120.f * 1.e-9f;
   __constant__ const float w = 80.f * 1.e-9f;
   __constant__ const float linelength = 1000.f * 1.e-9f;

   __constant__ const float radperdeg = 3.14159265358979323846f / 180.f;
   __constant__ const float thetar = 3.f * 3.14159265358979323846f / 180.f;
   __constant__ const float thetal = 3.f * 3.14159265358979323846f / 180.f;
   __constant__ const float radr = 20.f * 1.e-9f;
   __constant__ const float radl = 20.f * 1.e-9f;
   __constant__ const float layer1thickness = 80.f * 1.e-9f;
   __constant__ const float layer2thickness = 200.f * 1.e-9f;
   __constant__ const float beamsize = 0.5f * 1.e-9f;
   __constant__ const float deep = 15.f * 1.e-9f;

   __constant__ const double center[] = {
      0.0,
      0.0,
      0.0
   };

   __constant__ const float beamEeV = 500.f;
   __constant__ const float beamE = 1.60217653e-19f * 500.f;
   __constant__ const float binSizeEV = 10.f;

   //__device__ NullMaterialScatterModelT* NULL_MSM = nullptr;

   //__device__ RegionT* chamber = nullptr;

   __constant__ const double normalvector[] = { 0., 0., -1. };
   __constant__ const double layer1Pos[] = { 0., 0., 0. };

   //__device__ NormalMultiPlaneShapeT* layer1 = nullptr;
   //__device__ PlaneT* pl1 = nullptr;
   //__device__ RegionT* layer1Region = nullptr;

   __constant__ const double layer2Pos[] = { 0., 0., 80. * 1.e-9 };
   //__device__ NormalMultiPlaneShapeT* layer2 = nullptr;
   //__device__ PlaneT* pl2 = nullptr;
   //__device__ RegionT* layer2Region = nullptr;

   __constant__ const double layer3Pos[] = { 0., 0., 80. * 1.e-9 + 200. * 1.e-9 };
   //__device__ NormalMultiPlaneShapeT* layer3 = nullptr;
   //__device__ PlaneT* pl3 = nullptr;
   //__device__ RegionT* layer3Region = nullptr;

   __constant__ const double layer4Pos[] = { 0., 0., 80. * 1.e-9 + 200. * 1.e-9 + 15. * 1.e-9 };
   //__device__ NormalMultiPlaneShapeT* layer4 = nullptr;
   //__device__ PlaneT* pl4 = nullptr;
   //__device__ RegionT* layer4Region = nullptr;

   //__device__ RegionT* deepRegion = nullptr;

   __constant__ const float leftmostLineCenterx = -180.f * 1.e-9f * (3.f / 2.f);
   __constant__ const float xcenter = -180.f * 1.e-9f * (3.f / 2.f) + 0.f * 180.f * 1.e-9f;

   //__device__ NormalIntersectionShapeT* line = nullptr;
   //__device__ RegionT* lineRegion = nullptr;

   __device__ float* yvals = nullptr;
   __device__ float* xvals = nullptr;
   __device__ unsigned int yvalsSize = 0;
   __device__ unsigned int xvalsSize = 0;
#else
   const int nTrajectories = 100;

   const float pitchnm = 180.f;
   const int nlines = 3;
   const float hnm = 120.f;
   const float wnm = 80.f;
   const float linelengthnm = 1000.f;
   const float thetardeg = 3.f;
   const float thetaldeg = 3.f;
   const float radrnm = 20.f;
   const float radlnm = 20.f;
   const float layer1thicknessnm = 80.f;
   const float layer2thicknessnm = 200.f;

   const float beamEeVvals[] = { 500.f };
   const int beamEeVvalsLen = 1;
   const float beamsizenm = 0.5f;
   const float deepnm = 15.f;

   const bool trajImg = true;
   const int trajImgMaxTraj = 50;
   const float trajImgSize = 200e-9f;

   const bool VRML = false;
   const int VRMLImgMaxTraj = 0;

   //SEmaterialT* vacuum = nullptr;
   //ExpQMBarrierSMT* vacuumBarrier = nullptr;
   //ZeroCSDT* sZeroCSD = nullptr;

   //MONSEL_MaterialScatterModelT* vacuumMSM = nullptr;

   const float PMMAbreakE = 1.60217653e-19 * 45.f;
   const float PMMAdensity = 1190.f;
   const float PMMAworkfun = 5.5f;
   const float PMMAbandgap = 5.f;
   const float PMMAEFermi = -5.f;//-PMMAbandgap;
   const float PMMApotU = -5.5f - (-5.f);

   //SEmaterialT* PMMA = nullptr;

   //SelectableElasticSMT* PMMANISTMott = nullptr;

   //JoyLuoNieminenCSDT* PMMACSD = nullptr;
   //FittedInelSMT* PMMAfittedInel = nullptr;
   //GanachaudMokraniPolaronTrapSMT* PMMApolaron = nullptr;

   //ExpQMBarrierSMT* pmmaeqmbsm = nullptr;

   //MONSEL_MaterialScatterModelT* PMMAMSM = nullptr;

   //MONSEL_MaterialScatterModelT* PMMAMSMDeep = nullptr;

   //MONSEL_MaterialScatterModelT* ARCMSM = nullptr;

   const float glCdensity = 1800.f;
   const float glCworkfun = 5.0f;
   const float glCbandgap = 0.f;
   const float glCEFermi = 20.4f;
   const float glCpotU = -5.f - 20.4f;

   //SEmaterialT* glC = nullptr;

   //SelectableElasticSMT* glCNISTMott = nullptr;

   //TabulatedInelasticSMT* glCDS = nullptr;

   //ExpQMBarrierSMT* glceqmbsm = nullptr;

   //MONSEL_MaterialScatterModelT* glCMSM = nullptr;

   //MONSEL_MaterialScatterModelT* glCMSMDeep = nullptr;

   const float phononE = 0.063f;
   const float phononStrength = 3.f;

   const float Sidensity = 2330.f;
   const float Siworkfun = 4.85f;
   const float Sibandgap = 1.1f;
   const float SiEFermi = -1.1f;//-Sibandgap;
   const float SipotU = -44.85f - (-1.1f);//-Siworkfun - SiEFermi;

   //SEmaterialT* Si = nullptr;

   //SelectableElasticSMT* SiNISTMott = nullptr;

   //TabulatedInelasticSMT* SiDS = nullptr;

   //GanachaudMokraniPhononInelasticSMT* Siphonon = nullptr;

   //ExpQMBarrierSMT* sieqmbsm = nullptr;

   //MONSEL_MaterialScatterModelT* SiMSM = nullptr;

   //MONSEL_MaterialScatterModelT* SiMSMDeep = nullptr;

   //SphereT* sphere = nullptr;

   const float pitch = 180.f * 1.e-9f;
   const float h = 120.f * 1.e-9f;
   const float w = 80.f * 1.e-9f;
   const float linelength = 1000.f * 1.e-9f;

   const float radperdeg = 3.14159265358979323846f / 180.f;
   const float thetar = 3.f * 3.14159265358979323846f / 180.f;
   const float thetal = 3.f * 3.14159265358979323846f / 180.f;
   const float radr = 20.f * 1.e-9f;
   const float radl = 20.f * 1.e-9f;
   const float layer1thickness = 80.f * 1.e-9f;
   const float layer2thickness = 200.f * 1.e-9f;
   const float beamsize = 0.5f * 1.e-9f;
   const float deep = 15.f * 1.e-9f;

   const double center[] = {
      0.0,
      0.0,
      0.0
   };

   const float beamEeV = 500.f;
   const float beamE = 1.60217653e-19f * 500.f;
   const float binSizeEV = 10.f;

   //NullMaterialScatterModelT* NULL_MSM = nullptr;

   //RegionT* chamber = nullptr;

   const double normalvector[] = { 0., 0., -1. };
   const double layer1Pos[] = { 0., 0., 0. };

   //NormalMultiPlaneShapeT* layer1 = nullptr;
   //PlaneT* pl1 = nullptr;
   //RegionT* layer1Region = nullptr;

   const double layer2Pos[] = { 0., 0., 80. * 1.e-9 };
   //NormalMultiPlaneShapeT* layer2 = nullptr;
   //PlaneT* pl2 = nullptr;
   //RegionT* layer2Region = nullptr;

   const double layer3Pos[] = { 0., 0., 80. * 1.e-9 + 200. * 1.e-9 };
   //NormalMultiPlaneShapeT* layer3 = nullptr;
   //PlaneT* pl3 = nullptr;
   //RegionT* layer3Region = nullptr;

   const double layer4Pos[] = { 0., 0., 80. * 1.e-9 + 200. * 1.e-9 + 15. * 1.e-9 };
   //NormalMultiPlaneShapeT* layer4 = nullptr;
   //PlaneT* pl4 = nullptr;
   //RegionT* layer4Region = nullptr;

   //RegionT* deepRegion = nullptr;

   const float leftmostLineCenterx = -180.f * 1.e-9f * (3.f / 2.f);
   const float xcenter = -180.f * 1.e-9f * (3.f / 2.f) + 0.f * 180.f * 1.e-9f;

   //NormalIntersectionShapeT* line = nullptr;
   //RegionT* lineRegion = nullptr;

   float* yvals = nullptr;
   float* xvals = nullptr;
   unsigned int yvalsSize = 0;
   unsigned int xvalsSize = 0;
#endif

   __host__ __device__ void initRange()
   {
      VectorXf yvalstmp(128);
      for (int i = -64; i < 64; i += 1) {
         yvalstmp.push_back(i);
      }
      //VectorXf yvalstmp(1, -64);

      const float xbottom = wnm / 2.f;
      const float xtop = wnm / 2.f - hnm * ::tanf(thetar);
      const float xstart = xbottom - 100.5f;
      const float xstop = xbottom + 100.5f;
      const float xfinestart = xtop - 20.5f;
      float xfinestop;
      if (thetar < 0.f) xfinestop = xtop + 20.5f;
      else xfinestop = wnm / 2.f + 20.5f;

      VectorXf xvalstmp(128);
      float deltax = 5.f;
      float x = xstart;
      while (x < xfinestart) {
         xvalstmp.push_back(x);
         x += deltax;
      }
      x = xfinestart;
      deltax = 1.f;
      while (x < xfinestop) {
         xvalstmp.push_back(x);
         x += deltax;
      }
      x = xfinestop;
      deltax = 5.f;
      while (x < xstop) {
         xvalstmp.push_back(x);
         x += deltax;
      }
      xvalstmp.push_back(xstop);
      //VectorXf xvalstmp(2);
      //xvalstmp.push_back(xstart);
      //xvalstmp.push_back(xstart + 5.f);

      yvalsSize = yvalstmp.size();
      xvalsSize = xvalstmp.size();

      yvals = new float[yvalsSize];
      xvals = new float[xvalsSize];

      memcpy(yvals, yvalstmp.data(), yvalsSize * sizeof(yvals[0]));
      memcpy(xvals, xvalstmp.data(), xvalsSize * sizeof(xvals[0]));
   }

   __global__ void initCuda()
   {
      printf("LinesOnLayers: initCuda\n");
      for (int i = 0; i < 10; ++i) {
         printf("%.10e\n", Random::random());
      }

      initRange();

      printf("(%d, %d)", xvalsSize, yvalsSize);
   }

   __host__ __device__ void runSinglePixel(const unsigned int r, const unsigned int c, float* result)
   {
#if (!(defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0)))
      //srand(time(NULL));
#endif
      const float ynm = yvals[r];
      const float y = ynm * 1.e-9f;
      const float xnm = xvals[c];
      const float x = xnm * 1.e-9f;

      SEmaterialT vacuum;
      vacuum.setName("SE vacuum");
      ExpQMBarrierSMT vacuumBarrier(&vacuum);
      ZeroCSDT sZeroCSD;

      MONSEL_MaterialScatterModelT vacuumMSM(&vacuum, &vacuumBarrier, &sZeroCSD);

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      const ElementT* componentsCOH[] = { Element::dC, Element::dO, Element::dH };
#else
      const ElementT* componentsCOH[] = { &Element::C, &Element::O, &Element::H };
#endif

      CompositionT PMMAcomp;
      const ::Composition::data_type compositionCOH[] = { 5.f, 2.f, 8.f };
      PMMAcomp.defineByMoleFraction(componentsCOH, 3, compositionCOH, 3);
      SEmaterialT PMMA(PMMAcomp, PMMAdensity);
      PMMA.setName("PMMA");
      PMMA.setWorkfunction(ToSI::eV(PMMAworkfun));
      PMMA.setBandgap(ToSI::eV(PMMAbandgap));
      PMMA.setEnergyCBbottom(ToSI::eV(PMMApotU));

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      SelectableElasticSMT PMMANISTMott(PMMA, *NISTMottRS::d_Factory);
#else
      SelectableElasticSMT PMMANISTMott(PMMA, NISTMottRS::Factory);
#endif

      JoyLuoNieminenCSDT PMMACSD(PMMA, PMMAbreakE);
      FittedInelSMT PMMAfittedInel(PMMA, ToSI::eV(65.4f), PMMACSD);
      GanachaudMokraniPolaronTrapSMT PMMApolaron(2.e7f, 1.f / ToSI::eV(4.f));

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
      PMMAMSMDeep.setMinEforTracking(ToSI::eV(50.f));

      MONSEL_MaterialScatterModelT& ARCMSM = PMMAMSM;

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      const ElementT* glCComponents[] = { Element::dC };
#else
      const ElementT* glCComponents[] = { &Element::C };
#endif

      const ::Composition::data_type glCComposition[] = { 1.f };
      SEmaterialT glC(glCComponents, 1, glCComposition, 1, glCdensity, "glassy Carbon");
      glC.setWorkfunction(ToSI::eV(glCworkfun));
      glC.setEnergyCBbottom(ToSI::eV(glCpotU));
      glC.setBandgap(ToSI::eV(glCbandgap));
      const ::Composition::data_type glCCoreEnergy[] = { ToSI::eV(284.2f) };
      glC.setCoreEnergy(glCCoreEnergy, 1);

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      SelectableElasticSMT glCNISTMott(glC, *NISTMottRS::d_Factory);
#else
      SelectableElasticSMT glCNISTMott(glC, NISTMottRS::Factory);
#endif

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

      glCMSMDeep.setMinEforTracking(ToSI::eV(50.f));

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      const ElementT* SiComponent[] = { Element::dSi };
#else
      const ElementT* SiComponent[] = { &Element::Si };
#endif

      const ::Composition::data_type SiComposition[] = { 1. };
      SEmaterialT Si(SiComponent, 1, SiComposition, 1, Sidensity, "Silicon");
      Si.setWorkfunction(ToSI::eV(Siworkfun));
      Si.setEnergyCBbottom(ToSI::eV(SipotU));
      Si.setBandgap(ToSI::eV(Sibandgap));
      const ::Material::data_type SiCoreEnergy[] = { ToSI::eV(99.2f), ToSI::eV(99.8f), ToSI::eV(149.7f), ToSI::eV(1839.f) };
      Si.setCoreEnergy(SiCoreEnergy, 4);

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      SelectableElasticSMT SiNISTMott(Si, *NISTMottRS::d_Factory);
#else
      SelectableElasticSMT SiNISTMott(Si, NISTMottRS::Factory);
#endif

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

      TabulatedInelasticSMT SiDS(Si, 3, SiTables, ToSI::eV(13.54f));

      GanachaudMokraniPhononInelasticSMT Siphonon(phononStrength, ToSI::eV(phononE), 300.f, 11.7f, 1.f);

      ExpQMBarrierSMT sieqmbsm(&Si);

      MONSEL_MaterialScatterModelT SiMSM(&Si, &sieqmbsm, &sZeroCSD);
      SiMSM.addScatterMechanism(&SiNISTMott);
      SiMSM.addScatterMechanism(&SiDS);
      SiMSM.addScatterMechanism(&Siphonon);

      MONSEL_MaterialScatterModelT SiMSMDeep(&Si, &sieqmbsm, &sZeroCSD);
      SiMSMDeep.addScatterMechanism(&SiNISTMott);
      SiMSMDeep.addScatterMechanism(&SiDS);
      SiMSMDeep.addScatterMechanism(&Siphonon);

      SiMSMDeep.setMinEforTracking(ToSI::eV(50.f));

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
      const double egCenter[] = { x, y, -h - 20.f * 1.e-9f };
      eg.setCenter(egCenter);
      MonteCarloSST monte(&eg, &chamber, nullptr);

      const int nbins = (int)(beamEeV / binSizeEV);
      BackscatterStatsT back(monte, nbins); //printf("48\n");
      monte.addActionListener(back);

      monte.runMultipleTrajectories(nTrajectories);

      const HistogramT& hist = back.backscatterEnergyHistogram(); //printf("49\n");

      const float energyperbineV = beamEeV / hist.binCount();
      const float maxSEbin = 50.f / energyperbineV;
      int totalSE = 0;
      for (int j = 0; j < (int)maxSEbin; ++j) {
         totalSE = totalSE + hist.counts(j);
      }

      const float SEf = (float)totalSE / nTrajectories;
      const float bsf = back.backscatterFraction() - SEf;
      printf("%lf %lf %lf %lf %lf\n", beamEeV, xnm, ynm, bsf, SEf);
      monte.removeActionListener(back);

      result[r * xvalsSize + c] = SEf;
   }

   __global__ void
   //__launch_bounds__(256, 3)
   runCuda(float* result)
   {
      const unsigned int r = blockIdx.y*blockDim.y + threadIdx.y;
      const unsigned int c = blockIdx.x*blockDim.x + threadIdx.x;
      if (r >= yvalsSize || c >= xvalsSize) return;

      const unsigned int blockId = blockIdx.x + blockIdx.y * gridDim.x;
      const unsigned int threadId = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;
      printf("%d, %d (%d) began\n", r, c, threadId);

      runSinglePixel(r, c, result);

      printf("%d, %d (%d) ended\n", r, c, threadId);
   }
}