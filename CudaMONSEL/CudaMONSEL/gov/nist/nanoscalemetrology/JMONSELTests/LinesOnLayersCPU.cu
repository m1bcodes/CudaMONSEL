//#include "gov\nist\nanoscalemetrology\JMONSELTests\LinesOnLayersCPU.cuh"
//
//#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
//#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
//#include "gov\nist\microanalysis\NISTMonte\NullMaterialScatterModel.cuh"
//#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"
//#include "gov\nist\microanalysis\NISTMonte\GaussianBeam.cuh"
//#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
//#include "gov\nist\microanalysis\NISTMonte\Sphere.cuh"
//#include "gov\nist\microanalysis\NISTMonte\BackscatterStats.cuh"
//#include "gov\nist\microanalysis\Utility\Math2.cuh"
//#include "gov\nist\microanalysis\Utility\Histogram.cuh"
//
//#include "gov\nist\nanoscalemetrology\JMONSEL\SEmaterial.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\ExpQMBarrierSM.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\MONSEL_MaterialScatterModel.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\SelectableElasticSM.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\NISTMottRS.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\JoyLuoNieminenCSD.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\FittedInelSM.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\GanachaudMokraniPolaronTrapSM.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\TabulatedInelasticSM.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\GanachaudMokraniPhononInelasticSM.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\NormalMultiPlaneShape.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\NormalIntersectionShape.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\NShapes.cuh"
//
//#include "gov\nist\nanoscalemetrology\JMONSEL\NUTableInterpolation.cuh"
//
//#include "Amphibian\random.cuh"
//
//#include <fstream>
//
//#include <chrono>
//#include <time.h>
//
//#include "CudaUtil.h"
//#include "ImageUtil.h"
//
//namespace LinesOnLayers
//{
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//   __constant__ const int nTrajectories = 100;
//
//   __constant__ const float pitchnm = 180.f;
//   __constant__ const int nlines = 3;
//   __constant__ const float hnm = 120.f;
//   __constant__ const float wnm = 80.f;
//   __constant__ const float linelengthnm = 1000.f;
//   __constant__ const float thetardeg = 3.f;
//   __constant__ const float thetaldeg = 3.f;
//   __constant__ const float radrnm = 20.f;
//   __constant__ const float radlnm = 20.f;
//   __constant__ const float layer1thicknessnm = 80.f;
//   __constant__ const float layer2thicknessnm = 200.f;
//
//   __constant__ const float beamEeVvals[] = { 500.f };
//   __constant__ const int beamEeVvalsLen = 1;
//   __constant__ const float beamsizenm = 0.5f;
//   __constant__ const float deepnm = 15.f;
//
//   __constant__ const bool trajImg = true;
//   __constant__ const int trajImgMaxTraj = 50;
//   __constant__ const float trajImgSize = 200e-9f;
//
//   __constant__ const bool VRML = false;
//   __constant__ const int VRMLImgMaxTraj = 0;
//
//   //__device__ SEmaterialT* vacuum = nullptr;
//   //__device__ ExpQMBarrierSMT* vacuumBarrier = nullptr;
//   //__device__ ZeroCSDT* sZeroCSD = nullptr;
//
//   //__device__ MONSEL_MaterialScatterModelT* vacuumMSM = nullptr;
//
//   __constant__ const float PMMAbreakE = 1.60217653e-19 * 45.f;
//   __constant__ const float PMMAdensity = 1190.f;
//   __constant__ const float PMMAworkfun = 5.5f;
//   __constant__ const float PMMAbandgap = 5.f;
//   __constant__ const float PMMAEFermi = -5.f;//-PMMAbandgap;
//   __constant__ const float PMMApotU = -5.5f - (-5.f);
//
//   //__device__ SEmaterialT* PMMA = nullptr;
//
//   //__device__ SelectableElasticSMT* PMMANISTMott = nullptr;
//
//   //__device__ JoyLuoNieminenCSDT* PMMACSD = nullptr;
//   //__device__ FittedInelSMT* PMMAfittedInel = nullptr;
//   //__device__ GanachaudMokraniPolaronTrapSMT* PMMApolaron = nullptr;
//
//   //__device__ ExpQMBarrierSMT* pmmaeqmbsm = nullptr;
//
//   //__device__ MONSEL_MaterialScatterModelT* PMMAMSM = nullptr;
//
//   //__device__ MONSEL_MaterialScatterModelT* PMMAMSMDeep = nullptr;
//
//   //__device__ MONSEL_MaterialScatterModelT* ARCMSM = nullptr;
//
//   __constant__ const float glCdensity = 1800.f;
//   __constant__ const float glCworkfun = 5.0f;
//   __constant__ const float glCbandgap = 0.f;
//   __constant__ const float glCEFermi = 20.4f;
//   __constant__ const float glCpotU = -5.f - 20.4f;
//
//   //__device__ SEmaterialT* glC = nullptr;
//
//   //__device__ SelectableElasticSMT* glCNISTMott = nullptr;
//
//   //__device__ TabulatedInelasticSMT* glCDS = nullptr;
//
//   //__device__ ExpQMBarrierSMT* glceqmbsm = nullptr;
//
//   //__device__ MONSEL_MaterialScatterModelT* glCMSM = nullptr;
//
//   //__device__ MONSEL_MaterialScatterModelT* glCMSMDeep = nullptr;
//
//   __constant__ const float SiphononE = 0.063f;
//   __constant__ const float SiphononStrength = 3.f;
//   __constant__ const float Sidensity = 2330.f;
//   __constant__ const float Siworkfun = 4.85f;
//   __constant__ const float Sibandgap = 1.1f;
//   __constant__ const float SiEFermi = -1.1f;//-Sibandgap;
//   __constant__ const float SipotU = -4.85f - (-1.1f);//-Siworkfun - SiEFermi;
//   //__device__ SEmaterialT* Si = nullptr;
//
//   //__device__ SelectableElasticSMT* SiNISTMott = nullptr;
//
//   //__device__ TabulatedInelasticSMT* SiDS = nullptr;
//
//   //__device__ GanachaudMokraniPhononInelasticSMT* Siphonon = nullptr;
//
//   //__device__ ExpQMBarrierSMT* sieqmbsm = nullptr;
//
//   //__device__ MONSEL_MaterialScatterModelT* SiMSM = nullptr;
//
//   //__device__ MONSEL_MaterialScatterModelT* SiMSMDeep = nullptr;
//
//   __constant__ const float SiO2density = 2200.f;
//   __constant__ const float SiO2workfun = 10.f;
//   __constant__ const float SiO2phononStrength = 2.f; // Number of phonon modes
//   __constant__ const float SiO2phononE = 0.145f; // Phonon mode energy in eV
//   __constant__ const float SiO2bandgap = 8.9f; // width of band gap in eV
//   __constant__ const float SiO2EFermi = -8.9f; // This puts the Fermi level at the top of the valence band.
//   __constant__ const float SiO2potU = -10.f - (-8.9f);
//
//   //__device__ SphereT* sphere = nullptr;
//
//   __constant__ const float pitch = 180.f * 1.e-9f;
//   __constant__ const float h = 120.f * 1.e-9f;
//   __constant__ const float w = 80.f * 1.e-9f;
//   __constant__ const float linelength = 1000.f * 1.e-9f;
//
//   __constant__ const float radperdeg = 3.14159265358979323846f / 180.f;
//   __constant__ const float thetar = 3.f * 3.14159265358979323846f / 180.f;
//   __constant__ const float thetal = 3.f * 3.14159265358979323846f / 180.f;
//   __constant__ const float radr = 20.f * 1.e-9f;
//   __constant__ const float radl = 20.f * 1.e-9f;
//   __constant__ const float layer1thickness = 80.f * 1.e-9f;
//   __constant__ const float layer2thickness = 200.f * 1.e-9f;
//   __constant__ const float beamsize = 0.5f * 1.e-9f;
//   __constant__ const float deep = 15.f * 1.e-9f;
//
//   __constant__ const double center[] = {
//      0.0,
//      0.0,
//      0.0
//   };
//
//   __constant__ const float beamEeV = 500.f;
//   __constant__ const float beamE = 1.60217653e-19f * 500.f;
//   __constant__ const float binSizeEV = 10.f;
//   __constant__ const float cutoffFractionForSE = 0.1f;
//
//   //__device__ NullMaterialScatterModelT* NULL_MSM = nullptr;
//
//   //__device__ RegionT* chamber = nullptr;
//
//   __constant__ const double normalvector[] = { 0., 0., -1. };
//   __constant__ const double layer1Pos[] = { 0., 0., 0. };
//
//   //__device__ NormalMultiPlaneShapeT* layer1 = nullptr;
//   //__device__ PlaneT* pl1 = nullptr;
//   //__device__ RegionT* layer1Region = nullptr;
//
//   __constant__ const double layer2Pos[] = { 0., 0., 80. * 1.e-9 };
//   //__device__ NormalMultiPlaneShapeT* layer2 = nullptr;
//   //__device__ PlaneT* pl2 = nullptr;
//   //__device__ RegionT* layer2Region = nullptr;
//
//   __constant__ const double layer3Pos[] = { 0., 0., 80. * 1.e-9 + 200. * 1.e-9 };
//   //__device__ NormalMultiPlaneShapeT* layer3 = nullptr;
//   //__device__ PlaneT* pl3 = nullptr;
//   //__device__ RegionT* layer3Region = nullptr;
//
//   __constant__ const double layer4Pos[] = { 0., 0., 80. * 1.e-9 + 200. * 1.e-9 + 15. * 1.e-9 };
//   //__device__ NormalMultiPlaneShapeT* layer4 = nullptr;
//   //__device__ PlaneT* pl4 = nullptr;
//   //__device__ RegionT* layer4Region = nullptr;
//
//   //__device__ RegionT* deepRegion = nullptr;
//
//   __constant__ const float leftmostLineCenterx = -180.f * 1.e-9f * (3.f / 2.f);
//   __constant__ const float xcenter = -180.f * 1.e-9f * (3.f / 2.f) + 0.f * 180.f * 1.e-9f;
//
//   __constant__ const float stripWidth = 30.f * 1e-9f;
//
//   //__device__ NormalIntersectionShapeT* line = nullptr;
//   //__device__ RegionT* lineRegion = nullptr;
//
//   __constant__ const char* glCTables[] = {
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\IIMFPPennInterpglassyCSI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\interpNUSimReducedDeltaEglassyCSI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\interpsimTableThetaNUglassyCSI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\interpSimESE0NUglassyCSI.csv"
//   };
//
//   __constant__ const char* SiTables[] = {
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\IIMFPFullPennInterpSiSI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\interpNUSimReducedDeltaEFullPennSiSI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\interpNUThetaFullPennSiBGSI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\interpSimESE0NUSiBGSI.csv"
//   };
//
//   __constant__ const char* SiO2Tables[] = {
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\IIMFPPennInterpSiO2SI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\interpNUSimReducedDeltaESiO2SI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\interpsimTableThetaNUSiO2SI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\interpSimESE0NUSiO2SI.csv"
//   };
//
//   __device__ const ElementT* COHcomponents[3];
//   __constant__ const Composition::data_type COHcompositions[] = { 5.f, 2.f, 8.f };
//
//   __device__ const ElementT* glCComponents[1];
//   __constant__ const Composition::data_type glCComposition[] = { 1.f };
//
//   __device__ const ElementT* SiComponent[1];
//   __constant__ const Composition::data_type SiComposition[] = { 1. };
//
//   __device__ const ElementT* SiO2Components[2];
//
//   __device__ unsigned int ysize, xsize;
//   __device__ float xstartnm, xstopnm, ystartnm, ystopnm;
//   __device__ const NShapes::LineParams* lineParams[3];
//#else
//   const int nTrajectories = 10;
//
//   const float pitchnm = 180.f;
//   const int nlines = 3;
//   const float hnm = 120.f;
//   const float wnm = 80.f;
//   const float linelengthnm = 120.f;
//   const float thetardeg = 3.f;
//   const float thetaldeg = 3.f;
//   const float radrnm = 20.f;
//   const float radlnm = 20.f;
//   const float layer1thicknessnm = 80.f;
//   const float layer2thicknessnm = 200.f;
//
//   const float beamEeVvals[] = { 500.f };
//   const int beamEeVvalsLen = 1;
//   const float beamsizenm = 0.5f;
//   const float deepnm = 15.f;
//
//   const bool trajImg = true;
//   const int trajImgMaxTraj = 50;
//   const float trajImgSize = 200e-9f;
//
//   const bool VRML = false;
//   const int VRMLImgMaxTraj = 0;
//
//   //SEmaterialT* vacuum = nullptr;
//   //ExpQMBarrierSMT* vacuumBarrier = nullptr;
//   //ZeroCSDT* sZeroCSD = nullptr;
//
//   //MONSEL_MaterialScatterModelT* vacuumMSM = nullptr;
//
//   const float PMMAbreakE = ToSI::eV(45.f);
//   const float PMMAdensity = 1190.f;
//   const float PMMAworkfun = 5.5f;
//   const float PMMAbandgap = 5.f;
//   const float PMMAEFermi = -PMMAbandgap;
//   const float PMMApotU = -PMMAworkfun - PMMAEFermi;
//
//   //SEmaterialT* PMMA = nullptr;
//
//   //SelectableElasticSMT* PMMANISTMott = nullptr;
//
//   //JoyLuoNieminenCSDT* PMMACSD = nullptr;
//   //FittedInelSMT* PMMAfittedInel = nullptr;
//   //GanachaudMokraniPolaronTrapSMT* PMMApolaron = nullptr;
//
//   //ExpQMBarrierSMT* pmmaeqmbsm = nullptr;
//
//   //MONSEL_MaterialScatterModelT* PMMAMSM = nullptr;
//
//   //MONSEL_MaterialScatterModelT* PMMAMSMDeep = nullptr;
//
//   //MONSEL_MaterialScatterModelT* ARCMSM = nullptr;
//
//   const float glCdensity = 1800.f;
//   const float glCworkfun = 5.0f;
//   const float glCbandgap = 0.f;
//   const float glCEFermi = 20.4f;
//   const float glCpotU = -glCworkfun - glCEFermi;
//
//   //SEmaterialT* glC = nullptr;
//
//   //SelectableElasticSMT* glCNISTMott = nullptr;
//
//   //TabulatedInelasticSMT* glCDS = nullptr;
//
//   //ExpQMBarrierSMT* glceqmbsm = nullptr;
//
//   //MONSEL_MaterialScatterModelT* glCMSM = nullptr;
//
//   //MONSEL_MaterialScatterModelT* glCMSMDeep = nullptr;
//
//   const float SiphononE = 0.063f;
//   const float SiphononStrength = 3.f;
//   const float Sidensity = 2330.f;
//   const float Siworkfun = 4.85f;
//   const float Sibandgap = 1.1f;
//   const float SiEFermi = -Sibandgap;
//   const float SipotU = -Siworkfun - SiEFermi;//-Siworkfun - SiEFermi;
//
//   //SEmaterialT* Si = nullptr;
//
//   //SelectableElasticSMT* SiNISTMott = nullptr;
//
//   //TabulatedInelasticSMT* SiDS = nullptr;
//
//   //GanachaudMokraniPhononInelasticSMT* Siphonon = nullptr;
//
//   //ExpQMBarrierSMT* sieqmbsm = nullptr;
//
//   //MONSEL_MaterialScatterModelT* SiMSM = nullptr;
//
//   //MONSEL_MaterialScatterModelT* SiMSMDeep = nullptr;
//
//   const float SiO2density = 2200.f;
//   const float SiO2workfun = 10.f;
//   const float SiO2phononStrength = 2.f; // Number of phonon modes
//   const float SiO2phononE = 0.145f; // Phonon mode energy in eV
//   const float SiO2bandgap = 8.9f; // width of band gap in eV
//   const float SiO2EFermi = -SiO2bandgap; // This puts the Fermi level at the top of the valence band.
//   const float SiO2potU = -SiO2workfun - SiO2EFermi;
//
//   //SphereT* sphere = nullptr;
//
//   const float pitch = pitchnm * 1.e-9f;
//   const float h = hnm * 1.e-9f;
//   const float w = wnm * 1.e-9f;
//   const float linelength = linelengthnm * 1.e-9f;
//
//   const float radperdeg = Math2::PI / 180.f;
//   const float thetar = thetardeg * Math2::PI / 180.f;
//   const float thetal = thetaldeg * Math2::PI / 180.f;
//   const float radr = radrnm * 1.e-9f;
//   const float radl = radlnm * 1.e-9f;
//   const float layer1thickness = layer1thicknessnm * 1.e-9f;
//   const float layer2thickness = layer2thicknessnm * 1.e-9f;
//   const float beamsize = beamsizenm * 1.e-9f;
//   const float deep = deepnm * 1.e-9f;
//
//   const double center[] = {
//      0.0,
//      0.0,
//      0.0
//   };
//
//   const float beamEeV = 500.f;
//   const float beamE = ToSI::eV(beamEeV);
//   const float binSizeEV = 10.f;
//   const float cutoffFractionForSE = 0.1f;
//
//   //NullMaterialScatterModelT* NULL_MSM = nullptr;
//
//   //RegionT* chamber = nullptr;
//
//   const double normalvector[] = { 0., 0., -1. };
//   const double layer1Pos[] = { 0., 0., 0. };
//
//   //NormalMultiPlaneShapeT* layer1 = nullptr;
//   //PlaneT* pl1 = nullptr;
//   //RegionT* layer1Region = nullptr;
//
//   const double layer2Pos[] = { 0., 0., 80. * 1.e-9 };
//   //NormalMultiPlaneShapeT* layer2 = nullptr;
//   //PlaneT* pl2 = nullptr;
//   //RegionT* layer2Region = nullptr;
//
//   const double layer3Pos[] = { 0., 0., 80. * 1.e-9 + 200. * 1.e-9 };
//   //NormalMultiPlaneShapeT* layer3 = nullptr;
//   //PlaneT* pl3 = nullptr;
//   //RegionT* layer3Region = nullptr;
//
//   const double layer4Pos[] = { 0., 0., 80. * 1.e-9 + 200. * 1.e-9 + 15. * 1.e-9 };
//   //NormalMultiPlaneShapeT* layer4 = nullptr;
//   //PlaneT* pl4 = nullptr;
//   //RegionT* layer4Region = nullptr;
//
//   //RegionT* deepRegion = nullptr;
//
//   const float leftmostLineCenterx = -180.f * 1.e-9f * (3.f / 2.f);
//   const float xcenter = -180.f * 1.e-9f * (3.f / 2.f) + 0.f * 180.f * 1.e-9f;
//
//   const float stripWidth = 30.f * 1e-9f;
//
//   //NormalIntersectionShapeT* line = nullptr;
//   //RegionT* lineRegion = nullptr;
//
//   const char* glCTables[] = {
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\IIMFPPennInterpglassyCSI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\interpNUSimReducedDeltaEglassyCSI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\interpsimTableThetaNUglassyCSI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\interpSimESE0NUglassyCSI.csv"
//   };
//
//   const char* SiTables[] = {
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\IIMFPFullPennInterpSiSI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\interpNUSimReducedDeltaEFullPennSiSI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\interpNUThetaFullPennSiBGSI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\interpSimESE0NUSiBGSI.csv"
//   };
//
//   const char* SiO2Tables[] = {
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\IIMFPPennInterpSiO2SI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\interpNUSimReducedDeltaESiO2SI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\interpsimTableThetaNUSiO2SI.csv",
//      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\interpSimESE0NUSiO2SI.csv"
//   };
//
//   const ElementT* COHcomponents[] = { &Element::C, &Element::O, &Element::H };
//   const Composition::data_type COHcompositions[] = { 5.f, 2.f, 8.f };
//
//   const ElementT* glCComponents[] = { &Element::C };
//   const Composition::data_type glCComposition[] = { 1.f };
//
//   const ElementT* SiComponent[] = { &Element::Si };
//   const Composition::data_type SiComposition[] = { 1. };
//
//   const ElementT* SiO2Components[] = { &Element::Si, &Element::O };
//
//   unsigned int ysize, xsize;
//   float xstartnm, xstopnm, ystartnm, ystopnm;
//   const NShapes::LineParams* lineParams[3];
//#endif
//
//   __global__ void verifyNUTable1d(const char* fn)
//   {
//      const NUTableInterpolationT* table = NUTableInterpolation::getInstance(fn);
//      const VectorXf& data = table->gettable1d();
//      printf("GPU %s\n", fn);
//      for (auto v : data) {
//         printf("%.5e ", v);
//      }
//      printf("\n");
//   }
//
//   __global__ void verifyNUTable2d(const char* fn, const int r)
//   {
//      const NUTableInterpolationT* table = NUTableInterpolation::getInstance(fn);
//      const MatrixXf& data = table->gettable2d();
//      printf("GPU %s: row %d\n", fn, r);
//      for (auto v : data[r]) {
//         printf("%.5e ", v);
//      }
//      printf("\n");
//   }
//
//   __global__ void verifyNUTable3d(const char* fn, const int r, const int c)
//   {
//      const NUTableInterpolationT* table = NUTableInterpolation::getInstance(fn);
//      const Matrix3DXf& data = table->gettable3d();
//      printf("GPU %s: row %d, col %d\n", fn, r, c);
//      for (auto v : data[r][c]) {
//         printf("%.5e ", v);
//      }
//      printf("\n");
//   }
//
//   void loadNUTable()
//   {
//      NUTableInterpolation::getInstance(glCTables[0]);
//      NUTableInterpolation::getInstance(glCTables[1]);
//      NUTableInterpolation::getInstance(glCTables[2]);
//      NUTableInterpolation::getInstance(glCTables[3]);
//
//      NUTableInterpolation::getInstance(SiTables[0]);
//      NUTableInterpolation::getInstance(SiTables[1]);
//      NUTableInterpolation::getInstance(SiTables[2]);
//      NUTableInterpolation::getInstance(SiTables[3]);
//
//      NUTableInterpolation::getInstance(SiO2Tables[0]);
//      NUTableInterpolation::getInstance(SiO2Tables[1]);
//      NUTableInterpolation::getInstance(SiO2Tables[2]);
//      NUTableInterpolation::getInstance(SiO2Tables[3]);
//   }
//
//   void transferNUTableToCuda()
//   {
//      NUTableInterpolation::initFactory << <1, 1 >> >();
//      checkCudaErrors(cudaDeviceSynchronize());
//      checkCudaErrors(cudaGetLastError());
//
//      NUTableInterpolation::transferDataToCuda(glCTables[0]);
//      NUTableInterpolation::transferDataToCuda(glCTables[1]);
//      NUTableInterpolation::transferDataToCuda(glCTables[2]);
//      NUTableInterpolation::transferDataToCuda(glCTables[3]);
//
//      char* d_fn = nullptr;
//
//      checkCudaErrors(cudaMalloc((void **)&d_fn, strnlen_s(glCTables[0], 256) * sizeof(char)));
//      checkCudaErrors(cudaMemset(d_fn, 0, strnlen_s(glCTables[0], 256) * sizeof(char)));
//      checkCudaErrors(cudaMemcpy(d_fn, glCTables[0], strnlen_s(glCTables[0], 256) * sizeof(char), cudaMemcpyHostToDevice));
//      verifyNUTable1d << <1, 1 >> >(d_fn);
//      checkCudaErrors(cudaDeviceSynchronize());
//      checkCudaErrors(cudaGetLastError());
//      checkCudaErrors(cudaFree(d_fn));
//      const NUTableInterpolationT* table0 = NUTableInterpolation::getInstance(glCTables[0]);
//      const VectorXf& data0 = table0->gettable1d();
//      printf("CPU %s\n", glCTables[0]);
//      for (auto v : data0) {
//         printf("%.5e ", v);
//      }
//      printf("\n");
//
//      checkCudaErrors(cudaMalloc((void **)&d_fn, strnlen_s(glCTables[1], 256) * sizeof(char)));
//      checkCudaErrors(cudaMemset(d_fn, 0, strnlen_s(glCTables[1], 256) * sizeof(char)));
//      checkCudaErrors(cudaMemcpy(d_fn, glCTables[1], strnlen_s(glCTables[1], 256) * sizeof(char), cudaMemcpyHostToDevice));
//      verifyNUTable2d << <1, 1 >> >(d_fn, 0);
//      checkCudaErrors(cudaDeviceSynchronize());
//      checkCudaErrors(cudaGetLastError());
//      checkCudaErrors(cudaFree(d_fn));
//      const NUTableInterpolationT* table1 = NUTableInterpolation::getInstance(glCTables[1]);
//      const MatrixXf& data1 = table1->gettable2d();
//      printf("CPU %s: row %d\n", glCTables[1], 0);
//      for (auto v : data1[0]) {
//         printf("%.5e ", v);
//      }
//      printf("\n");
//
//      checkCudaErrors(cudaMalloc((void **)&d_fn, strnlen_s(glCTables[2], 256) * sizeof(char)));
//      checkCudaErrors(cudaMemset(d_fn, 0, strnlen_s(glCTables[2], 256) * sizeof(char)));
//      checkCudaErrors(cudaMemcpy(d_fn, glCTables[2], strnlen_s(glCTables[2], 256) * sizeof(char), cudaMemcpyHostToDevice));
//      verifyNUTable3d << <1, 1 >> >(d_fn, 50, 50);
//      checkCudaErrors(cudaDeviceSynchronize());
//      checkCudaErrors(cudaGetLastError());
//      checkCudaErrors(cudaFree(d_fn));
//      const NUTableInterpolationT* table2 = NUTableInterpolation::getInstance(glCTables[2]);
//      const Matrix3DXf& data2 = table2->gettable3d();
//      printf("CPU %s: row %d, col %d\n", glCTables[2], 50, 50);
//      for (auto v : data2[50][50]) {
//         printf("%.5e ", v);
//      }
//      printf("\n");
//
//      NUTableInterpolation::transferDataToCuda(SiTables[0]);
//      NUTableInterpolation::transferDataToCuda(SiTables[1]);
//      NUTableInterpolation::transferDataToCuda(SiTables[2]);
//      NUTableInterpolation::transferDataToCuda(SiTables[3]);
//   }
//
//   __host__ __device__ void initRange()
//   {
//      //VectorXf yvalstmp(128);
//      //for (int i = -64; i < 64; i += 1) {
//      //   yvalstmp.push_back(i);
//      //}
//      ////VectorXf yvalstmp(1, -64);
//
//      //const float xbottom = wnm / 2.f;
//      //const float xtop = wnm / 2.f - hnm * ::tanf(thetar);
//      const float xbottom = 80.f / 2.f;
//      const float xtop = 80.f / 2.f - hnm * ::tanf(thetar);
//      xstartnm = xbottom - 200.f;
//      xstopnm = xbottom + 200.f;
//      //const float xfinestart = xtop - 20.5f;
//      //const float xfinestop = (thetar < 0.f) ? xtop + 20.5f : wnm / 2.f + 20.5f;
//
//      ystartnm = -128.f;
//      ystopnm = 128;
//
//      xsize = 512;
//      ysize = 512;
//
//      //VectorXf xvalstmp(80);
//      //float deltax = 5.f;
//      //float x = xstart;
//      //while (x < xfinestart) {
//      //   xvalstmp.push_back(x);
//      //   x += deltax;
//      //}
//      //x = xfinestart;
//      //deltax = 1.f;
//      //while (x < xfinestop) {
//      //   xvalstmp.push_back(x);
//      //   x += deltax;
//      //}
//      //x = xfinestop;
//      //deltax = 5.f;
//      //while (x < xstop) {
//      //   xvalstmp.push_back(x);
//      //   x += deltax;
//      //}
//      //xvalstmp.push_back(xstop);
//      ////VectorXf xvalstmp(2);
//      ////xvalstmp.push_back(xstart);
//      ////xvalstmp.push_back(xstart + 5.f);
//
//      if (ystopnm - ystartnm < 0) printf("initRange(): ystopnm - ystartnm < 0\n");
//      if (xstopnm - xstartnm < 0) printf("initRange(): xstopnm - xstartnm < 0\n");
//   }
//
//   __global__ void initCuda()
//   {
//      printf("LinesOnLayers: initCuda\n");
//      for (int i = 0; i < 10; ++i) {
//         printf("%.10e\n", Random::random());
//      }
//
//      initRange();
//
//      printf("(%d, %d)", xsize, ysize);
//   }
//
//   __host__ __device__ void runSinglePixel(const unsigned int r, const unsigned int c, float* result)
//   {
//      const float deltay = (ystopnm - ystartnm) / ysize;
//      const float deltax = (xstopnm - xstartnm) / xsize;
//
//      const float ynm = (ystartnm + r * deltay);
//      const float y = ynm * 1.e-9f;
//      const float xnm = (xstartnm + c * deltax);
//      const float x = xnm * 1.e-9f;
//
//      SEmaterialT vacuum; // TODO: move this to global
//      vacuum.setName("SE vacuum");
//      ExpQMBarrierSMT vacuumBarrier(&vacuum); // TODO: move this to global
//      ZeroCSDT sZeroCSD; // TODO: move this to global
//
//      MONSEL_MaterialScatterModelT vacuumMSM(&vacuum, &vacuumBarrier, &sZeroCSD); // TODO: move this to global
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      COHcomponents[0] = Element::dC;
//      COHcomponents[1] = Element::dO;
//      COHcomponents[2] = Element::dH;
//#endif
//
//      CompositionT PMMAcomp;
//      PMMAcomp.defineByMoleFraction(COHcomponents, 3, COHcompositions, 3);
//      SEmaterialT PMMA(PMMAcomp, PMMAdensity); // TODO: move this to global
//      PMMA.setName("PMMA");
//      PMMA.setWorkfunction(ToSI::eV(PMMAworkfun));
//      PMMA.setBandgap(ToSI::eV(PMMAbandgap));
//      PMMA.setEnergyCBbottom(ToSI::eV(PMMApotU));
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      SelectableElasticSMT PMMANISTMott(PMMA, *NISTMottRS::dFactory);
//#else
//      SelectableElasticSMT PMMANISTMott(PMMA, NISTMottRS::Factory);
//#endif
//
//      JoyLuoNieminenCSDT PMMACSD(PMMA, PMMAbreakE);
//      FittedInelSMT PMMAfittedInel(PMMA, ToSI::eV(65.4f), PMMACSD); // TODO: move this to global
//      GanachaudMokraniPolaronTrapSMT PMMApolaron(2.e7f, 1.f / ToSI::eV(4.f)); // TODO: move this to global
//
//      ExpQMBarrierSMT pmmaeqmbsm(&PMMA); // TODO: move this to global
//
//      MONSEL_MaterialScatterModelT PMMAMSM(&PMMA, &pmmaeqmbsm, &sZeroCSD);
//      PMMAMSM.addScatterMechanism(&PMMANISTMott);
//      PMMAMSM.addScatterMechanism(&PMMAfittedInel);
//      PMMAMSM.addScatterMechanism(&PMMApolaron);
//
//      PMMAMSM.setCSD(&PMMACSD);
//
//      MONSEL_MaterialScatterModelT PMMAMSMDeep(&PMMA, &pmmaeqmbsm, &sZeroCSD);
//      PMMAMSMDeep.addScatterMechanism(&PMMANISTMott);
//      PMMAMSMDeep.addScatterMechanism(&PMMAfittedInel);
//      PMMAMSMDeep.addScatterMechanism(&PMMApolaron);
//
//      PMMAMSMDeep.setCSD(&PMMACSD);
//      PMMAMSMDeep.setMinEforTracking(ToSI::eV(cutoffEnergyForSE));
//
//      MONSEL_MaterialScatterModelT& ARCMSM = PMMAMSM;
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      glCComponents[0] = Element::dC;
//#endif
//
//      SEmaterialT glC(glCComponents, 1, glCComposition, 1, glCdensity, "glassy Carbon");
//      glC.setWorkfunction(ToSI::eV(glCworkfun));
//      glC.setEnergyCBbottom(ToSI::eV(glCpotU));
//      glC.setBandgap(ToSI::eV(glCbandgap));
//      const Composition::data_type glCCoreEnergy[] = { ToSI::eV(284.2f) };
//      glC.setCoreEnergy(glCCoreEnergy, 1);
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      SelectableElasticSMT glCNISTMott(glC, *NISTMottRS::dFactory);
//#else
//      SelectableElasticSMT glCNISTMott(glC, NISTMottRS::Factory);
//#endif
//      TabulatedInelasticSMT glCDS(glC, 3, glCTables);
//
//      ExpQMBarrierSMT glceqmbsm(&glC);
//
//      MONSEL_MaterialScatterModelT glCMSM(&glC, &glceqmbsm, &sZeroCSD);
//      glCMSM.addScatterMechanism(&glCNISTMott);
//      glCMSM.addScatterMechanism(&glCDS);
//
//      MONSEL_MaterialScatterModelT glCMSMDeep(&glC, &glceqmbsm, &sZeroCSD);
//      glCMSMDeep.addScatterMechanism(&glCNISTMott);
//      glCMSMDeep.addScatterMechanism(&glCDS);
//
//      glCMSMDeep.setMinEforTracking(ToSI::eV(cutoffEnergyForSE));
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      SiComponent[0] = Element::dSi;
//#endif
//
//      SEmaterialT Si(SiComponent, 1, SiComposition, 1, Sidensity, "Silicon");
//      Si.setWorkfunction(ToSI::eV(Siworkfun));
//      Si.setEnergyCBbottom(ToSI::eV(SipotU));
//      Si.setBandgap(ToSI::eV(Sibandgap));
//      const Material::data_type SiCoreEnergy[] = { ToSI::eV(99.2f), ToSI::eV(99.8f), ToSI::eV(149.7f), ToSI::eV(1839.f) };
//      Si.setCoreEnergy(SiCoreEnergy, 4);
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      SelectableElasticSMT SiNISTMott(Si, *NISTMottRS::dFactory);
//#else
//      SelectableElasticSMT SiNISTMott(Si, NISTMottRS::Factory);
//#endif
//      TabulatedInelasticSMT SiDS(Si, 3, SiTables, ToSI::eV(13.54f));
//
//      GanachaudMokraniPhononInelasticSMT Siphonon(SiphononStrength, ToSI::eV(SiphononE), 300.f, 11.7f, 1.f);
//
//      ExpQMBarrierSMT sieqmbsm(&Si);
//
//      MONSEL_MaterialScatterModelT SiMSM(&Si, &sieqmbsm, &sZeroCSD);
//      SiMSM.addScatterMechanism(&SiNISTMott);
//      SiMSM.addScatterMechanism(&SiDS);
//      SiMSM.addScatterMechanism(&Siphonon);
//
//      MONSEL_MaterialScatterModelT SiMSMDeep(&Si, &sieqmbsm, &sZeroCSD);
//      SiMSMDeep.addScatterMechanism(&SiNISTMott);
//      SiMSMDeep.addScatterMechanism(&SiDS);
//      SiMSMDeep.addScatterMechanism(&Siphonon);
//
//      SiMSMDeep.setMinEforTracking(ToSI::eV(cutoffEnergyForSE));
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      const float SiWeight = Element::dSi->getAtomicWeight();
//      const float OxWeight = 2.f * Element::dO->getAtomicWeight();
//#else
//      const float SiWeight = Element::Si.getAtomicWeight();
//      const float OxWeight = 2.f * Element::O.getAtomicWeight();
//#endif
//      const float totalWeight = SiWeight + OxWeight;
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      SiO2Components[0] = Element::dSi;
//      SiO2Components[1] = Element::dO;
//#endif
//
//      const float SiO2massFrac[] = { SiWeight / totalWeight, OxWeight / totalWeight };
//      SEmaterialT SiO2(SiO2Components, 2, SiO2massFrac, 2, SiO2density, "Silicon Dioxide");
//
//      SiO2.setWorkfunction(ToSI::eV(SiO2workfun));
//      SiO2.setBandgap(ToSI::eV(SiO2bandgap));
//      SiO2.setEnergyCBbottom(ToSI::eV(SiO2potU));
//      const float SiO2ce[] = { ToSI::eV(41.6), ToSI::eV(99.2), ToSI::eV(99.8), ToSI::eV(149.7), ToSI::eV(543.1), ToSI::eV(1839.) };
//      SiO2.setCoreEnergy(SiO2ce, 6);
//
//      // Create scatter mechanisms
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      SelectableElasticSMT SiO2NISTMott((MaterialT&)SiO2, *NISTMottRS::dFactory);
//#else
//      SelectableElasticSMT SiO2NISTMott((MaterialT&)SiO2, NISTMottRS::Factory);
//#endif
//
//      TabulatedInelasticSMT SiO2DS(SiO2, 3, SiO2Tables, ToSI::eV(20. + SiO2bandgap));
//      GanachaudMokraniPhononInelasticSMT SiO2phonon(SiO2phononStrength, ToSI::eV(SiO2phononE), 300., 3.82, 1.);
//      GanachaudMokraniPolaronTrapSMT SiO2polaron(1.0e9, 1. / ToSI::eV(1.));
//      ExpQMBarrierSMT SiO2ClassicalBarrier(&SiO2); // The default.No need to actually execute this line.
//      //SiO2CSD = mon.ZeroCSD(); // The default.No need to actually execute this line.
//      // Make a material scatter model
//      // MSM to be used in thin layer(includes SE generation)
//      MONSEL_MaterialScatterModelT SiO2MSM(&SiO2, &SiO2ClassicalBarrier, &sZeroCSD);
//      SiO2MSM.addScatterMechanism(&SiO2NISTMott);
//      SiO2MSM.addScatterMechanism(&SiO2DS);
//      SiO2MSM.addScatterMechanism(&SiO2phonon);
//      // SiO2MSM.addScatterMechanism(SiO2polaron);
//      //SiO2MSM.setCSD(SiO2CSD);
//      //SiO2MSM.setBarrierSM(SiO2ClassicalBarrier);
//
//      // MSM to be used deep inside(drops electrons with E < 50 eV);
//      //MONSEL_MaterialScatterModelT SiO2MSMDeep(&SiO2, &SiO2ClassicalBarrier, &sZeroCSD);
//      //SiO2MSMDeep.addScatterMechanism(&SiO2NISTMott);
//      //SiO2MSMDeep.addScatterMechanism(&SiO2DS);
//      //SiO2MSMDeep.addScatterMechanism(&SiO2phonon);
//      ////SiO2MSMDeep.setCSD(SiO2CSD);
//      ////SiO2MSMDeep.setBarrierSM(SiO2ClassicalBarrier);
//      //SiO2MSMDeep.setMinEforTracking(ToSI::eV(50.));
//
//      SphereT sphere(center, MonteCarloSS::ChamberRadius);
//
//      NullMaterialScatterModelT NULL_MSM;
//
//      RegionT chamber(nullptr, &NULL_MSM, &sphere);
//      chamber.updateMaterial(*(chamber.getScatterModel()), vacuumMSM);
//
//      NShapes::HorizontalStrip cs0(stripWidth);
//      NormalMultiPlaneShapeT* layer0 = cs0.get();
//      double dist0[3] = { 0., stripWidth / 2, 0. };
//      layer0->translate(dist0);
//      RegionT layer0Region(&chamber, &ARCMSM, (NormalShapeT*)layer0);
//
//      NShapes::HorizontalStrip cs1(stripWidth);
//      NormalMultiPlaneShapeT* layer1 = cs1.get();
//      dist0[1] += stripWidth;
//      layer1->translate(dist0);
//      RegionT layer1Region(&chamber, &SiO2MSM, (NormalShapeT*)layer1);
//
//      NShapes::HorizontalStrip cs2(stripWidth);
//      NormalMultiPlaneShapeT* layer2 = cs2.get();
//      dist0[1] += stripWidth;
//      layer2->translate(dist0);
//      RegionT layer2Region(&chamber, &ARCMSM, (NormalShapeT*)layer2);
//
//      //NormalMultiPlaneShapeT layer1;
//      //PlaneT pl1(normalvector, layer1Pos);
//      //layer1.addPlane(&pl1);
//      //RegionT layer1Region(&chamber, &ARCMSM, (NormalShapeT*)&layer1);
//
//      //NormalMultiPlaneShapeT layer2;
//      //PlaneT pl2(normalvector, layer2Pos);
//      //layer2.addPlane(&pl2);
//      //RegionT layer2Region(&layer1Region, &glCMSM, (NormalShapeT*)&layer2);
//
//      //NormalMultiPlaneShapeT layer3;
//      //PlaneT pl3(normalvector, layer3Pos);
//      //layer3.addPlane(&pl3);
//      //RegionT layer3Region(&layer2Region, &SiMSM, (NormalShapeT*)&layer3);
//
//      //NormalMultiPlaneShapeT layer4;
//      //PlaneT pl4(normalvector, layer4Pos);
//      //RegionT layer4Region(&layer3Region, &SiMSM, (NormalShapeT*)&layer4);
//
//      //RegionT deepRegion(&layer3Region, &SiMSMDeep, (NormalShapeT*)&layer4);
//
//      //NShapes::Line line(-h, w, linelength, thetal, thetar, radl, radr);
//      //const double pivot[3] = { 0.f, 0.f, 0.f };
//      //line.get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
//      //const double dist1[3] = { 0.f, 0.f, linelength / 2. };
//      //line.get()->translate(dist1);
//
//      //RegionT lineRegion(&chamber, &PMMAMSM, (NormalIntersectionShapeT*)line.get());
//
//      NShapes::Line* lines[3];
//      RegionT* regions[3];
//      for (int i = 0; i < nlines; ++i) {
//         //NShapes::Line line(-h, w, linelength, thetal, thetar, radl, radr);
//         lines[i] = new NShapes::Line(-lineParams[i]->h, lineParams[i]->w, lineParams[i]->linelength, lineParams[i]->thetal, lineParams[i]->thetar, lineParams[i]->radl, lineParams[i]->radr);
//         const double pivot[3] = { 0.f, 0.f, 0.f };
//         lines[i]->get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
//         //lines[i]->get()->rotate(pivot, -Math2::PI / 2.f + 20.f / (Math2::PI / 2.f), Math2::PI / 2.f - 20.f / (Math2::PI / 2.f), Math2::PI / 2.f - 20.f / (Math2::PI / 2.f));
//         //const double dist1[3] = { 0.f, 0.f, linelength / 2. };
//         const double offset[3] = { lineParams[i]->x, 0.f, linelength / 2. };
//         lines[i]->get()->translate(offset);
//         lines[i]->addRestrainingPlanes();
//         regions[i] = new RegionT(&chamber, &PMMAMSM, (NormalIntersectionShapeT*)lines[i]->get());
//      }
//
//      //float curx = -w / 2.f;
//      //float curh, curw, curlinelength, curthetal, curthetar, curradl, curradr;
//
//      //NShapes::Line* lines[2];
//      //RegionT* regions[2];
//
//      //for (int i = 0; i < nlines; ++i) {
//      //   curh = getAdjustedValPlusMinus(h, 0.4);
//      //   curw = getAdjustedValPlusMinus(w, 0.4);
//      //   curlinelength = getAdjustedValPlusMinus(linelength, 0.4);
//      //   curthetal = getAdjustedValPlusMinus(thetal, 0.4);
//      //   curthetar = getAdjustedValPlusMinus(thetar, 0.4);
//      //   curradl = getAdjustedValPlusMinus(radl, 0.4);
//      //   curradr = getAdjustedValPlusMinus(radr, 0.4);
//
//      //   //NShapes::Line line(-h, w, linelength, thetal, thetar, radl, radr);
//      //   lines[i] = new NShapes::Line(-curh, curw, curlinelength, curthetal, curthetar, curradl, curradr);
//      //   const double pivot[3] = { 0.f, 0.f, 0.f };
//      //   lines[i]->get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
//      //   //const double dist1[3] = { 0.f, 0.f, linelength / 2. };
//      //   const double dist1[3] = { curx, 0.f, linelength / 2. };
//      //   lines[i]->get()->translate(dist1);
//      //   regions[i] = new RegionT(&chamber, &PMMAMSM, (NormalIntersectionShapeT*)lines[i]->get());
//
//      //   //curx += w * (1. + Random::random());
//      //   curx += w * 2.;
//      //}
//
//      const double egCenter[] = { x, y, -h - 20.f * 1.e-9f };
//      GaussianBeamT eg(beamsize, beamE, egCenter);
//      MonteCarloSST monte(&eg, &chamber, nullptr);
//
//      const int nbins = (int)(beamEeV / binSizeEV);
//      BackscatterStatsT back(monte, nbins);
//      monte.addActionListener(back);
//
//      //#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      //      monte.runMultipleTrajectories(nTrajectories);
//      //#else
//      //      try {
//      //         monte.runMultipleTrajectories(nTrajectories);
//      //      }
//      //      catch (std::exception ex) {
//      //         printf("runSinglePixel: %s\n", ex.what());
//      //      }
//      //#endif
//
//      monte.runMultipleTrajectories(nTrajectories);
//
//      const HistogramT& hist = back.backscatterEnergyHistogram();
//
//      const float energyperbineV = beamEeV / hist.binCount();
//      const float maxSEbin = beamEeV * cutoffFractionForSE / energyperbineV;
//      int totalSE = 0;
//      for (int j = 0; j < (int)maxSEbin; ++j) {
//         totalSE = totalSE + hist.counts(j);
//      }
//
//      const float SEf = (float)totalSE / nTrajectories;
//      const float bsf = back.backscatterFraction() - SEf;
//      //printf("%lf %lf %lf %lf %lf\n", beamEeV, xnm, ynm, bsf, SEf);
//      monte.removeActionListener(back);
//
//      result[r * xsize + c] = SEf;
//
//      for (int i = 0; i < nlines; ++i) {
//         delete lines[i];
//         delete regions[i];
//      }
//   }
//
//   Shapes* s0;
//   Shapes* s1;
//   Shapes* s2;
//   Shapes* s3;
//   Shapes* s4;
//   Shapes* s5;
//   Shapes* s6;
//   Shapes* s7;
//   Shapes* s8;
//   Shapes* s9;
//   Shapes* s10;
//
//   void createShapes()
//   {
//      s0 = new Shapes(stripWidth, stripWidth, stripWidth, lineParams[0], lineParams[1], lineParams[2]);
//      s1 = new Shapes(stripWidth, stripWidth, stripWidth, lineParams[0], lineParams[1], lineParams[2]);
//      s2 = new Shapes(stripWidth, stripWidth, stripWidth, lineParams[0], lineParams[1], lineParams[2]);
//      s3 = new Shapes(stripWidth, stripWidth, stripWidth, lineParams[0], lineParams[1], lineParams[2]);
//      s4 = new Shapes(stripWidth, stripWidth, stripWidth, lineParams[0], lineParams[1], lineParams[2]);
//      s5 = new Shapes(stripWidth, stripWidth, stripWidth, lineParams[0], lineParams[1], lineParams[2]);
//      s6 = new Shapes(stripWidth, stripWidth, stripWidth, lineParams[0], lineParams[1], lineParams[2]);
//      s7 = new Shapes(stripWidth, stripWidth, stripWidth, lineParams[0], lineParams[1], lineParams[2]);
//      s8 = new Shapes(stripWidth, stripWidth, stripWidth, lineParams[0], lineParams[1], lineParams[2]);
//      s9 = new Shapes(stripWidth, stripWidth, stripWidth, lineParams[0], lineParams[1], lineParams[2]);
//      s10 = new Shapes(stripWidth, stripWidth, stripWidth, lineParams[0], lineParams[1], lineParams[2]);
//   }
//
//   void destroyShapes()
//   {
//      delete s0;
//      delete s1;
//      delete s2;
//      delete s3;
//      delete s4;
//      delete s5;
//      delete s6;
//      delete s7;
//      delete s8;
//      delete s9;
//      delete s10;
//   }
//
//   void runSinglePixelCPU(const unsigned int id, const unsigned int r, const unsigned int c, float* result)
//   {
//      Shapes* shapes = nullptr;
//      switch (id) {
//      case 0: shapes = s0; break;
//      case 1: shapes = s1; break;
//      case 2: shapes = s2; break;
//      case 3: shapes = s3; break;
//      case 4: shapes = s4; break;
//      case 5: shapes = s5; break;
//      case 6: shapes = s6; break;
//      case 7: shapes = s7; break;
//      case 8: shapes = s8; break;
//      case 9: shapes = s9; break;
//      case 10: shapes = s10; break;
//      }
//
//      const float deltay = (ystopnm - ystartnm) / ysize;
//      const float deltax = (xstopnm - xstartnm) / xsize;
//
//      const float ynm = (ystartnm + r * deltay);
//      const float y = ynm * 1.e-9f;
//      const float xnm = (xstartnm + c * deltax);
//      const float x = xnm * 1.e-9f;
//
//      SEmaterialT vacuum; // TODO: move this to global
//      vacuum.setName("SE vacuum");
//      ExpQMBarrierSMT vacuumBarrier(&vacuum); // TODO: move this to global
//      ZeroCSDT sZeroCSD; // TODO: move this to global
//
//      MONSEL_MaterialScatterModelT vacuumMSM(&vacuum, &vacuumBarrier, &sZeroCSD); // TODO: move this to global
//
//      CompositionT PMMAcomp;
//      PMMAcomp.defineByMoleFraction(COHcomponents, 3, COHcompositions, 3);
//      SEmaterialT PMMA(PMMAcomp, PMMAdensity); // TODO: move this to global
//      PMMA.setName("PMMA");
//      PMMA.setWorkfunction(ToSI::eV(PMMAworkfun));
//      PMMA.setBandgap(ToSI::eV(PMMAbandgap));
//      PMMA.setEnergyCBbottom(ToSI::eV(PMMApotU));
//
//      SelectableElasticSMT PMMANISTMott(PMMA, NISTMottRS::Factory);
//
//      JoyLuoNieminenCSDT PMMACSD(PMMA, PMMAbreakE);
//      FittedInelSMT PMMAfittedInel(PMMA, ToSI::eV(65.4f), PMMACSD); // TODO: move this to global
//      GanachaudMokraniPolaronTrapSMT PMMApolaron(2.e7f, 1.f / ToSI::eV(4.f)); // TODO: move this to global
//
//      ExpQMBarrierSMT pmmaeqmbsm(&PMMA); // TODO: move this to global
//
//      MONSEL_MaterialScatterModelT PMMAMSM(&PMMA, &pmmaeqmbsm, &sZeroCSD);
//      PMMAMSM.addScatterMechanism(&PMMANISTMott);
//      PMMAMSM.addScatterMechanism(&PMMAfittedInel);
//      PMMAMSM.addScatterMechanism(&PMMApolaron);
//
//      PMMAMSM.setCSD(&PMMACSD);
//
//      MONSEL_MaterialScatterModelT PMMAMSMDeep(&PMMA, &pmmaeqmbsm, &sZeroCSD);
//      PMMAMSMDeep.addScatterMechanism(&PMMANISTMott);
//      PMMAMSMDeep.addScatterMechanism(&PMMAfittedInel);
//      PMMAMSMDeep.addScatterMechanism(&PMMApolaron);
//
//      PMMAMSMDeep.setCSD(&PMMACSD);
//      PMMAMSMDeep.setMinEforTracking(ToSI::eV(cutoffEnergyForSE));
//
//      MONSEL_MaterialScatterModelT& ARCMSM = PMMAMSM;
//
//      SEmaterialT glC(glCComponents, 1, glCComposition, 1, glCdensity, "glassy Carbon");
//      glC.setWorkfunction(ToSI::eV(glCworkfun));
//      glC.setEnergyCBbottom(ToSI::eV(glCpotU));
//      glC.setBandgap(ToSI::eV(glCbandgap));
//      const Composition::data_type glCCoreEnergy[] = { ToSI::eV(284.2f) };
//      glC.setCoreEnergy(glCCoreEnergy, 1);
//
//      SelectableElasticSMT glCNISTMott(glC, NISTMottRS::Factory);
//      TabulatedInelasticSMT glCDS(glC, 3, glCTables);
//
//      ExpQMBarrierSMT glceqmbsm(&glC);
//
//      MONSEL_MaterialScatterModelT glCMSM(&glC, &glceqmbsm, &sZeroCSD);
//      glCMSM.addScatterMechanism(&glCNISTMott);
//      glCMSM.addScatterMechanism(&glCDS);
//
//      MONSEL_MaterialScatterModelT glCMSMDeep(&glC, &glceqmbsm, &sZeroCSD);
//      glCMSMDeep.addScatterMechanism(&glCNISTMott);
//      glCMSMDeep.addScatterMechanism(&glCDS);
//
//      glCMSMDeep.setMinEforTracking(ToSI::eV(cutoffEnergyForSE));
//
//      SEmaterialT Si(SiComponent, 1, SiComposition, 1, Sidensity, "Silicon");
//      Si.setWorkfunction(ToSI::eV(Siworkfun));
//      Si.setEnergyCBbottom(ToSI::eV(SipotU));
//      Si.setBandgap(ToSI::eV(Sibandgap));
//      const Material::data_type SiCoreEnergy[] = { ToSI::eV(99.2f), ToSI::eV(99.8f), ToSI::eV(149.7f), ToSI::eV(1839.f) };
//      Si.setCoreEnergy(SiCoreEnergy, 4);
//
//      SelectableElasticSMT SiNISTMott(Si, NISTMottRS::Factory);
//      TabulatedInelasticSMT SiDS(Si, 3, SiTables, ToSI::eV(13.54f));
//
//      GanachaudMokraniPhononInelasticSMT Siphonon(SiphononStrength, ToSI::eV(SiphononE), 300.f, 11.7f, 1.f);
//
//      ExpQMBarrierSMT sieqmbsm(&Si);
//
//      MONSEL_MaterialScatterModelT SiMSM(&Si, &sieqmbsm, &sZeroCSD);
//      SiMSM.addScatterMechanism(&SiNISTMott);
//      SiMSM.addScatterMechanism(&SiDS);
//      SiMSM.addScatterMechanism(&Siphonon);
//
//      MONSEL_MaterialScatterModelT SiMSMDeep(&Si, &sieqmbsm, &sZeroCSD);
//      SiMSMDeep.addScatterMechanism(&SiNISTMott);
//      SiMSMDeep.addScatterMechanism(&SiDS);
//      SiMSMDeep.addScatterMechanism(&Siphonon);
//
//      SiMSMDeep.setMinEforTracking(ToSI::eV(cutoffEnergyForSE));
//
//      const float SiWeight = Element::Si.getAtomicWeight();
//      const float OxWeight = 2.f * Element::O.getAtomicWeight();
//      const float totalWeight = SiWeight + OxWeight;
//
//      const float SiO2massFrac[] = { SiWeight / totalWeight, OxWeight / totalWeight };
//      SEmaterialT SiO2(SiO2Components, 2, SiO2massFrac, 2, SiO2density, "Silicon Dioxide");
//
//      SiO2.setWorkfunction(ToSI::eV(SiO2workfun));
//      SiO2.setBandgap(ToSI::eV(SiO2bandgap));
//      SiO2.setEnergyCBbottom(ToSI::eV(SiO2potU));
//      const float SiO2ce[] = { ToSI::eV(41.6), ToSI::eV(99.2), ToSI::eV(99.8), ToSI::eV(149.7), ToSI::eV(543.1), ToSI::eV(1839.) };
//      SiO2.setCoreEnergy(SiO2ce, 6);
//
//      // Create scatter mechanisms
//      SelectableElasticSMT SiO2NISTMott((MaterialT&)SiO2, NISTMottRS::Factory);
//
//      TabulatedInelasticSMT SiO2DS(SiO2, 3, SiO2Tables, ToSI::eV(20. + SiO2bandgap));
//      GanachaudMokraniPhononInelasticSMT SiO2phonon(SiO2phononStrength, ToSI::eV(SiO2phononE), 300., 3.82, 1.);
//      GanachaudMokraniPolaronTrapSMT SiO2polaron(1.0e9, 1. / ToSI::eV(1.));
//      ExpQMBarrierSMT SiO2ClassicalBarrier(&SiO2); // The default.No need to actually execute this line.
//      //SiO2CSD = mon.ZeroCSD(); // The default.No need to actually execute this line.
//      // Make a material scatter model
//      // MSM to be used in thin layer(includes SE generation)
//      MONSEL_MaterialScatterModelT SiO2MSM(&SiO2, &SiO2ClassicalBarrier, &sZeroCSD);
//      SiO2MSM.addScatterMechanism(&SiO2NISTMott);
//      SiO2MSM.addScatterMechanism(&SiO2DS);
//      SiO2MSM.addScatterMechanism(&SiO2phonon);
//      // SiO2MSM.addScatterMechanism(SiO2polaron);
//      //SiO2MSM.setCSD(SiO2CSD);
//      //SiO2MSM.setBarrierSM(SiO2ClassicalBarrier);
//
//      // MSM to be used deep inside(drops electrons with E < 50 eV);
//      //MONSEL_MaterialScatterModelT SiO2MSMDeep(&SiO2, &SiO2ClassicalBarrier, &sZeroCSD);
//      //SiO2MSMDeep.addScatterMechanism(&SiO2NISTMott);
//      //SiO2MSMDeep.addScatterMechanism(&SiO2DS);
//      //SiO2MSMDeep.addScatterMechanism(&SiO2phonon);
//      ////SiO2MSMDeep.setCSD(SiO2CSD);
//      ////SiO2MSMDeep.setBarrierSM(SiO2ClassicalBarrier);
//      //SiO2MSMDeep.setMinEforTracking(ToSI::eV(50.));
//
//      SphereT sphere(center, MonteCarloSS::ChamberRadius);
//
//      NullMaterialScatterModelT NULL_MSM;
//
//      RegionT chamber(nullptr, &NULL_MSM, &sphere);
//      chamber.updateMaterial(*(chamber.getScatterModel()), vacuumMSM);
//
//      //NShapes::HorizontalStrip cs0(stripWidth);
//      //NormalMultiPlaneShapeT* layer0 = cs0.get();
//      //double dist0[3] = { 0., stripWidth / 2, 0. };
//      //layer0->translate(dist0);
//      //RegionT layer0Region(&chamber, &ARCMSM, (NormalShapeT*)layer0);
//      RegionT layer0Region(&chamber, &ARCMSM, (NormalShapeT*)shapes->hs0.get());
//
//      //NShapes::HorizontalStrip cs1(stripWidth);
//      //NormalMultiPlaneShapeT* layer1 = cs1.get();
//      //dist0[1] += stripWidth;
//      //layer1->translate(dist0);
//      //RegionT layer1Region(&chamber, &SiO2MSM, (NormalShapeT*)layer1);
//      RegionT layer1Region(&chamber, &ARCMSM, (NormalShapeT*)shapes->hs1.get());
//
//      //NShapes::HorizontalStrip cs2(stripWidth);
//      //NormalMultiPlaneShapeT* layer2 = cs2.get();
//      //dist0[1] += stripWidth;
//      //layer2->translate(dist0);
//      //RegionT layer2Region(&chamber, &ARCMSM, (NormalShapeT*)layer2);
//      RegionT layer2Region(&chamber, &ARCMSM, (NormalShapeT*)shapes->hs2.get());
//
//      //NormalMultiPlaneShapeT layer1;
//      //PlaneT pl1(normalvector, layer1Pos);
//      //layer1.addPlane(&pl1);
//      //RegionT layer1Region(&chamber, &ARCMSM, (NormalShapeT*)&layer1);
//
//      //NormalMultiPlaneShapeT layer2;
//      //PlaneT pl2(normalvector, layer2Pos);
//      //layer2.addPlane(&pl2);
//      //RegionT layer2Region(&layer1Region, &glCMSM, (NormalShapeT*)&layer2);
//
//      //NormalMultiPlaneShapeT layer3;
//      //PlaneT pl3(normalvector, layer3Pos);
//      //layer3.addPlane(&pl3);
//      //RegionT layer3Region(&layer2Region, &SiMSM, (NormalShapeT*)&layer3);
//
//      //NormalMultiPlaneShapeT layer4;
//      //PlaneT pl4(normalvector, layer4Pos);
//      //RegionT layer4Region(&layer3Region, &SiMSM, (NormalShapeT*)&layer4);
//
//      //RegionT deepRegion(&layer3Region, &SiMSMDeep, (NormalShapeT*)&layer4);
//
//      //NShapes::Line line(-h, w, linelength, thetal, thetar, radl, radr);
//      //const double pivot[3] = { 0.f, 0.f, 0.f };
//      //line.get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
//      //const double dist1[3] = { 0.f, 0.f, linelength / 2. };
//      //line.get()->translate(dist1);
//
//      //RegionT lineRegion(&chamber, &PMMAMSM, (NormalIntersectionShapeT*)line.get());
//
//      //NShapes::Line* lines[3];
//      RegionT* regions[3];
//      //for (int i = 0; i < nlines; ++i) {
//      //   //NShapes::Line line(-h, w, linelength, thetal, thetar, radl, radr);
//      //   lines[i] = new NShapes::Line(-lineParams[i]->h, lineParams[i]->w, lineParams[i]->linelength, lineParams[i]->thetal, lineParams[i]->thetar, lineParams[i]->radl, lineParams[i]->radr);
//      //   const double pivot[3] = { 0.f, 0.f, 0.f };
//      //   lines[i]->get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
//      //   //lines[i]->get()->rotate(pivot, -Math2::PI / 2.f + 20.f / (Math2::PI / 2.f), Math2::PI / 2.f - 20.f / (Math2::PI / 2.f), Math2::PI / 2.f - 20.f / (Math2::PI / 2.f));
//      //   //const double dist1[3] = { 0.f, 0.f, linelength / 2. };
//      //   const double offset[3] = { lineParams[i]->x, 0.f, linelength / 2. };
//      //   lines[i]->get()->translate(offset);
//      //   lines[i]->addRestrainingPlanes();
//      //   regions[i] = new RegionT(&chamber, &PMMAMSM, (NormalIntersectionShapeT*)lines[i]->get());
//      //}
//      regions[0] = new RegionT(&chamber, &PMMAMSM, (NormalIntersectionShapeT*)shapes->l0.get());
//      regions[1] = new RegionT(&chamber, &PMMAMSM, (NormalIntersectionShapeT*)shapes->l1.get());
//      regions[2] = new RegionT(&chamber, &PMMAMSM, (NormalIntersectionShapeT*)shapes->l2.get());
//
//      //float curx = -w / 2.f;
//      //float curh, curw, curlinelength, curthetal, curthetar, curradl, curradr;
//
//      //NShapes::Line* lines[2];
//      //RegionT* regions[2];
//
//      //for (int i = 0; i < nlines; ++i) {
//      //   curh = getAdjustedValPlusMinus(h, 0.4);
//      //   curw = getAdjustedValPlusMinus(w, 0.4);
//      //   curlinelength = getAdjustedValPlusMinus(linelength, 0.4);
//      //   curthetal = getAdjustedValPlusMinus(thetal, 0.4);
//      //   curthetar = getAdjustedValPlusMinus(thetar, 0.4);
//      //   curradl = getAdjustedValPlusMinus(radl, 0.4);
//      //   curradr = getAdjustedValPlusMinus(radr, 0.4);
//
//      //   //NShapes::Line line(-h, w, linelength, thetal, thetar, radl, radr);
//      //   lines[i] = new NShapes::Line(-curh, curw, curlinelength, curthetal, curthetar, curradl, curradr);
//      //   const double pivot[3] = { 0.f, 0.f, 0.f };
//      //   lines[i]->get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
//      //   //const double dist1[3] = { 0.f, 0.f, linelength / 2. };
//      //   const double dist1[3] = { curx, 0.f, linelength / 2. };
//      //   lines[i]->get()->translate(dist1);
//      //   regions[i] = new RegionT(&chamber, &PMMAMSM, (NormalIntersectionShapeT*)lines[i]->get());
//
//      //   //curx += w * (1. + Random::random());
//      //   curx += w * 2.;
//      //}
//
//      const double egCenter[] = { x, y, -h - 20.f * 1.e-9f };
//      GaussianBeamT eg(beamsize, beamE, egCenter);
//      MonteCarloSST monte(&eg, &chamber, nullptr);
//
//      const int nbins = (int)(beamEeV / binSizeEV);
//      BackscatterStatsT back(monte, nbins);
//      monte.addActionListener(back);
//
//      //#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      //      monte.runMultipleTrajectories(nTrajectories);
//      //#else
//      //      try {
//      //         monte.runMultipleTrajectories(nTrajectories);
//      //      }
//      //      catch (std::exception ex) {
//      //         printf("runSinglePixel: %s\n", ex.what());
//      //      }
//      //#endif
//
//      monte.runMultipleTrajectories(nTrajectories);
//
//      const HistogramT& hist = back.backscatterEnergyHistogram();
//
//      const float energyperbineV = beamEeV / hist.binCount();
//      const float maxSEbin = beamEeV * cutoffFractionForSE / energyperbineV;
//      int totalSE = 0;
//      for (int j = 0; j < (int)maxSEbin; ++j) {
//         totalSE = totalSE + hist.counts(j);
//      }
//
//      const float SEf = (float)totalSE / nTrajectories;
//      const float bsf = back.backscatterFraction() - SEf;
//      //printf("%lf %lf %lf %lf %lf\n", beamEeV, xnm, ynm, bsf, SEf);
//      monte.removeActionListener(back);
//
//      result[r * xsize + c] = SEf;
//
//      for (int i = 0; i < nlines; ++i) {
//         //delete lines[i];
//         delete regions[i];
//      }
//   }
//
//   __global__ void
//      //__launch_bounds__(256, 3)
//      runCuda(float* result)
//   {
//      const unsigned int r = blockIdx.y*blockDim.y + threadIdx.y;
//      const unsigned int c = blockIdx.x*blockDim.x + threadIdx.x;
//      if (r >= ysize || c >= xsize) return;
//
//      const unsigned int blockId = blockIdx.x + blockIdx.y * gridDim.x;
//      const unsigned int threadId = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;
//      printf("%d, %d (%d) began\n", r, c, threadId);
//
//      runSinglePixel(r, c, result);
//
//      printf("%d, %d (%d) ended\n", r, c, threadId);
//   }
//
//   void runSinglePixelThread(int id, const unsigned int r, const unsigned int c, float* result)
//   {
//      try {
//         //runSinglePixel(r, c, result);
//         runSinglePixelCPU(id, r, c, result);
//      }
//      catch (std::exception ex) {
//         printf("%s\n", ex.what());
//      }
//   }
//
//   inline float getAdjustedValPlusMinus(float val, float fraction)
//   {
//      return val * (1.f + (1.f - Random::random() * 2.f) * fraction); // val +/- val*fraction
//   }
//
//   void lineProjection()
//   {
//      NShapes::HorizontalStrip hs0(stripWidth);
//      NormalMultiPlaneShapeT* layer0 = hs0.get();
//      double offset[3] = { 0., stripWidth / 2, 0. };
//      layer0->translate(offset);
//
//      NShapes::HorizontalStrip hs1(stripWidth);
//      NormalMultiPlaneShapeT* layer1 = hs1.get();
//      offset[1] += stripWidth;
//      layer1->translate(offset);
//
//      NShapes::HorizontalStrip hs2(stripWidth);
//      NormalMultiPlaneShapeT* layer2 = hs2.get();
//      offset[1] += stripWidth;
//      layer2->translate(offset);
//
//      hs0.calcGroundtruth();
//      hs1.calcGroundtruth();
//      hs2.calcGroundtruth();
//
//      const double p[3] = { xstartnm * 1.e-9, ystartnm * 1.e-9, 0. };
//      const double n[3] = { 0., 0., -1. };
//      PlaneT projectionPlane(n, p); // projection plane
//      const double axis0[3] = { 1., 0., 0. }; // absolute (projection) vector from plane origin
//      const double axis1[3] = { 0., 1., 0. }; // absolute (projection) vector from plane origin
//
//      const float xlenperpix = (xstopnm - xstartnm) / xsize * 1.e-9;
//      const float ylenperpix = (ystopnm - ystartnm) / ysize * 1.e-9;
//
//      char* gt = new char[ysize * xsize];
//      memset(gt, 0, sizeof(gt[0]) * ysize * xsize);
//
//      hs0.calcRasterization(projectionPlane, axis0, axis1, xlenperpix, ylenperpix, gt, xsize, ysize);
//      hs1.calcRasterization(projectionPlane, axis0, axis1, xlenperpix, ylenperpix, gt, xsize, ysize);
//      hs2.calcRasterization(projectionPlane, axis0, axis1, xlenperpix, ylenperpix, gt, xsize, ysize);
//
//      float curx = -w;
//
//      for (int i = 0; i < nlines; ++i) {
//         lineParams[i] = new NShapes::LineParams(
//            //NShapes::LineParams lp(
//            getAdjustedValPlusMinus(h, 0.15f),
//            w,
//            getAdjustedValPlusMinus(linelength, 0.f),
//            getAdjustedValPlusMinus(thetal, 0.f),
//            getAdjustedValPlusMinus(thetar, 0.2f),
//            getAdjustedValPlusMinus(radl, 0.f),
//            getAdjustedValPlusMinus(radr, 0.f),
//            curx
//            );
//         //NShapes::Line line(-h, w, linelength, thetal, thetar, radl, radr);
//         NShapes::Line line(-lineParams[i]->h, lineParams[i]->w, lineParams[i]->linelength, lineParams[i]->thetal, lineParams[i]->thetar, lineParams[i]->radl, lineParams[i]->radr);
//         //NShapes::Line line(-lp.h, lp.w, lp.linelength, lp.thetal, lp.thetar, lp.radl, lp.radr);
//         const double pivot[3] = { 0.f, 0.f, 0.f };
//         line.get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
//         //line.get()->rotate(pivot, -Math2::PI / 2.f + 20.f / (Math2::PI / 2.f), Math2::PI / 2.f - 20.f / (Math2::PI / 2.f), Math2::PI / 2.f - 20.f / (Math2::PI / 2.f));
//         //const double dist1[3] = { 0.f, 0.f, linelength / 2. };
//         //const double dist1[3] = { lp.x, 0.f, linelength / 2. };
//         const double offset[3] = { lineParams[i]->x, 0.f, linelength / 2. };
//         line.get()->translate(offset);
//         line.calcGroundtruth(); // get points/line segments that need to be projected
//         line.calcRasterization(projectionPlane, axis0, axis1, xlenperpix, ylenperpix, gt, xsize, ysize); // needs to be last since bottom needs to be removed due to same material
//
//         curx += w * 1.5;
//      }
//
//#if (!(defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0)))
//      ImageUtil::saveImage("gt.bmp", gt, xsize, ysize);
//#endif
//
//      delete[] gt;
//
//      //NShapes::TestProjection();
//   }
//
//   // note must call initRange() first
//   __host__ __device__ Shapes::Shapes(const float stripWidth0, const float stripWidth1, const float stripWidth2, const NShapes::LineParams* lp0, const NShapes::LineParams* lp1, const NShapes::LineParams* lp2) :
//      hs0(stripWidth0),
//      hs1(stripWidth1),
//      hs2(stripWidth2),
//      l0(-lp0->h, lp0->w, lp0->linelength, lp0->thetal, lp0->thetar, lp0->radl, lp0->radr),
//      l1(-lp1->h, lp1->w, lp1->linelength, lp1->thetal, lp1->thetar, lp1->radl, lp1->radr),
//      l2(-lp2->h, lp2->w, lp2->linelength, lp2->thetal, lp2->thetar, lp2->radl, lp2->radr)
//   {
//      const double stripOffset0[3] = { 0., stripWidth / 2, 0. };
//      hs0.get()->translate(stripOffset0);
//
//      const double stripOffset1[3] = { 0., stripWidth * 3. / 2., 0. };
//      hs1.get()->translate(stripOffset1);
//
//      const double stripOffset2[3] = { 0., stripWidth * 5. / 2., 0. };
//      hs2.get()->translate(stripOffset2);
//
//      const double pivot[3] = { 0.f, 0.f, 0.f };
//      const float startx = -w;
//      const float separation = w * 1.5f;
//
//      l0.get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
//      const double lineOffset0[3] = { startx + separation * 0.0f, 0.f, linelength / 2. };
//      l0.get()->translate(lineOffset0);
//
//      l1.get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
//      const double lineOffset1[3] = { startx + separation * 1.f, 0.f, linelength / 2. };
//      l1.get()->translate(lineOffset1);
//
//      l2.get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
//      const double lineOffset2[3] = { startx + separation * 2.f, 0.f, linelength / 2. };
//      l2.get()->translate(lineOffset2);
//   }
//}