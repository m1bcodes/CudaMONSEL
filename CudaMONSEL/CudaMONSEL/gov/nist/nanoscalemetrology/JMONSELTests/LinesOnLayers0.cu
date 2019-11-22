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
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalDifferenceShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NShapes.cuh"

#include "gov\nist\nanoscalemetrology\JMONSEL\NUTableInterpolation.cuh"

#include "Amphibian\random.cuh"

#include <fstream>
#include <string>

#include <chrono>
#include <time.h>

#include "CudaUtil.h"
#include "ImageUtil.h"

namespace LinesOnLayers
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __device__ unsigned int nTrajectories = 100;

   __constant__ const float pitchnm = 180.f;
   __device__ unsigned int nlines = 3;
   __device__ unsigned int linemat = 0;
   __constant__ const float hnm = 120.f;
   __constant__ const float wnm = 80.f;
   __device__ float linelengthnm = 1000.f;
   __constant__ const float thetardeg = 3.f;
   __constant__ const float thetaldeg = 3.f;
   __constant__ const float radrnm = 20.f;
   __constant__ const float radlnm = 20.f;
   __constant__ const float layerthicknessnm = 30.f;
   __constant__ const float layer1thicknessnm = 80.f;
   __constant__ const float layer2thicknessnm = 200.f;

   __constant__ const float beamEeVvals[] = { 500.f };
   __constant__ const int beamEeVvalsLen = 1;
   __device__ float beamsizenm = 0.5f;
   __device__ float beamznm = -120.f - 20.f;
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

   __constant__ const float SiphononE = 0.063f;
   __constant__ const float SiphononStrength = 3.f;
   __constant__ const float Sidensity = 2330.f;
   __constant__ const float Siworkfun = 4.85f;
   __constant__ const float Sibandgap = 1.1f;
   __constant__ const float SiEFermi = -1.1f;//-Sibandgap;
   __constant__ const float SipotU = -4.85f - (-1.1f);//-Siworkfun - SiEFermi;
   //__device__ SEmaterialT* Si = nullptr;

   //__device__ SelectableElasticSMT* SiNISTMott = nullptr;

   //__device__ TabulatedInelasticSMT* SiDS = nullptr;

   //__device__ GanachaudMokraniPhononInelasticSMT* Siphonon = nullptr;

   //__device__ ExpQMBarrierSMT* sieqmbsm = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* SiMSM = nullptr;

   //__device__ MONSEL_MaterialScatterModelT* SiMSMDeep = nullptr;

   __constant__ const float SiO2density = 2200.f;
   __constant__ const float SiO2workfun = 10.f;
   __constant__ const float SiO2phononStrength = 2.f; // Number of phonon modes
   __constant__ const float SiO2phononE = 0.145f; // Phonon mode energy in eV
   __constant__ const float SiO2bandgap = 8.9f; // width of band gap in eV
   __constant__ const float SiO2EFermi = -8.9f; // This puts the Fermi level at the top of the valence band.
   __constant__ const float SiO2potU = -10.f - (-8.9f);

   //__device__ SphereT* sphere = nullptr;

   __constant__ const float pitch = 180.f * 1.e-9f;
   __constant__ const float h = 120.f * 1.e-9f;
   __constant__ const float w = 80.f * 1.e-9f;
   __device__ float linelength = 1000.f * 1.e-9f;

   __constant__ const float radperdeg = 3.14159265358979323846f / 180.f;
   __constant__ const float thetar = 3.f * 3.14159265358979323846f / 180.f;
   __constant__ const float thetal = 3.f * 3.14159265358979323846f / 180.f;
   __constant__ const float radr = 20.f * 1.e-9f;
   __constant__ const float radl = 20.f * 1.e-9f;
   __device__ float layerthickness = 80.f * 1.e-9f;
   __constant__ const float layer1thickness = 80.f * 1.e-9f;
   __constant__ const float layer2thickness = 200.f * 1.e-9f;
   __device__ float beamsize = 0.5f * 1.e-9f;
   __device__ float beamz = -140.f * 1.e-9f;
   __constant__ const float deep = 15.f * 1.e-9f;

   __constant__ const double center[] = {
      0.0,
      0.0,
      0.0
   };

   __device__ float beamEeV = 500.f;
   __device__ float beamE = 1.60217653e-19f * 500.f;
   __constant__ const float binSizeEV = 10.f;
   __constant__ const float cutoffEnergyForSE = 50.f;
   __device__ float beamphideg = 1.f;
   __device__ float beamthetadeg = 0.f;
   __device__ float beamfocallength = 20.f;;

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

   __constant__ const char* glCTables[] = {
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\IIMFPPennInterpglassyCSI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\interpNUSimReducedDeltaEglassyCSI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\interpsimTableThetaNUglassyCSI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\interpSimESE0NUglassyCSI.csv"
      };

   __constant__ const char* SiTables[] = {
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\IIMFPFullPennInterpSiSI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\interpNUSimReducedDeltaEFullPennSiSI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\interpNUThetaFullPennSiBGSI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\interpSimESE0NUSiBGSI.csv"
   };

   __constant__ const char* SiO2Tables[] = {
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\IIMFPPennInterpSiO2SI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\interpNUSimReducedDeltaESiO2SI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\interpsimTableThetaNUSiO2SI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\interpSimESE0NUSiO2SI.csv"
      };

   __device__ const ElementT* COHcomponents[3];
   __constant__ const Composition::data_type COHcompositions[] = { 5.f, 2.f, 8.f };

   __device__ const ElementT* glCComponents[1];
   __constant__ const Composition::data_type glCComposition[] = { 1.f };

   __device__ const ElementT* SiComponent[1];
   __constant__ const Composition::data_type SiComposition[] = { 1. };

   __device__ const ElementT* SiO2Components[2];

   __device__ unsigned int ysize, xsize;
   __device__ float xstartnm, xstopnm, ystartnm, ystopnm;

   __device__ NShapes::LineParams** lineParams;

   __device__ unsigned int nhstrips;
   __device__ NShapes::HorizontalStripParams** hstripParams;
#else
   unsigned int nTrajectories = 100;
   //unsigned int nTrajectories = 250;

   //const float wnm = 80.f;
   //const float thetardeg = 3.f;
   //const float thetaldeg = 3.f;
   //const float radrnm = 20.f;
   //const float radlnm = 20.f;

   const float pitchnm = 180.f;
   unsigned int nlines = 3;
   unsigned int linemat = 0;
   const float hnm = 200.f;
   const float wnm = 40.f;
   float linelengthnm = 120.f;
   const float thetardeg = 1.f;
   const float thetaldeg = 1.f;
   const float radrnm = wnm / 3.f;
   const float radlnm = wnm / 3.f;
   const float layerthicknessnm = 20.f;
   const float layer1thicknessnm = 80.f;
   const float layer2thicknessnm = 200.f;

   const float beamEeVvals[] = { 500.f };
   const int beamEeVvalsLen = 1;
   float beamsizenm = 0.1f;
   float beamznm = -hnm - 20.f;
   //float beamsizenm = 0.f;
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

   const float PMMAbreakE = ToSI::eV(45.f);
   const float PMMAdensity = 1190.f;
   const float PMMAworkfun = 5.5f;
   const float PMMAbandgap = 5.f;
   const float PMMAEFermi = -PMMAbandgap;
   const float PMMApotU = -PMMAworkfun - PMMAEFermi;

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
   const float glCpotU = -glCworkfun - glCEFermi;

   //SEmaterialT* glC = nullptr;

   //SelectableElasticSMT* glCNISTMott = nullptr;

   //TabulatedInelasticSMT* glCDS = nullptr;

   //ExpQMBarrierSMT* glceqmbsm = nullptr;

   //MONSEL_MaterialScatterModelT* glCMSM = nullptr;

   //MONSEL_MaterialScatterModelT* glCMSMDeep = nullptr;

   const float SiphononE = 0.063f;
   const float SiphononStrength = 3.f;
   const float Sidensity = 2330.f;
   const float Siworkfun = 4.85f;
   const float Sibandgap = 1.1f;
   const float SiEFermi = -Sibandgap;
   const float SipotU = -Siworkfun - SiEFermi;//-Siworkfun - SiEFermi;

   //SEmaterialT* Si = nullptr;

   //SelectableElasticSMT* SiNISTMott = nullptr;

   //TabulatedInelasticSMT* SiDS = nullptr;

   //GanachaudMokraniPhononInelasticSMT* Siphonon = nullptr;

   //ExpQMBarrierSMT* sieqmbsm = nullptr;

   //MONSEL_MaterialScatterModelT* SiMSM = nullptr;

   //MONSEL_MaterialScatterModelT* SiMSMDeep = nullptr;

   const float SiO2density = 2200.f;
   const float SiO2workfun = 10.f;
   const float SiO2phononStrength = 2.f; // Number of phonon modes
   const float SiO2phononE = 0.145f; // Phonon mode energy in eV
   const float SiO2bandgap = 8.9f; // width of band gap in eV
   const float SiO2EFermi = -SiO2bandgap; // This puts the Fermi level at the top of the valence band.
   const float SiO2potU = -SiO2workfun - SiO2EFermi;

   //SphereT* sphere = nullptr;

   const float pitch = pitchnm * 1.e-9f;
   const float h = hnm * 1.e-9f;
   const float w = wnm * 1.e-9f;
   float linelength = linelengthnm * 1.e-9f;

   const float radperdeg = Math2::PI / 180.f;
   const float thetar = thetardeg * Math2::PI / 180.f;
   const float thetal = thetaldeg * Math2::PI / 180.f;
   const float radr = radrnm * 1.e-9f;
   const float radl = radlnm * 1.e-9f;
   float layerthickness = layerthicknessnm * 1.e-9f;
   const float layer1thickness = layer1thicknessnm * 1.e-9f;
   const float layer2thickness = layer2thicknessnm * 1.e-9f;
   float beamsize = beamsizenm * 1.e-9f;
   float beamz = beamznm * 1.e-9f;
   const float deep = deepnm * 1.e-9f;

   const double center[] = {
      0.0,
      0.0,
      0.0
   };

   float beamEeV = 3000.f;
   //float beamEeV = 0.f;
   float beamE = ToSI::eV(beamEeV);
   const float binSizeEV = 10.f;
   const float cutoffEnergyForSE = 50.f;
   float beamphideg = 0.f;
   float beamthetadeg = 20.f;
   float beamfocallength = 20.f;

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

   const char* glCTables[] = {
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\IIMFPPennInterpglassyCSI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\interpNUSimReducedDeltaEglassyCSI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\interpsimTableThetaNUglassyCSI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\interpSimESE0NUglassyCSI.csv"
   };

   const char* SiTables[] = {
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\IIMFPFullPennInterpSiSI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\interpNUSimReducedDeltaEFullPennSiSI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\interpNUThetaFullPennSiBGSI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\interpSimESE0NUSiBGSI.csv"
   };

   const char* SiO2Tables[] = {
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\IIMFPPennInterpSiO2SI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\interpNUSimReducedDeltaESiO2SI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\interpsimTableThetaNUSiO2SI.csv",
      "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiO2Tables\\interpSimESE0NUSiO2SI.csv"
   };

   const ElementT* COHcomponents[] = { &Element::C, &Element::O, &Element::H };
   const Composition::data_type COHcompositions[] = { 5.f, 2.f, 8.f };

   const ElementT* glCComponents[] = { &Element::C };
   const Composition::data_type glCComposition[] = { 1.f };

   const ElementT* SiComponent[] = { &Element::Si };
   const Composition::data_type SiComposition[] = { 1. };

   const ElementT* SiO2Components[] = { &Element::Si, &Element::O };

   unsigned int ysize, xsize;
   float xstartnm, xstopnm, ystartnm, ystopnm;

   NShapes::LineParams** lineParams;

   NShapes::HorizontalStripParams** hstripParams;
   unsigned int nhstrips;
#endif

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
      NUTableInterpolation::getInstance(glCTables[0]);
      NUTableInterpolation::getInstance(glCTables[1]);
      NUTableInterpolation::getInstance(glCTables[2]);
      NUTableInterpolation::getInstance(glCTables[3]);

      NUTableInterpolation::getInstance(SiTables[0]);
      NUTableInterpolation::getInstance(SiTables[1]);
      NUTableInterpolation::getInstance(SiTables[2]);
      NUTableInterpolation::getInstance(SiTables[3]);

      NUTableInterpolation::getInstance(SiO2Tables[0]);
      NUTableInterpolation::getInstance(SiO2Tables[1]);
      NUTableInterpolation::getInstance(SiO2Tables[2]);
      NUTableInterpolation::getInstance(SiO2Tables[3]);
   }

   void transferNUTableToCuda()
   {
      NUTableInterpolation::initFactory << <1, 1 >> >();
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());

      NUTableInterpolation::transferDataToCuda(glCTables[0]);
      NUTableInterpolation::transferDataToCuda(glCTables[1]);
      NUTableInterpolation::transferDataToCuda(glCTables[2]);
      NUTableInterpolation::transferDataToCuda(glCTables[3]);

      char* d_fn = nullptr;

      checkCudaErrors(cudaMalloc((void **)&d_fn, strnlen_s(glCTables[0], 256) * sizeof(char)));
      checkCudaErrors(cudaMemset(d_fn, 0, strnlen_s(glCTables[0], 256) * sizeof(char)));
      checkCudaErrors(cudaMemcpy(d_fn, glCTables[0], strnlen_s(glCTables[0], 256) * sizeof(char), cudaMemcpyHostToDevice));
      verifyNUTable1d << <1, 1 >> >(d_fn);
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      checkCudaErrors(cudaFree(d_fn));
      const NUTableInterpolationT* table0 = NUTableInterpolation::getInstance(glCTables[0]);
      const VectorXf& data0 = table0->gettable1d();
      printf("CPU %s\n", glCTables[0]);
      for (auto v : data0) {
         printf("%.5e ", v);
      }
      printf("\n");

      checkCudaErrors(cudaMalloc((void **)&d_fn, strnlen_s(glCTables[1], 256) * sizeof(char)));
      checkCudaErrors(cudaMemset(d_fn, 0, strnlen_s(glCTables[1], 256) * sizeof(char)));
      checkCudaErrors(cudaMemcpy(d_fn, glCTables[1], strnlen_s(glCTables[1], 256) * sizeof(char), cudaMemcpyHostToDevice));
      verifyNUTable2d << <1, 1 >> >(d_fn, 0);
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      checkCudaErrors(cudaFree(d_fn));
      const NUTableInterpolationT* table1 = NUTableInterpolation::getInstance(glCTables[1]);
      const MatrixXf& data1 = table1->gettable2d();
      printf("CPU %s: row %d\n", glCTables[1], 0);
      for (auto v : data1[0]) {
         printf("%.5e ", v);
      }
      printf("\n");

      checkCudaErrors(cudaMalloc((void **)&d_fn, strnlen_s(glCTables[2], 256) * sizeof(char)));
      checkCudaErrors(cudaMemset(d_fn, 0, strnlen_s(glCTables[2], 256) * sizeof(char)));
      checkCudaErrors(cudaMemcpy(d_fn, glCTables[2], strnlen_s(glCTables[2], 256) * sizeof(char), cudaMemcpyHostToDevice));
      verifyNUTable3d << <1, 1 >> >(d_fn, 50, 50);
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      checkCudaErrors(cudaFree(d_fn));
      const NUTableInterpolationT* table2 = NUTableInterpolation::getInstance(glCTables[2]);
      const Matrix3DXf& data2 = table2->gettable3d();
      printf("CPU %s: row %d, col %d\n", glCTables[2], 50, 50);
      for (auto v : data2[50][50]) {
         printf("%.5e ", v);
      }
      printf("\n");

      NUTableInterpolation::transferDataToCuda(SiTables[0]);
      NUTableInterpolation::transferDataToCuda(SiTables[1]);
      NUTableInterpolation::transferDataToCuda(SiTables[2]);
      NUTableInterpolation::transferDataToCuda(SiTables[3]);
   }

   __host__ __device__ void initImgRange()
   {
      //VectorXf yvalstmp(128);
      //for (int i = -64; i < 64; i += 1) {
      //   yvalstmp.push_back(i);
      //}
      ////VectorXf yvalstmp(1, -64);

      //const float xbottom = wnm / 2.f;
      //const float xtop = wnm / 2.f - hnm * ::tanf(thetar);
      //xstartnm = xbottom - 200.f;
      //xstopnm = xbottom + 200.f;
      //const float xfinestart = xtop - 20.5f;
      //const float xfinestop = (thetar < 0.f) ? xtop + 20.5f : wnm / 2.f + 20.5f;

      xstartnm = ToSI::GIGA * (lineParams[0]->x - lineParams[0]->w / 2. - lineParams[0]->radl) - 10. * Random::random();
      xstopnm = ToSI::GIGA * (lineParams[nlines - 1]->x + lineParams[nlines - 1]->w / 2. + lineParams[nlines - 1]->radr) + 10. * Random::random();

      float minz = 0.f;
      for (int i = 0; i < nlines; ++i) {
         minz = -lineParams[i]->h < minz ? -lineParams[i]->h : minz;
         ystartnm = ToSI::GIGA * minz - Random::random() * 20.f;
      }
      ystopnm = ToSI::GIGA * (hstripParams[nhstrips - 1]->y + hstripParams[nhstrips - 1]->w / 2.) + (Random::random() - 0.5) * 2. * 10.f;

      xsize = 512;
      ysize = 512;

      //VectorXf xvalstmp(80);
      //float deltax = 5.f;
      //float x = xstart;
      //while (x < xfinestart) {
      //   xvalstmp.push_back(x);
      //   x += deltax;
      //}
      //x = xfinestart;
      //deltax = 1.f;
      //while (x < xfinestop) {
      //   xvalstmp.push_back(x);
      //   x += deltax;
      //}
      //x = xfinestop;
      //deltax = 5.f;
      //while (x < xstop) {
      //   xvalstmp.push_back(x);
      //   x += deltax;
      //}
      //xvalstmp.push_back(xstop);
      ////VectorXf xvalstmp(2);
      ////xvalstmp.push_back(xstart);
      ////xvalstmp.push_back(xstart + 5.f);

      if (ystopnm - ystartnm < 0) printf("initRange(): ystopnm - ystartnm < 0\n");
      if (xstopnm - xstartnm < 0) printf("initRange(): xstopnm - xstartnm < 0\n");
   }

   __global__ void initCuda()
   {
      printf("LinesOnLayers: initCuda\n");
      for (int i = 0; i < 10; ++i) {
         printf("%.10e\n", Random::random());
      }

      initImgRange();

      printf("(%d, %d)", xsize, ysize);
   }

   enum MaterialTypes
   {
      PMMA = 0,
      Si = 1,
      SiO2 = 2,
   };

   __host__ __device__ void runSinglePixel(const unsigned int r, const unsigned int c, float* bse, float* fse, float* totalse)
   {
      const float deltay = (ystopnm - ystartnm) / ysize;
      const float deltax = (xstopnm - xstartnm) / xsize;

      const float ynm = (ystartnm + r * deltay); 
      const float y = ynm * 1.e-9f;
      const float xnm = (xstartnm + c * deltax);
      const float x = xnm * 1.e-9f;

      SEmaterialT vacuum; // TODO: move this to global
      vacuum.setName("SE vacuum");
      ExpQMBarrierSMT vacuumBarrier(&vacuum); // TODO: move this to global
      ZeroCSDT sZeroCSD; // TODO: move this to global

      MONSEL_MaterialScatterModelT vacuumMSM(&vacuum, &vacuumBarrier, &sZeroCSD); // TODO: move this to global

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      COHcomponents[0] = Element::dC;
      COHcomponents[1] = Element::dO;
      COHcomponents[2] = Element::dH;
#endif

      CompositionT PMMAcomp;
      PMMAcomp.defineByMoleFraction(COHcomponents, 3, COHcompositions, 3);
      SEmaterialT PMMA(PMMAcomp, PMMAdensity); // TODO: move this to global
      PMMA.setName("PMMA");
      PMMA.setWorkfunction(ToSI::eV(PMMAworkfun));
      PMMA.setBandgap(ToSI::eV(PMMAbandgap));
      PMMA.setEnergyCBbottom(ToSI::eV(PMMApotU));

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      SelectableElasticSMT PMMANISTMott(PMMA, *NISTMottRS::dFactory);
#else
      SelectableElasticSMT PMMANISTMott(PMMA, NISTMottRS::Factory);
#endif

      JoyLuoNieminenCSDT PMMACSD(PMMA, PMMAbreakE);
      FittedInelSMT PMMAfittedInel(PMMA, ToSI::eV(65.4f), PMMACSD); // TODO: move this to global
      GanachaudMokraniPolaronTrapSMT PMMApolaron(2.e7f, 1.f / ToSI::eV(4.f)); // TODO: move this to global

      ExpQMBarrierSMT pmmaeqmbsm(&PMMA); // TODO: move this to global

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
      PMMAMSMDeep.setMinEforTracking(ToSI::eV(cutoffEnergyForSE));

      MONSEL_MaterialScatterModelT& ARCMSM = PMMAMSM;

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      glCComponents[0] = Element::dC;
#endif

      SEmaterialT glC(glCComponents, 1, glCComposition, 1, glCdensity, "glassy Carbon");
      glC.setWorkfunction(ToSI::eV(glCworkfun));
      glC.setEnergyCBbottom(ToSI::eV(glCpotU));
      glC.setBandgap(ToSI::eV(glCbandgap));
      const Composition::data_type glCCoreEnergy[] = { ToSI::eV(284.2f) };
      glC.setCoreEnergy(glCCoreEnergy, 1);

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      SelectableElasticSMT glCNISTMott(glC, *NISTMottRS::dFactory);
#else
      SelectableElasticSMT glCNISTMott(glC, NISTMottRS::Factory);
#endif
      TabulatedInelasticSMT glCDS(glC, 3, glCTables);

      ExpQMBarrierSMT glceqmbsm(&glC);

      MONSEL_MaterialScatterModelT glCMSM(&glC, &glceqmbsm, &sZeroCSD);
      glCMSM.addScatterMechanism(&glCNISTMott);
      glCMSM.addScatterMechanism(&glCDS);

      MONSEL_MaterialScatterModelT glCMSMDeep(&glC, &glceqmbsm, &sZeroCSD);
      glCMSMDeep.addScatterMechanism(&glCNISTMott);
      glCMSMDeep.addScatterMechanism(&glCDS);

      glCMSMDeep.setMinEforTracking(ToSI::eV(cutoffEnergyForSE));

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      SiComponent[0] = Element::dSi;
#endif

      SEmaterialT Si(SiComponent, 1, SiComposition, 1, Sidensity, "Silicon");
      Si.setWorkfunction(ToSI::eV(Siworkfun));
      Si.setEnergyCBbottom(ToSI::eV(SipotU));
      Si.setBandgap(ToSI::eV(Sibandgap));
      const Material::data_type SiCoreEnergy[] = { ToSI::eV(99.2f), ToSI::eV(99.8f), ToSI::eV(149.7f), ToSI::eV(1839.f) };
      Si.setCoreEnergy(SiCoreEnergy, 4);

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      SelectableElasticSMT SiNISTMott(Si, *NISTMottRS::dFactory);
#else
      SelectableElasticSMT SiNISTMott(Si, NISTMottRS::Factory);
#endif
      TabulatedInelasticSMT SiDS(Si, 3, SiTables, ToSI::eV(13.54f));

      GanachaudMokraniPhononInelasticSMT Siphonon(SiphononStrength, ToSI::eV(SiphononE), 300.f, 11.7f, 1.f);

      ExpQMBarrierSMT sieqmbsm(&Si);

      MONSEL_MaterialScatterModelT SiMSM(&Si, &sieqmbsm, &sZeroCSD);
      SiMSM.addScatterMechanism(&SiNISTMott);
      SiMSM.addScatterMechanism(&SiDS);
      SiMSM.addScatterMechanism(&Siphonon);

      MONSEL_MaterialScatterModelT SiMSMDeep(&Si, &sieqmbsm, &sZeroCSD);
      SiMSMDeep.addScatterMechanism(&SiNISTMott);
      SiMSMDeep.addScatterMechanism(&SiDS);
      SiMSMDeep.addScatterMechanism(&Siphonon);

      SiMSMDeep.setMinEforTracking(ToSI::eV(cutoffEnergyForSE));

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      const float SiWeight = Element::dSi->getAtomicWeight();
      const float OxWeight = 2.f * Element::dO->getAtomicWeight();
#else
      const float SiWeight = Element::Si.getAtomicWeight();
      const float OxWeight = 2.f * Element::O.getAtomicWeight();
#endif
      const float totalWeight = SiWeight + OxWeight;

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      SiO2Components[0] = Element::dSi;
      SiO2Components[1] = Element::dO;
#endif

      const float SiO2massFrac[] = { SiWeight / totalWeight, OxWeight / totalWeight };
      SEmaterialT SiO2(SiO2Components, 2, SiO2massFrac, 2, SiO2density, "Silicon Dioxide");

      SiO2.setWorkfunction(ToSI::eV(SiO2workfun));
      SiO2.setBandgap(ToSI::eV(SiO2bandgap));
      SiO2.setEnergyCBbottom(ToSI::eV(SiO2potU));
      const float SiO2ce[] = { ToSI::eV(41.6), ToSI::eV(99.2), ToSI::eV(99.8), ToSI::eV(149.7), ToSI::eV(543.1), ToSI::eV(1839.) };
      SiO2.setCoreEnergy(SiO2ce, 6);

      // Create scatter mechanisms
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      SelectableElasticSMT SiO2NISTMott((MaterialT&)SiO2, *NISTMottRS::dFactory);
#else
      SelectableElasticSMT SiO2NISTMott((MaterialT&)SiO2, NISTMottRS::Factory);
#endif

      TabulatedInelasticSMT SiO2DS(SiO2, 3, SiO2Tables, ToSI::eV(20. + SiO2bandgap));
      GanachaudMokraniPhononInelasticSMT SiO2phonon(SiO2phononStrength, ToSI::eV(SiO2phononE), 300., 3.82, 1.);
      GanachaudMokraniPolaronTrapSMT SiO2polaron(1.0e9, 1. / ToSI::eV(1.));
      ExpQMBarrierSMT SiO2ClassicalBarrier(&SiO2); // The default.No need to actually execute this line.
      //SiO2CSD = mon.ZeroCSD(); // The default.No need to actually execute this line.
      // Make a material scatter model
      // MSM to be used in thin layer(includes SE generation)
      MONSEL_MaterialScatterModelT SiO2MSM(&SiO2, &SiO2ClassicalBarrier, &sZeroCSD);
      SiO2MSM.addScatterMechanism(&SiO2NISTMott);
      SiO2MSM.addScatterMechanism(&SiO2DS);
      SiO2MSM.addScatterMechanism(&SiO2phonon);
      // SiO2MSM.addScatterMechanism(SiO2polaron);
      //SiO2MSM.setCSD(SiO2CSD);
      //SiO2MSM.setBarrierSM(SiO2ClassicalBarrier);

      // MSM to be used deep inside(drops electrons with E < 50 eV);
      //MONSEL_MaterialScatterModelT SiO2MSMDeep(&SiO2, &SiO2ClassicalBarrier, &sZeroCSD);
      //SiO2MSMDeep.addScatterMechanism(&SiO2NISTMott);
      //SiO2MSMDeep.addScatterMechanism(&SiO2DS);
      //SiO2MSMDeep.addScatterMechanism(&SiO2phonon);
      ////SiO2MSMDeep.setCSD(SiO2CSD);
      ////SiO2MSMDeep.setBarrierSM(SiO2ClassicalBarrier);
      //SiO2MSMDeep.setMinEforTracking(ToSI::eV(50.));

      SphereT sphere(center, MonteCarloSS::ChamberRadius);

      NullMaterialScatterModelT NULL_MSM;

      RegionT chamber(nullptr, &NULL_MSM, &sphere);
      chamber.updateMaterial(*(chamber.getScatterModel()), vacuumMSM);

      NShapes::HorizontalStrip** strips = new NShapes::HorizontalStrip*[nhstrips];
      RegionT** stripRegions = new RegionT*[nhstrips];
      for (int i = 0; i < nhstrips; ++i) {
         strips[i] = new NShapes::HorizontalStrip(hstripParams[i]->w, hstripParams[i]->fadebot);
         NormalMultiPlaneShapeT* layer = strips[i]->get();
         double dist[3] = { 0., hstripParams[i]->y, 0. };
         layer->translate(dist);

         //stripRegions[i] = new RegionT(&chamber, &PMMAMSM, (NormalShapeT*)layer);

         switch (hstripParams[i]->material)
         {
         case MaterialTypes::SiO2:
            stripRegions[i] = new RegionT(&chamber, &SiO2MSM, (NormalShapeT*)layer);
            break;
         case MaterialTypes::Si:
            stripRegions[i] = new RegionT(&chamber, &SiMSM, (NormalShapeT*)layer);
            break;
         case MaterialTypes::PMMA:
            stripRegions[i] = new RegionT(&chamber, &PMMAMSM, (NormalShapeT*)layer);
            break;
         }
      }

      //NormalMultiPlaneShapeT layer1;
      //PlaneT pl1(normalvector, layer1Pos);
      //layer1.addPlane(&pl1);
      //RegionT layer1Region(&chamber, &ARCMSM, (NormalShapeT*)&layer1);

      //NormalMultiPlaneShapeT layer2;
      //PlaneT pl2(normalvector, layer2Pos);
      //layer2.addPlane(&pl2);
      //RegionT layer2Region(&layer1Region, &glCMSM, (NormalShapeT*)&layer2);

      //NormalMultiPlaneShapeT layer3;
      //PlaneT pl3(normalvector, layer3Pos);
      //layer3.addPlane(&pl3);
      //RegionT layer3Region(&layer2Region, &SiMSM, (NormalShapeT*)&layer3);

      //NormalMultiPlaneShapeT layer4;
      //PlaneT pl4(normalvector, layer4Pos);
      //RegionT layer4Region(&layer3Region, &SiMSM, (NormalShapeT*)&layer4);

      //RegionT deepRegion(&layer3Region, &SiMSMDeep, (NormalShapeT*)&layer4);

      //NShapes::Line line(-h, w, linelength, thetal, thetar, radl, radr);
      //const double pivot[3] = { 0.f, 0.f, 0.f };
      //line.get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
      //const double dist1[3] = { 0.f, 0.f, linelength / 2. };
      //line.get()->translate(dist1);

      //RegionT lineRegion(&chamber, &PMMAMSM, (NormalIntersectionShapeT*)line.get());

      //NShapes::Line* lines[3];
      //RegionT* regions[3];
      //for (int i = 0; i < nlines; ++i) {
      //   //NShapes::Line line(-h, w, linelength, thetal, thetar, radl, radr);
      //   lines[i] = new NShapes::Line(-lineParams[i]->h, lineParams[i]->w, lineParams[i]->linelength, lineParams[i]->thetal, lineParams[i]->thetar, lineParams[i]->radl, lineParams[i]->radr);
      //   const double pivot[3] = { 0.f, 0.f, 0.f };
      //   lines[i]->get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
      //   //lines[i]->get()->rotate(pivot, -Math2::PI / 2.f + 20.f / (Math2::PI / 2.f), Math2::PI / 2.f - 20.f / (Math2::PI / 2.f), Math2::PI / 2.f - 20.f / (Math2::PI / 2.f));
      //   //const double dist1[3] = { 0.f, 0.f, linelength / 2. };
      //   const double offset[3] = { lineParams[i]->x, 0.f, linelength / 2. };
      //   lines[i]->get()->translate(offset);
      //   //lines[i]->addRestrainingPlanes();
      //   regions[i] = new RegionT(&chamber, &PMMAMSM, (NormalIntersectionShapeT*)lines[i]->get());
      //}

      NShapes::Line** lines = new NShapes::Line*[nlines];
      RegionT** lineRegions = new RegionT*[nlines];
      for (int i = 0; i < nlines; ++i) {
         //NShapes::Line line(-h, w, linelength, thetal, thetar, radl, radr);
         lines[i] = new NShapes::Line(-lineParams[i]->h, lineParams[i]->w, lineParams[i]->linelength, lineParams[i]->thetal, lineParams[i]->thetar, lineParams[i]->radl, lineParams[i]->radr);
         const double pivot[3] = { 0.f, 0.f, 0.f };
         lines[i]->get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
         //lines[i]->get()->rotate(pivot, -Math2::PI / 2.f + 20.f / (Math2::PI / 2.f), Math2::PI / 2.f - 20.f / (Math2::PI / 2.f), Math2::PI / 2.f - 20.f / (Math2::PI / 2.f));
         //const double dist1[3] = { 0.f, 0.f, linelength / 2. };
         const double offset[3] = { lineParams[i]->x, 0.f, linelength / 2. };
         lines[i]->get()->translate(offset);
         //lines[i]->addRestrainingPlanes();
         switch (lineParams[i]->material) {
         case MaterialTypes::PMMA:
            lineRegions[i] = new RegionT(&chamber, &PMMAMSM, (NormalIntersectionShapeT*)lines[i]->get());
            break;
         case MaterialTypes::Si:
            lineRegions[i] = new RegionT(&chamber, &SiMSM, (NormalIntersectionShapeT*)lines[i]->get());
            break;
         case MaterialTypes::SiO2:
            lineRegions[i] = new RegionT(&chamber, &SiO2MSM, (NormalIntersectionShapeT*)lines[i]->get());
            break;
         }
         //regions[i] = new RegionT(&chamber, curmat, (NormalIntersectionShapeT*)lines[i]->get());
         //regions[i] = new RegionT(&chamber, &SiMSM, (NormalIntersectionShapeT*)lines[i]->get());
      }

      //NShapes::Line line(-h, w, linelength, thetal, thetar, radl, radr);
      //const double pivot[3] = { 0.f, 0.f, 0.f };
      //line.get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
      //const double offset[3] = { lineParams[0]->x, 0.f, linelength / 2. };
      //line.get()->translate(offset);
      //RegionT region(&chamber, &PMMAMSM, (NormalIntersectionShapeT*)line.get());

      //NShapes::Washer washer(w / 10., w / 5.);
      //RegionT region(&chamber, &PMMAMSM, (NormalDifferenceShapeT*)washer.get());
      
      const double egCenter[] = { x, y, beamz};
      GaussianBeamT eg(beamsize, beamE, Math2::toRadians(beamthetadeg), Math2::toRadians(beamphideg), egCenter, beamfocallength);
      MonteCarloSST monte(&eg, &chamber, nullptr);

      const int nbins = (int)(beamEeV / binSizeEV);
      BackscatterStatsT back(monte, nbins);
      monte.addActionListener(back);

//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      monte.runMultipleTrajectories(nTrajectories);
//#else
//      try {
//         monte.runMultipleTrajectories(nTrajectories);
//      }
//      catch (std::exception ex) {
//         printf("runSinglePixel: %s\n", ex.what());
//      }
//#endif

      monte.runMultipleTrajectories(nTrajectories);

      const HistogramT& histBSE = back.backscatterEnergyHistogram();
      const float energyperbineVBSE = beamEeV / histBSE.binCount();
      const float maxSEbinBSE = cutoffEnergyForSE / energyperbineVBSE;
      int totalBSE = 0;
      for (int j = 0; j < (int)maxSEbinBSE; ++j) {
         totalBSE += histBSE.counts(j);
      }
      const float BSEf = (float)totalBSE / nTrajectories;

      const HistogramT& histFSE = back.forwardscatterEnergyHistogram();
      const float energyperbineVFSE = beamEeV / histFSE.binCount();
      const float maxSEbinFSE = cutoffEnergyForSE / energyperbineVFSE;
      int totalFSE = 0;
      for (int j = 0; j < (int)maxSEbinFSE; ++j) {
         totalFSE += histFSE.counts(j);
      }
      const float FSEf = (float)totalFSE / nTrajectories;

      const float bsf = back.backscatterFraction() - BSEf;
      //printf("%lf %lf %lf %lf %lf\n", beamEeV, xnm, ynm, bsf, SEf);
      monte.removeActionListener(back);

      bse[r * xsize + c] = BSEf;
      fse[r * xsize + c] = FSEf;
      totalse[r * xsize + c] = BSEf + FSEf;

      // clean up
      for (int i = 0; i < nhstrips; ++i) {
         delete strips[i];
         delete stripRegions[i];
      }
      delete strips;
      delete stripRegions;

      for (int i = 0; i < nlines; ++i) {
         delete lines[i];
         delete lineRegions[i];
      }
      delete lines;
      delete lineRegions;
   }

   __global__ void
      //__launch_bounds__(256, 3)
      runCuda(float* bse, float* fse, float* totalse)
   {
      const unsigned int r = blockIdx.y*blockDim.y + threadIdx.y;
      const unsigned int c = blockIdx.x*blockDim.x + threadIdx.x;
      if (r >= ysize || c >= xsize) return;

      const unsigned int blockId = blockIdx.x + blockIdx.y * gridDim.x;
      const unsigned int threadId = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;
      printf("%d, %d (%d) began\n", r, c, threadId);

      runSinglePixel(r, c, bse, fse, totalse);

      printf("%d, %d (%d) ended\n", r, c, threadId);
   }

   void runSinglePixelThread(int id, const unsigned int r, const unsigned int c, float* bse, float* fse, float* totalse)
   {
      try {
         runSinglePixel(r, c, bse, fse, totalse);
      }
      catch (std::exception ex) {
         printf("%s\n", ex.what());
      }
   }

   inline float getAdjustedValPlusMinus(float val, float fraction)
   {
      return val * (1.f + (1.f - Random::random() * 2.f) * fraction); // val +/- val*fraction
   }

   //__host__ __device__ void ResetParams()
   //{
   //   strip0Width = 30.f * 1.e-9f;
   //   strip1Width = 30.f * 1.e-9f;
   //   strip2Width = 30.f * 1.e-9f;
   //   fade2 = false;

   //   nlines = 3;
   //   linemat = 0;

   //   nTrajectories = 100;

   //   beamEeV = 3000.f;
   //}

   __host__ __device__ void setSimParams()
   {
      // hrizontal strips
      //nhstrips = 3 + Random::randomInt(2);
      nhstrips = 1;
      hstripParams = new NShapes::HorizontalStripParams*[nhstrips];
      layerthickness = 30.f * 1e-9f;
      // generate strips
      //unsigned int mat = Random::randomInt(3);
      unsigned int mat = 0; // PMMA

      for (int i = 0; i < nhstrips; ++i) {
         //switch (i) {
         //case MaterialTypes::PMMA: mat = MaterialTypes::PMMA; break;
         //case MaterialTypes::Si: mat = MaterialTypes::Si; break;
         //case MaterialTypes::SiO2: mat = MaterialTypes::SiO2; break;
         //}

         hstripParams[i] = new NShapes::HorizontalStripParams(layerthickness * (1 + (Random::random() - .5f)), i == nhstrips - 1 ? Random::random() > .5f : false, false, mat, 0.f, 0.f);
      }

      // translate
      float cury = hstripParams[0]->w / 2.f;
      hstripParams[0]->y = cury;
      for (int i = 1; i < nhstrips; ++i) {
         cury += hstripParams[i - 1]->w / 2.f + hstripParams[i]->w / 2.f;
         hstripParams[i]->y = cury;
      }

      //// TODO: delete tmp code
      //{
      //   hstripParams[1]->y = hstripParams[0]->y;
      //   hstripParams[1]->z = hstripParams[0]->z;

      //   hstripParams[0]->y = (-100.f + (1.f - Random::random()) * 10.f) * ToSI::NANO;
      //   hstripParams[0]->z = linelength;
      //}
      // feature lines
      nlines = 3 + Random::randomInt(3);
      //nlines = 2;
      //linemat = Random::randomInt(3);
      linemat = mat;
      lineParams = new NShapes::LineParams*[nlines];
      // generate line shape
      const float h0 = h * (.5f + Random::random()); // 50% to 150%
      const float w0 = w + w * 0.1f * 2.f * (Random::random() - 0.5f); // +/- 10%
      for (int i = 0; i < nlines; ++i) {
         const float curh = h0 * (1.f + Random::random() * 0.01f); // 100% to 105%
         const float curw = w0 * (1.f + Random::random() * 0.05f); // 100% to 105%
         const float curl = linelength;
         const float curtl = thetal + Math2::toRadians(Random::random() / 2.f); // + 0 to .5 deg
         const float curtr = thetar + Math2::toRadians(Random::random() / 2.f); // + 0 to .5 deg
         const float currl = curw / 3.f + curw / 3.f * (Random::random() - .5f);
         const float currr = curw / 3.f + curw / 3.f * (Random::random() - .5f);

         lineParams[i] = new NShapes::LineParams(curh, curw, curl, curtl, curtr, currl, currr, linemat, 0.f);
      }
      // update line position
      float curx = -10.f * ToSI::NANO;
      lineParams[0]->x = curx;
      for (int i = 1; i < nlines; ++i) {
         curx += lineParams[i - 1]->w / 2.f + lineParams[i - 1]->radr + lineParams[i]->w / 2.f + lineParams[i]->radl;
         lineParams[i]->x = curx;
      }

      //nTrajectories += 250;
      nTrajectories = 50 + Random::random() * 150;
      //nTrajectories = 100;

      beamEeV = 100.f + 400.f * Random::random();
      beamE = ToSI::eV(beamEeV);

      //beamsizenm += 0.1f;
      beamsizenm = .25f + .5f * Random::random();
      beamsize = beamsizenm * ToSI::NANO; // 0.1 nm to 0.5 nm
      beamthetadeg = -5.f  + Random::random() + 10.f; // polar
      beamphideg = -10.f  + 20.f * Random::random(); // azimuth
      beamz = -h0 - 20.f / ToSI::GIGA;
      beamfocallength = (-beamz - h0) / ::cos(Math2::toRadians(beamthetadeg));

      printf("nlines: %d\n", nlines);
      printf("linemat: %d\n", linemat);
      printf("nTrajectories: %d\n", nTrajectories);
      printf("beamEeV: %.5e\n", beamEeV);
      printf("beamsizenm: %.5e\n", beamsizenm);
      printf("beamz: %.5e\n", beamz);
      printf("beamfocallength: %.5e\n", beamfocallength);
   }

   //void setSimParamsFromCSVFile(const char* fname)
   //{
   //   std::fstream fin(fname, std::ios::in);

   //   std::string line, word, temp;

   //   while (fin >> temp) {
   //      getline(fin, line);

   //      std::stringstream s(line);
   //      int c = 0;
   //      while (getline(s, word, ',')) {
   //         float val = std::stof(word.c_str());

   //         switch (c) {
   //         case 0:
   //         {
   //            break;
   //         }
   //         default:
   //         {
   //            break;
   //         }
   //         }

   //         ++c;
   //      }
   //   }
   //   fin.close();
   //}

   __host__ __device__ void destroySimParams()
   {
      for (int i = 0; i < nlines; ++i) {
         delete lineParams[i];
      }
      delete[] lineParams;

      for (int i = 0; i < nhstrips; ++i) {
         delete hstripParams[i];
      }
      delete[] hstripParams;
   }

   void lineProjection(const unsigned int idx, char* gt)
   {
      const double p[3] = { xstartnm * 1.e-9, ystartnm * 1.e-9, beamz };
      static const double n[3] = { 0., 0., 1. };
      PlaneT projectionPlane(n, p); // image plane (where light source will be moving on)
      static const double axis0[3] = { 1., 0., 0. }; // relative spanning vector from plane origin
      static const double axis1[3] = { 0., 1., 0. }; // relative spanning vector from plane origin

      const float xlenperpix = (xstopnm - xstartnm) / xsize * 1.e-9; // resolution m / pix
      const float ylenperpix = (ystopnm - ystartnm) / ysize * 1.e-9; // resolution m / pix

      static const double st = ::sinf(Math2::toRadians(beamthetadeg));
      static const double beamdir[3] = {
         ::cosf(Math2::toRadians(beamphideg)) * st,
         ::sinf(Math2::toRadians(beamphideg)) * st,
         ::cosf(Math2::toRadians(beamthetadeg))
      };
      
      // horizontal strips
      for (int i = 0; i < nhstrips; ++i) {
         NShapes::HorizontalStrip hs(hstripParams[i]->w);
         NormalMultiPlaneShapeT* strip = hs.get();
         double offset[3] = { 0., hstripParams[i]->y, hstripParams[i]->z };
         strip->translate(offset);

         hs.calcGroundtruth();
         hs.calcRasterization(projectionPlane, axis0, axis1, beamdir, xlenperpix, ylenperpix, gt, xsize, ysize);
      }

      // lines
      for (int i = 0; i < nlines; ++i) {
         //NShapes::Line line(-h, w, linelength, thetal, thetar, radl, radr);
         NShapes::Line line(-lineParams[i]->h, lineParams[i]->w, lineParams[i]->linelength, lineParams[i]->thetal, lineParams[i]->thetar, lineParams[i]->radl, lineParams[i]->radr);
         //NShapes::Line line(-lp.h, lp.w, lp.linelength, lp.thetal, lp.thetar, lp.radl, lp.radr);
         const double pivot[3] = { 0.f, 0.f, 0.f };
         line.get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
         //line.get()->rotate(pivot, -Math2::PI / 2.f + 20.f / (Math2::PI / 2.f), Math2::PI / 2.f - 20.f / (Math2::PI / 2.f), Math2::PI / 2.f - 20.f / (Math2::PI / 2.f));
         //const double dist1[3] = { 0.f, 0.f, linelength / 2. };
         //const double dist1[3] = { lp.x, 0.f, linelength / 2. };
         const double offset[3] = { lineParams[i]->x, 0.f, linelength / 2. };
         line.get()->translate(offset);

         line.calcGroundtruth(); // get points/line segments that need to be projected
         line.calcRasterization(projectionPlane, axis0, axis1, beamdir, xlenperpix, ylenperpix, gt, xsize, ysize); // calculation for line needs to be last since bottom needs to be removed due to same material
      }

      // corrections, separated from above so that it is more robust to changing of code order
      // corrections are needed bcuz in the case where two regions have the same material, the groundtruths has to change accordingly
      // horizontal strips
      for (int i = 0; i < nhstrips; ++i) {
         NShapes::HorizontalStrip hs(hstripParams[i]->w);
         NormalMultiPlaneShapeT* layer = hs.get();
         double offset[3] = { 0., hstripParams[i]->y, hstripParams[i]->z };
         layer->translate(offset);

         hs.calcGroundtruth();

         // boundary corrections
         if (i > 0 && hstripParams[i]->material == hstripParams[i - 1]->material) { // bottom needs to be removed due to same material
            hs.calcRasterizationCorrection(projectionPlane, axis0, axis1, beamdir, xlenperpix, ylenperpix, gt, xsize, ysize);
         }
      }

      // lines
      for (int i = 0; i < nlines; ++i) {
         //NShapes::Line line(-h, w, linelength, thetal, thetar, radl, radr);
         NShapes::Line line(-lineParams[i]->h, lineParams[i]->w, lineParams[i]->linelength, lineParams[i]->thetal, lineParams[i]->thetar, lineParams[i]->radl, lineParams[i]->radr);
         //NShapes::Line line(-lp.h, lp.w, lp.linelength, lp.thetal, lp.thetar, lp.radl, lp.radr);
         const double pivot[3] = { 0.f, 0.f, 0.f };
         line.get()->rotate(pivot, -Math2::PI / 2.f, Math2::PI / 2.f, Math2::PI / 2.f);
         //line.get()->rotate(pivot, -Math2::PI / 2.f + 20.f / (Math2::PI / 2.f), Math2::PI / 2.f - 20.f / (Math2::PI / 2.f), Math2::PI / 2.f - 20.f / (Math2::PI / 2.f));
         //const double dist1[3] = { 0.f, 0.f, linelength / 2. };
         //const double dist1[3] = { lp.x, 0.f, linelength / 2. };
         const double offset[3] = { lineParams[i]->x, 0.f, linelength / 2. };
         line.get()->translate(offset);

         line.calcGroundtruth(); // get points/line segments that need to be projected

         // boundary corrections
         if (lineParams[i]->material == hstripParams[0]->material) { // bottom needs to be removed due to same material
            line.calcRasterizationCorrection(projectionPlane, axis0, axis1, beamdir, xlenperpix, ylenperpix, gt, xsize, ysize);
         }
      }

      //NShapes::TestProjection();
   }

   void writeSerializedParams(const char fname[])
   {
      std::fstream fout(fname, std::ios::out | std::ios::app);

      fout << xstartnm << ",";
      fout << xstopnm << ",";

      fout << ystartnm << ",";
      fout << ystopnm << ",";

      fout << xsize << ",";
      fout << ysize << ",";

      fout << nhstrips << ",";
      for (int i = 0; i < nhstrips; ++i) {
         fout << hstripParams[i]->w << ",";
         fout << hstripParams[i]->fadetop << ",";
         fout << hstripParams[i]->fadebot << ",";

         fout << hstripParams[i]->material << ",";
         fout << hstripParams[i]->y << ",";
         fout << hstripParams[i]->z << ",";
      }

      fout << nlines << ",";
      for (int i = 0; i < nlines; ++i) {
         fout << lineParams[i]->h << ",";
         fout << lineParams[i]->w << ",";
         fout << lineParams[i]->linelength << ",";
         fout << lineParams[i]->thetal << ",";
         fout << lineParams[i]->thetar << ",";
         fout << lineParams[i]->radl << ",";
         fout << lineParams[i]->radr << ",";

         fout << lineParams[i]->material << ",";
         fout << lineParams[i]->x << ",";
      }

      fout << nTrajectories << ",";
      fout << beamEeV << ",";
      fout << beamsizenm << ",";
      fout << beamz << ",";
      fout << beamphideg << ",";
      fout << beamthetadeg << ",";
      fout << beamfocallength << ",";

      fout << "\n";

      fout.close();
   }

   void readSerializedParams(const char fname[])
   {
      std::fstream fin(fname, std::ios::in);
      std::string line, word;
      while (fin >> line) {
         std::stringstream s(line);

         getline(s, word, ',');
         xstartnm = std::stof(word.c_str()); //printf("%.5e\n", xstartnm);

         getline(s, word, ',');
         xstopnm = std::stof(word.c_str()); //printf("%.5e\n", xstopnm);

         getline(s, word, ',');
         ystartnm = std::stof(word.c_str()); //printf("%.5e\n", ystartnm);

         getline(s, word, ',');
         ystopnm = std::stof(word.c_str()); //printf("%.5e\n", ystopnm);

         getline(s, word, ',');
         xsize = std::stoi(word.c_str()); //printf("%d\n", xsize);

         getline(s, word, ',');
         ysize = std::stoi(word.c_str()); //printf("%d\n", ysize);

         getline(s, word, ',');
         nhstrips = std::stoi(word.c_str()); //printf("%d\n", nhstrips);
         hstripParams = new NShapes::HorizontalStripParams*[nhstrips];
         for (int i = 0; i < nhstrips; ++i) {
            getline(s, word, ',');
            const float w = std::stof(word.c_str()); //printf("%.5e\n", w);
            getline(s, word, ',');
            const int ft = std::stoi(word.c_str()); //printf("%d\n", ft);
            getline(s, word, ',');
            const int fb = std::stoi(word.c_str()); //printf("%d\n", fb);

            getline(s, word, ',');
            const unsigned int mat = std::stoi(word.c_str()); //printf("%d\n", mat);
            getline(s, word, ',');
            const float y = std::stof(word.c_str()); //printf("%.5e\n", y);
            getline(s, word, ',');
            const float z = std::stof(word.c_str()); //printf("%.5e\n", z);

            hstripParams[i] = new NShapes::HorizontalStripParams(w, ft, fb, mat, y, z);
         }

         getline(s, word, ',');
         nlines = std::stoi(word.c_str()); //printf("%d\n", nlines);
         lineParams = new NShapes::LineParams*[nlines];
         for (int i = 0; i < nlines; ++i) {
            getline(s, word, ',');
            const float h = std::stof(word.c_str()); //printf("%.5e\n", h);
            getline(s, word, ',');
            const float w = std::stof(word.c_str()); //printf("%.5e\n", w);
            getline(s, word, ',');
            const float l = std::stof(word.c_str()); //printf("%.5e\n", l);
            getline(s, word, ',');
            const float tl = std::stof(word.c_str()); //printf("%.5e\n", tl);
            getline(s, word, ',');
            const float tr = std::stof(word.c_str()); //printf("%.5e\n", tr);
            getline(s, word, ',');
            const float rl = std::stof(word.c_str()); //printf("%.5e\n", rl);
            getline(s, word, ',');
            const float rr = std::stof(word.c_str()); //printf("%.5e\n", rr);

            getline(s, word, ',');
            const unsigned int mat = std::stoi(word.c_str()); //printf("%d\n", mat);
            getline(s, word, ',');
            const float x = std::stof(word.c_str()); //printf("%.5e\n", x);

            lineParams[i] = new NShapes::LineParams(h, w, l, tl, tr, rl, rr, mat, x);
         }

         getline(s, word, ',');
         nTrajectories = std::stoi(word.c_str()); //printf("%d\n", nTrajectories);
         getline(s, word, ',');
         beamEeV = std::stof(word.c_str()); //printf("%.5e\n", beamEeV);
         getline(s, word, ',');
         beamsizenm = std::stof(word.c_str()); //printf("%.5e\n", beamsizenm);
         getline(s, word, ',');
         beamz = std::stof(word.c_str()); //printf("%.5e\n", beamz);
         getline(s, word, ',');
         beamphideg = std::stof(word.c_str()); //printf("%.5e\n", beamphideg);
         getline(s, word, ',');
         beamthetadeg = std::stof(word.c_str()); //printf("%.5e\n", beamthetadeg);
         getline(s, word, ',');
         beamfocallength = std::stof(word.c_str()); //printf("%.5e\n", beamfocallength);
      }

      fin.close();
   }
}