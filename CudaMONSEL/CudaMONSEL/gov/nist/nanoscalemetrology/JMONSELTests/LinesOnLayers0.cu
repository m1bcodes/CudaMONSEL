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
   void copyDataToCuda()
   {
      NUTableInterpolation::initFactory << <1, 1 >> >();
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());

      StringT tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\";
      const StringT IIMFPPennInterpglassy = tablePath + "IIMFPPennInterpglassyCSI.csv";
      const StringT SimReducedDeltaEglassy = tablePath + "interpNUSimReducedDeltaEglassyCSI.csv";
      const StringT simTableThetaNUglassy = tablePath + "interpsimTableThetaNUglassyCSI.csv";
      const StringT SimESE0NUglassy = tablePath + "interpSimESE0NUglassyCSI.csv";
      NUTableInterpolation::copyDataToCuda(IIMFPPennInterpglassy.c_str());
      NUTableInterpolation::copyDataToCuda(SimReducedDeltaEglassy.c_str());
      NUTableInterpolation::copyDataToCuda(simTableThetaNUglassy.c_str());
      NUTableInterpolation::copyDataToCuda(SimESE0NUglassy.c_str());

      tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\";
      const StringT IIMFPFullPennInterpSiSI = tablePath + "IIMFPFullPennInterpSiSI.csv";
      const StringT interpNUSimReducedDeltaEFullPennSiSI = tablePath + "interpNUSimReducedDeltaEFullPennSiSI.csv";
      const StringT interpNUThetaFullPennSiBGSI = tablePath + "interpNUThetaFullPennSiBGSI.csv";
      const StringT interpSimESE0NUSiBGSI = tablePath + "interpSimESE0NUSiBGSI.csv";
      NUTableInterpolation::copyDataToCuda(IIMFPFullPennInterpSiSI.c_str());
      NUTableInterpolation::copyDataToCuda(interpNUSimReducedDeltaEFullPennSiSI.c_str());
      NUTableInterpolation::copyDataToCuda(interpNUThetaFullPennSiBGSI.c_str());
      NUTableInterpolation::copyDataToCuda(interpSimESE0NUSiBGSI.c_str());
   }

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

   __device__ SEmaterialT* vacuum = nullptr;
   __device__ ExpQMBarrierSMT* vacuumBarrier = nullptr;
   __device__ ZeroCSDT* sZeroCSD = nullptr;

   __device__ MONSEL_MaterialScatterModelT* vacuumMSM = nullptr;

   __constant__ const double PMMAbreakE = 1.60217653e-19 * 45.;
   __constant__ const double PMMAdensity = 1190.;
   __constant__ const double PMMAworkfun = 5.5;
   __constant__ const double PMMAbandgap = 5.;
   __constant__ const double PMMAEFermi = -5.;//-PMMAbandgap;
   __constant__ const double PMMApotU = -5.5 - (-5.);

   __device__ SEmaterialT* PMMA = nullptr;

   __device__ SelectableElasticSMT* PMMANISTMott = nullptr;

   __device__ JoyLuoNieminenCSDT* PMMACSD = nullptr;
   __device__ FittedInelSMT* PMMAfittedInel = nullptr;
   __device__ GanachaudMokraniPolaronTrapSMT* PMMApolaron = nullptr;

   __device__ ExpQMBarrierSMT* pmmaeqmbsm = nullptr;

   __device__ MONSEL_MaterialScatterModelT* PMMAMSM = nullptr;

   __device__ MONSEL_MaterialScatterModelT* PMMAMSMDeep = nullptr;

   __device__ MONSEL_MaterialScatterModelT* ARCMSM = nullptr;

   __constant__ const double glCdensity = 1800.;
   __constant__ const double glCworkfun = 5.0;
   __constant__ const double glCbandgap = 0.;
   __constant__ const double glCEFermi = 20.4;
   __constant__ const double glCpotU = -5. - 20.4;

   __device__ SEmaterialT* glC = nullptr;

   __device__ SelectableElasticSMT* glCNISTMott = nullptr;

   __device__ TabulatedInelasticSMT* glCDS = nullptr;

   __device__ ExpQMBarrierSMT* glceqmbsm = nullptr;

   __device__ MONSEL_MaterialScatterModelT* glCMSM = nullptr;

   __device__ MONSEL_MaterialScatterModelT* glCMSMDeep = nullptr;

   __constant__ const double phononE = 0.063;
   __constant__ const double phononStrength = 3.;

   __constant__ const double Sidensity = 2330.;
   __constant__ const double Siworkfun = 4.85;
   __constant__ const double Sibandgap = 1.1;
   __constant__ const double SiEFermi = -1.1;//-Sibandgap;
   __constant__ const double SipotU = -44.85 - (-1.1);//-Siworkfun - SiEFermi;

   __device__ SEmaterialT* Si = nullptr;

   __device__ SelectableElasticSMT* SiNISTMott = nullptr;

   __device__ TabulatedInelasticSMT* SiDS = nullptr;

   __device__ GanachaudMokraniPhononInelasticSMT* Siphonon = nullptr;

   __device__ ExpQMBarrierSMT* sieqmbsm = nullptr;

   __device__ MONSEL_MaterialScatterModelT* SiMSM = nullptr;

   __device__ MONSEL_MaterialScatterModelT* SiMSMDeep = nullptr;

   __device__ SphereT* sphere = nullptr;
   __device__ GaussianBeamT* eg = nullptr;

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

   __device__ NullMaterialScatterModelT* NULL_MSM = nullptr;

   __device__ RegionT* chamber = nullptr;

   __constant__ const double normalvector[] = { 0., 0., -1. };
   __constant__ const double layer1Pos[] = { 0., 0., 0. };

   __device__ NormalMultiPlaneShapeT* layer1 = nullptr;
   __device__ PlaneT* pl1 = nullptr;
   __device__ RegionT* layer1Region = nullptr;

   __constant__ const double layer2Pos[] = { 0., 0., 80 * 1.e-9 };
   __device__ NormalMultiPlaneShapeT* layer2 = nullptr;
   __device__ PlaneT* pl2 = nullptr;
   __device__ RegionT* layer2Region = nullptr;

   __constant__ const double layer3Pos[] = { 0., 0., 80 * 1.e-9 + 200 * 1.e-9 };
   __device__ NormalMultiPlaneShapeT* layer3 = nullptr;
   __device__ PlaneT* pl3 = nullptr;
   __device__ RegionT* layer3Region = nullptr;

   __constant__ const double layer4Pos[] = { 0., 0., 80 * 1.e-9 + 200 * 1.e-9 + 15 * 1.e-9 };
   __device__ NormalMultiPlaneShapeT* layer4 = nullptr;
   __device__ PlaneT* pl4 = nullptr;
   __device__ RegionT* layer4Region = nullptr;

   __device__ RegionT* deepRegion = nullptr;

   __constant__ const double leftmostLineCenterx = -180. * 1.e-9 * (3. / 2.);
   __constant__ const double xcenter = -180. * 1.e-9 * (3. / 2.) + 0 * 180 * 1.e-9;

   __device__ NormalIntersectionShapeT* line = nullptr;
   __device__ RegionT* lineRegion = nullptr;

   __device__ MonteCarloSST* monte = nullptr;

   __global__ void initCuda()
   {
      printf("LinesOnLayers: initCuda");
   }

   __global__ void runCuda()
   {
      printf("LinesOnLayers: runCuda");
      for (int i = 0; i < 10; ++i) {
         printf("\n %.10e", Random::random());
      }

      vacuum = new SEmaterialT(); printf("0");
      vacuum->setName("SE vacuum");
      vacuumBarrier = new ExpQMBarrierSMT(vacuum); printf("1");
      sZeroCSD = new ZeroCSDT(); printf("2");

      vacuumMSM = new MONSEL_MaterialScatterModelT(vacuum, vacuumBarrier, sZeroCSD); printf("3");

      const ElementT* componentsCOH[] = { Element::dC, Element::dO, Element::dH };
      const double compositionCOH[] = { 5. / 15., 2. / 15., 8. / 15. };
      PMMA = new SEmaterialT(componentsCOH, 3, compositionCOH, 3, PMMAdensity, "PMMA"); printf("4");
      PMMA->setWorkfunction(ToSI::eV(PMMAworkfun));
      PMMA->setBandgap(ToSI::eV(PMMAbandgap));
      PMMA->setEnergyCBbottom(ToSI::eV(PMMApotU));

      PMMANISTMott = new SelectableElasticSMT(*PMMA, *NISTMottRS::d_Factory); printf("5");

      PMMACSD = new JoyLuoNieminenCSDT(*PMMA, PMMAbreakE); printf("6");
      PMMAfittedInel = new FittedInelSMT(*PMMA, ToSI::eV(65.4), *PMMACSD); printf("7");
      PMMApolaron = new GanachaudMokraniPolaronTrapSMT(2.e7, 1. / ToSI::eV(4.)); printf("8");

      pmmaeqmbsm = new ExpQMBarrierSMT(PMMA); printf("9");

      PMMAMSM = new MONSEL_MaterialScatterModelT(PMMA, pmmaeqmbsm, sZeroCSD); printf("10");
      PMMAMSM->addScatterMechanism(PMMANISTMott);
      PMMAMSM->addScatterMechanism(PMMAfittedInel);
      PMMAMSM->addScatterMechanism(PMMApolaron);

      PMMAMSM->setCSD(PMMACSD);

      PMMAMSMDeep = new MONSEL_MaterialScatterModelT(PMMA, pmmaeqmbsm, sZeroCSD); printf("11");
      PMMAMSMDeep->addScatterMechanism(PMMANISTMott);
      PMMAMSMDeep->addScatterMechanism(PMMAfittedInel);
      PMMAMSMDeep->addScatterMechanism(PMMApolaron);

      PMMAMSMDeep->setCSD(PMMACSD);
      PMMAMSMDeep->setMinEforTracking(ToSI::eV(50.));

      ARCMSM = PMMAMSM; printf("12");

      const ElementT* glCComponents[] = { Element::dC };
      const double glCComposition[] = { 1. };
      glC = new SEmaterialT(glCComponents, 1, glCComposition, 1, glCdensity, "glassy Carbon"); printf("13");
      glC->setWorkfunction(ToSI::eV(glCworkfun));
      glC->setEnergyCBbottom(ToSI::eV(glCpotU));
      glC->setBandgap(ToSI::eV(glCbandgap));
      const double glCCoreEnergy[] = { ToSI::eV(284.2) };
      glC->setCoreEnergy(glCCoreEnergy, 1);

      glCNISTMott = new SelectableElasticSMT(*glC, *NISTMottRS::d_Factory); printf("14");

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

      glCDS = new TabulatedInelasticSMT(*glC, 3, glCTables); printf("15");

      glceqmbsm = new ExpQMBarrierSMT(glC); printf("16");

      glCMSM = new MONSEL_MaterialScatterModelT(glC, glceqmbsm, sZeroCSD); printf("17");
      glCMSM->addScatterMechanism(glCNISTMott);
      glCMSM->addScatterMechanism(glCDS);

      glCMSMDeep = new MONSEL_MaterialScatterModelT(glC, glceqmbsm, sZeroCSD); printf("18");
      glCMSMDeep->addScatterMechanism(glCNISTMott);
      glCMSMDeep->addScatterMechanism(glCDS);

      glCMSMDeep->setMinEforTracking(ToSI::eV(50.));

      const ElementT* SiComponent[] = { Element::dSi };
      const double SiComposition[] = { 1. };
      Si = new SEmaterialT(SiComponent, 1, SiComposition, 1, Sidensity, "Silicon"); printf("19");
      Si->setWorkfunction(ToSI::eV(Siworkfun));
      Si->setEnergyCBbottom(ToSI::eV(SipotU));
      Si->setBandgap(ToSI::eV(Sibandgap));
      const double SiCoreEnergy[] = { ToSI::eV(99.2), ToSI::eV(99.8), ToSI::eV(149.7), ToSI::eV(1839.) };
      Si->setCoreEnergy(SiCoreEnergy, 4);

      SiNISTMott = new SelectableElasticSMT(*Si, *NISTMottRS::d_Factory); printf("20");

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

      SiDS = new TabulatedInelasticSMT(*Si, 3, SiTables, ToSI::eV(13.54)); printf("21");

      Siphonon = new GanachaudMokraniPhononInelasticSMT(phononStrength, ToSI::eV(phononE), 300., 11.7, 1.); printf("22");

      sieqmbsm = new ExpQMBarrierSMT(Si); printf("23");

      SiMSM = new MONSEL_MaterialScatterModelT(Si, sieqmbsm, sZeroCSD); printf("24");
      SiMSM->addScatterMechanism(SiNISTMott);
      SiMSM->addScatterMechanism(SiDS);
      SiMSM->addScatterMechanism(Siphonon);

      SiMSMDeep = new MONSEL_MaterialScatterModelT(Si, sieqmbsm, sZeroCSD); printf("25");
      SiMSMDeep->addScatterMechanism(SiNISTMott);
      SiMSMDeep->addScatterMechanism(SiDS);
      SiMSMDeep->addScatterMechanism(Siphonon);

      SiMSMDeep->setMinEforTracking(ToSI::eV(50.));

      const double beamEeV = beamEeVvals[0];
      const double beamE = ToSI::eV(beamEeV);
      sphere = new SphereT(center, MonteCarloSS::ChamberRadius); printf("26");
      eg = new GaussianBeamT(beamsize, beamE, center); printf("27");

      NULL_MSM = new NullMaterialScatterModelT(); printf("28");
      chamber = new RegionT(nullptr, NULL_MSM, sphere); printf("29");
      chamber->updateMaterial(*(chamber->getScatterModel()), *vacuumMSM);

      layer1 = new NormalMultiPlaneShapeT(); printf("30");
      pl1 = new PlaneT(normalvector, layer1Pos); printf("31");
      layer1->addPlane(pl1);
      layer1Region = new RegionT(chamber, ARCMSM, (NormalShapeT*)layer1); printf("32");

      layer2 = new NormalMultiPlaneShapeT(); printf("33");
      pl2 = new PlaneT(normalvector, layer2Pos); printf("34");
      layer2->addPlane(pl2);
      layer2Region = new RegionT(layer1Region, glCMSM, (NormalShapeT*)layer2); printf("35");

      layer3 = new NormalMultiPlaneShapeT(); printf("36");
      pl3 = new PlaneT(normalvector, layer2Pos); printf("37");
      layer3->addPlane(pl3);
      layer3Region = new RegionT(layer2Region, SiMSM, (NormalShapeT*)layer3); printf("38");

      layer4 = new NormalMultiPlaneShapeT(); printf("39");
      pl4 = new PlaneT(normalvector, layer4Pos); printf("40");
      layer4Region = new RegionT(layer3Region, SiMSM, (NormalShapeT*)layer4); printf("41");

      deepRegion = new RegionT(layer3Region, SiMSMDeep, (NormalShapeT*)layer4); printf("42");

      //for (int i = 0; i < nlines; ++i) {
      //   double xcenter = leftmostLineCenterx + i*pitch;
      //   NormalIntersectionShapeT* line = (NormalIntersectionShapeT*)NShapes::createLine(-h, w, linelength, thetal, thetar, radl, radr);
      //   const double newLinePos[] = { xcenter, 0., 0. };
      //   line->translate(newLinePos);
      //   RegionT lineRegion(&chamber, &PMMAMSM, line);
      //}

      line = (NormalIntersectionShapeT*)NShapes::createLine(-h, w, linelength, thetal, thetar, radl, radr); printf("43");
      lineRegion = new RegionT(chamber, PMMAMSM, line); printf("44");

      //VectorXd yvals = { 0. };
      VectorXd yvals(128); printf("45");
      //for (int i = -100; i < 100; ++i) {
      //   yvals.push_back(i);
      //}
      for (int i = -64; i < 64; i += 1) {
         yvals.push_back(i);
      }

      const double xbottom = wnm / 2.;
      const double xtop = wnm / 2. - hnm * ::tan(thetar);
      const double xstart = xbottom - 100.5;
      const double xstop = xbottom + 100.5;
      const double xfinestart = xtop - 20.5;
      double xfinestop;
      if (thetar < 0.) xfinestop = xtop + 20.5;
      else xfinestop = wnm / 2. + 20.5;

      VectorXd xvals(128); printf("46");
      double deltax = 5.;
      double x = xstart;
      while (x < xfinestart) {
         xvals.push_back(x);
         x += deltax;
      }
      x = xfinestart;
      deltax = 1;
      while (x < xfinestop) {
         xvals.push_back(x);
         x += deltax;
      }
      x = xfinestop;
      deltax = 5.;
      while (x < xstop) {
         xvals.push_back(x);
         x += deltax;
      }
      xvals.push_back(xstop);

      const double binSizeEV = 10.;

      monte = new MonteCarloSST(eg, chamber, eg->createElectron()); printf("47");

      printf("\n# Trajectories at each landing position: %d", nTrajectories);
      printf("\n# Pitch of lines (nm): %.10e", pitchnm);
      printf("\n# lines: %d", nlines);
      printf("\nLine height (nm): %.10e", hnm);
      printf("\nLine bottom width (nm): %.10e", wnm);
      printf("\nLine length (nm): %.10e", linelengthnm);
      printf("\nLeft and right sidewall angles (deg): %.10e %.10e", thetaldeg, thetardeg);
      printf("\nLeft and right top corner radii (nm): %.10e %.10e", radlnm, radrnm);
      printf("\nThicknesses of 1st and second layers (nm): %.10e %.10e", layer1thicknessnm, layer2thicknessnm);
      printf("\nBeam landing energies (eV): ");

      for (int i = 0; i < beamEeVvalsLen; i++) {
         printf("\n%.10e", beamEeVvals[i]);
      }
      printf("\nBeam size (standard deviation, in nm): %.10e", beamsizenm);

      printf("\n");
      printf("\nbeamE (eV)\t x(nm)\t y (nm)\t BSE yield\t SE yield");

      for (auto ynm : yvals) {
         //double ynm = yvals[0];
         const double y = ynm * 1.e-9;
         for (auto xnm : xvals) {
            x = xnm*1.e-9;
            const double egCenter[] = { x, y, -h - 20.*1.e-9 };
            eg->setCenter(egCenter);

            const int nbins = (int)(beamEeV / binSizeEV);
            BackscatterStatsT* back = new BackscatterStatsT(*monte, nbins); printf("48");
            monte->addActionListener(*back);

            monte->runMultipleTrajectories(nTrajectories);

            const HistogramT& hist = back->backscatterEnergyHistogram(); printf("49");

            const double energyperbineV = beamEeV / hist.binCount();
            const double maxSEbin = 50. / energyperbineV;
            int totalSE = 0;
            for (int j = 0; j < (int)maxSEbin; ++j) {
               totalSE = totalSE + hist.counts(j);
            }

            const double SEf = (float)totalSE / nTrajectories;
            const double bsf = back->backscatterFraction() - SEf;
            printf("\n %.10e %.10e %.10e %.10e %.10e", beamEeV, xnm, ynm, bsf, SEf);

            monte->removeActionListener(*back);
            delete back;
         }
      }
      printf("\n");
   }
}