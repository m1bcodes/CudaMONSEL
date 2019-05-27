#include "gov\nist\nanoscalemetrology\JMONSELTests\LinesOnLayers.cuh"

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

#include <fstream>

#include <chrono>

namespace LinesOnLayers
{
   static const char DefaultOutput[] = "results//LinesOnLayers";
   static const char PathSep[] = "//";

   void run()
   {
      //String dest = DefaultOutput;
      //new File(dest).mkdirs();
      //String filename = DefaultOutput + PathSep + "output.txt";

      StringT output = "Output will be to file: ";

      int seed = rand();
      srand(seed);
      output += "\n Random number seed: ";
      output += std::to_string(seed);
      for (int i = 0; i < 10; ++i) {
         double r = Math2::random();
         output += "\n " + std::to_string(r);
      }

      int nTrajectories = 100;

      double pitchnm = 180;
      int nlines = 3;
      double hnm = 120;
      double wnm = 80;
      double linelengthnm = 1000;
      double thetardeg = 3;
      double thetaldeg = 3;
      double radrnm = 20;
      double radlnm = 20;
      double layer1thicknessnm = 80;
      double layer2thicknessnm = 200;

      double beamEeVvals[] = { 500. };
      int beamEeVvalsLen = 1;
      double beamsizenm = 0.5;
      double deepnm = 15;

      bool trajImg = true;
      int trajImgMaxTraj = 50;
      double trajImgSize = 200e-9;

      bool VRML = false;
      int VRMLImgMaxTraj = 0;

      SEmaterialT vacuum;
      vacuum.setName("SE vacuum");
      ExpQMBarrierSMT vacuumBarrier(&vacuum);

      MONSEL_MaterialScatterModelT vacuumMSM(&vacuum, &vacuumBarrier);

      double breakE = ToSI::eV(45.);
      double density = 1190.;
      double workfun = 5.5;
      double bandgap = 5.;
      double EFermi = -bandgap;
      double potU = -workfun - EFermi;

      const ElementT& C = Element::C;
      const ElementT& Ox = Element::O;
      const ElementT& H = Element::H;
      CompositionT PMMAcomp;
      const ElementT* componentsCOH[] = { &C, &Ox, &H };
      const double compositionCOH[] = { 5, 2, 8 };
      PMMAcomp.defineByMoleFraction(componentsCOH, 3, compositionCOH, 3);
      SEmaterialT PMMA(PMMAcomp, density);
      PMMA.setName("PMMA");
      double PMMAWorkfunction = ToSI::eV(workfun);
      PMMA.setWorkfunction(PMMAWorkfunction);
      PMMA.setBandgap(ToSI::eV(bandgap));
      PMMA.setEnergyCBbottom(ToSI::eV(potU));

      SelectableElasticSMT PMMANISTMott(PMMA, NISTMottRS::Factory);
      JoyLuoNieminenCSDT PMMACSD(PMMA, breakE);
      FittedInelSMT PMMAfittedInel(PMMA, ToSI::eV(65.4), PMMACSD);
      GanachaudMokraniPolaronTrapSMT PMMApolaron(2.e7, 1. / ToSI::eV(4.));

      ExpQMBarrierSMT pmmaeqmbsm(&PMMA);

      MONSEL_MaterialScatterModelT PMMAMSM(&PMMA, &pmmaeqmbsm);
      PMMAMSM.addScatterMechanism(&PMMANISTMott);
      PMMAMSM.addScatterMechanism(&PMMAfittedInel);
      PMMAMSM.addScatterMechanism(&PMMApolaron);

      PMMAMSM.setCSD(&PMMACSD);

      MONSEL_MaterialScatterModelT PMMAMSMDeep(&PMMA, &pmmaeqmbsm);
      PMMAMSMDeep.addScatterMechanism(&PMMANISTMott);
      PMMAMSMDeep.addScatterMechanism(&PMMAfittedInel);
      PMMAMSMDeep.addScatterMechanism(&PMMApolaron);

      PMMAMSMDeep.setCSD(&PMMACSD);
      PMMAMSMDeep.setMinEforTracking(ToSI::eV(50.));

      MONSEL_MaterialScatterModelT ARCMSM = PMMAMSM;

      density = 1800.;
      workfun = 5.0;
      bandgap = 0.;
      EFermi = 20.4;
      potU = -workfun - EFermi;
      const ElementT* glCComponents[] = { &Element::C };
      const double glCComposition[] = { 1. };
      SEmaterialT glC(glCComponents, 1, glCComposition, 1, density, "glassy Carbon");
      double glCWorkfunction = ToSI::eV(workfun);
      glC.setWorkfunction(glCWorkfunction);
      glC.setEnergyCBbottom(ToSI::eV(potU));
      glC.setBandgap(ToSI::eV(bandgap));
      const double glCCoreEnergy[] = { ToSI::eV(284.2) };
      glC.setCoreEnergy(glCCoreEnergy, 1);

      StringT tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\";
      StringT glCTables[] = {
         tablePath + "IIMFPPennInterpglassyCSI.csv",
         tablePath + "interpNUSimReducedDeltaEglassyCSI.csv",
         tablePath + "interpsimTableThetaNUglassyCSI.csv",
         tablePath + "interpSimESE0NUglassyCSI.csv"
      };

      SelectableElasticSMT glCNISTMott(glC, NISTMottRS::Factory);
      TabulatedInelasticSMT glCDS(glC, 3, glCTables);

      ExpQMBarrierSMT glceqmbsm(&glC);

      MONSEL_MaterialScatterModelT glCMSM(&glC, &glceqmbsm);
      glCMSM.addScatterMechanism(&glCNISTMott);
      glCMSM.addScatterMechanism(&glCDS);

      MONSEL_MaterialScatterModelT glCMSMDeep(&glC, &glceqmbsm);
      glCMSMDeep.addScatterMechanism(&glCNISTMott);
      glCMSMDeep.addScatterMechanism(&glCDS);

      glCMSMDeep.setMinEforTracking(ToSI::eV(50.));

      double phononE = 0.063;
      double phononStrength = 3.;

      density = 2330.;
      workfun = 4.85;
      bandgap = 1.1;
      EFermi = -bandgap;
      potU = -workfun - EFermi;
      const ElementT* SiComponent[] = { &Element::Si };
      const double SiComposition[] = { 1. };
      SEmaterialT Si(SiComponent, 1, SiComposition, 1, density, "Silicon");
      double SiWorkfunction = ToSI::eV(workfun);
      Si.setWorkfunction(SiWorkfunction);
      Si.setEnergyCBbottom(ToSI::eV(potU));
      Si.setBandgap(ToSI::eV(bandgap));
      double SiCoreEnergy[] = { ToSI::eV(99.2), ToSI::eV(99.8), ToSI::eV(149.7), ToSI::eV(1839.) };
      Si.setCoreEnergy(SiCoreEnergy, 4);

      tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\";
      StringT SiTables[] = {
         tablePath + "IIMFPFullPennInterpSiSI.csv",
         tablePath + "interpNUSimReducedDeltaEFullPennSiSI.csv",
         tablePath + "interpNUThetaFullPennSiBGSI.csv",
         tablePath + "interpSimESE0NUSiBGSI.csv"
      };

      SelectableElasticSMT SiNISTMott(Si, NISTMottRS::Factory);
      TabulatedInelasticSMT SiDS(Si, 3, SiTables, ToSI::eV(13.54));

      GanachaudMokraniPhononInelasticSMT Siphonon(phononStrength, ToSI::eV(phononE), 300., 11.7, 1.);

      ExpQMBarrierSMT sieqmbsm(&Si);

      MONSEL_MaterialScatterModelT SiMSM(&Si, &sieqmbsm);
      SiMSM.addScatterMechanism(&SiNISTMott);
      SiMSM.addScatterMechanism(&SiDS);
      SiMSM.addScatterMechanism(&Siphonon);

      MONSEL_MaterialScatterModelT SiMSMDeep(&Si, &sieqmbsm);
      SiMSMDeep.addScatterMechanism(&SiNISTMott);
      SiMSMDeep.addScatterMechanism(&SiDS);
      SiMSMDeep.addScatterMechanism(&Siphonon);

      SiMSMDeep.setMinEforTracking(ToSI::eV(50.));

      double meterspernm = 1.e-9;
      double pitch = pitchnm*meterspernm;
      double h = hnm*meterspernm;
      double w = wnm*meterspernm;
      double linelength = linelengthnm*meterspernm;

      double radperdeg = Math2::PI / 180.;
      double thetar = thetardeg*radperdeg;
      double thetal = thetaldeg*radperdeg;
      double radr = radrnm*meterspernm;
      double radl = radlnm*meterspernm;
      double layer1thickness = layer1thicknessnm*meterspernm;
      double layer2thickness = layer2thicknessnm*meterspernm;
      double beamsize = beamsizenm*meterspernm;
      double deep = deepnm*meterspernm;

      NullMaterialScatterModelT NULL_MSM;
      const double center[] = {
         0.0,
         0.0,
         0.0
      };

      double beamEeV = beamEeVvals[0];
      double beamE = ToSI::eV(beamEeV);
      SphereT sphere(center, MonteCarloSS::ChamberRadius);
      GaussianBeamT eg(beamsize, beamE, center);

      RegionT chamber(nullptr, &NULL_MSM, &sphere);
      chamber.updateMaterial(*chamber.getScatterModel(), vacuumMSM);

      const double normalvector[] = { 0., 0., -1. };
      const double layer1Pos[] = { 0., 0., 0. };
      NormalMultiPlaneShapeT layer1;
      PlaneT pl1(normalvector, 3, layer1Pos, 3);
      layer1.addPlane(pl1);
      RegionT layer1Region(&chamber, &ARCMSM, (NormalShapeT*)&layer1);

      double layer2Pos[] = { 0., 0., layer1thickness };
      NormalMultiPlaneShapeT layer2;
      PlaneT pl2(normalvector, 3, layer2Pos, 3);
      layer2.addPlane(pl2);
      RegionT layer2Region(&layer1Region, &glCMSM, (NormalShapeT*)&layer2);

      double layer3Pos[] = { 0., 0., layer1thickness + layer2thickness };
      NormalMultiPlaneShapeT layer3;
      PlaneT pl3(normalvector, 3, layer2Pos, 3);
      layer3.addPlane(pl3);
      RegionT layer3Region(&layer2Region, &SiMSM, (NormalShapeT*)&layer3);

      double layer4Pos[] = { 0., 0., layer1thickness + layer2thickness + deep };
      NormalMultiPlaneShapeT layer4;
      PlaneT pl4(normalvector, 3, layer4Pos, 3);
      RegionT layer4Region(&layer3Region, &SiMSM, (NormalShapeT*)&layer4);

      RegionT deepRegion(&layer3Region, &SiMSMDeep, (NormalShapeT*)&layer4);

      double leftmostLineCenterx = -pitch*(nlines / 2);
      //for (int i = 0; i < nlines; ++i) {
      //   double xcenter = leftmostLineCenterx + i*pitch;
      //   NormalIntersectionShapeT* line = (NormalIntersectionShapeT*)NShapes::createLine(-h, w, linelength, thetal, thetar, radl, radr);
      //   const double newLinePos[] = { xcenter, 0., 0. };
      //   line->translate(newLinePos);
      //   RegionT lineRegion(&chamber, &PMMAMSM, line);
      //}

      double xcenter = leftmostLineCenterx + 0*pitch;
      NormalIntersectionShapeT* line = (NormalIntersectionShapeT*)NShapes::createLine(-h, w, linelength, thetal, thetar, radl, radr);
      RegionT lineRegion(&chamber, &PMMAMSM, line);

//      VectorXd yvals = { 0. };
      VectorXd yvals;
      for (int i = -300; i < 100; ++i) {
         yvals.push_back(i);
      }
      for (int i = 100; i < 200; ++i) {
         yvals.push_back(100 + i*5);
      }

      double xbottom = wnm / 2.;
      double xtop = wnm / 2. - hnm * ::tan(thetar);
      double xstart = xbottom - 300.5;
      double xstop = xbottom + 100.5;
      double xfinestart = xtop - 20.5;
      double xfinestop;
      if (thetar < 0.) {
         xfinestop = xtop + 20.5;
      }
      else {
         xfinestop = wnm / 2. + 20.5;
      }

      VectorXd xvals;
      double deltax = 5.;
      double x = xstart;
      while (x < xfinestart) {
         xvals.push_back(x);
         x += deltax;
      }
      x = xfinestart;
      deltax = 1.;
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

      double binSizeEV = 10.;

      MonteCarloSS::MonteCarloSS monte(&eg, &chamber, eg.createElectron());

      output += "\n# Trajectories at each landing position: " + std::to_string(nTrajectories);
      output += "\n# Pitch of lines (nm): " + std::to_string(pitchnm);
      output += "\n# lines: " + std::to_string(nlines);
      output += "\nLine height (nm): " + std::to_string(hnm);
      output += "\nLine bottom width (nm): " + std::to_string(wnm);
      output += "\nLine length (nm): " + std::to_string(linelengthnm);
      output += "\nLeft and right sidewall angles (deg): " + std::to_string(thetaldeg) + " " + std::to_string(thetardeg);
      output += "\nLeft and right top corner radii (nm): " + std::to_string(radlnm) + " " + std::to_string(radrnm);
      output += "\nThicknesses of 1st and second layers (nm): " + std::to_string(layer1thicknessnm) + " " + std::to_string(layer2thicknessnm);
      output += "\nBeam landing energies (eV): ";

      for (int i = 0; i < beamEeVvalsLen; i++) {
         output += std::to_string(beamEeVvals[i]);
      }
      output += "\nBeam size (standard deviation, in nm): " + std::to_string(beamsizenm);

      output += "\n";
      output += "\nbeamE (eV)	x(nm)	y (nm)	BSE yield	SE yield";

      auto start = std::chrono::system_clock::now();
      for (auto ynm : yvals) {
         //double ynm = yvals[0];
         double y = ynm*meterspernm;
         for (auto xnm : xvals) {
            x = xnm*meterspernm;
            double egCenter[] = { x, y, -h - 20.*meterspernm };
            eg.setCenter(egCenter);

            int nbins = (int)(beamEeV / binSizeEV);
            BackscatterStatsT back(monte, nbins);
            monte.addActionListener(back);

            try {
               monte.runMultipleTrajectories(nTrajectories);

               const HistogramT& hist = back.backscatterEnergyHistogram();

               double energyperbineV = beamEeV / hist.binCount();
               double maxSEbin = 50. / energyperbineV;
               int totalSE = 0;
               for (int j = 0; j < (int)maxSEbin; ++j) {
                  totalSE = totalSE + hist.counts(j);
               }

               double SEf = (float)totalSE / nTrajectories;
               double bsf = back.backscatterFraction() - SEf;
               output += "\n" + std::to_string(beamEeV) + " " + std::to_string(xnm) + " " + std::to_string(ynm) + " " + std::to_string(bsf) + " " + std::to_string(SEf);

               StringT tmp("\n" + std::to_string(beamEeV) + " " + std::to_string(xnm) + " " + std::to_string(ynm) + " " + std::to_string(bsf) + " " + std::to_string(SEf));
               printf("%s\n", tmp.c_str());

               monte.removeActionListener(back);
            }
            catch (std::exception&) {
               printf("wtfweewew\n");
            }
         }
      }
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;
      std::time_t end_time = std::chrono::system_clock::to_time_t(end);
      std::cout << "finished computation at " << std::ctime(&end_time) << "elapsed time: " << elapsed_seconds.count() << "s\n";
      output += "\n" + std::to_string(elapsed_seconds.count());

      std::ofstream myfile;
      myfile.open("output.txt");
      myfile << output.c_str();
      myfile.close();
   }
}
