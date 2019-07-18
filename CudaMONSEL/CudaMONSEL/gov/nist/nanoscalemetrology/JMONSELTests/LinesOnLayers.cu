//#include "gov\nist\nanoscalemetrology\JMONSELTests\LinesOnLayers.cuh"
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
//
//#include "CudaUtil.h"
//
//namespace LinesOnLayers
//{
//   static const char DefaultOutput[] = "results//LinesOnLayers";
//   static const char PathSep[] = "//";
//
//   void run()
//   {
//      //String dest = DefaultOutput;
//      //new File(dest).mkdirs();
//      //String filename = DefaultOutput + PathSep + "output.txt";
//
//      std::string output = "Output will be to file: ";
//
//      const int seed = rand();
//      srand(seed);
//      output += "\n Random number seed: ";
//      output += std::to_string(seed);
//      for (int i = 0; i < 10; ++i) {
//         output += "\n " + std::to_string(Random::random());
//      }
//
//      const int nTrajectories = 100;
//
//      const double pitchnm = 180;
//      const int nlines = 3;
//      const double hnm = 120;
//      const double wnm = 80;
//      const double linelengthnm = 1000;
//      const double thetardeg = 3;
//      const double thetaldeg = 3;
//      const double radrnm = 20;
//      const double radlnm = 20;
//      const double layer1thicknessnm = 80;
//      const double layer2thicknessnm = 200;
//
//      const double beamEeVvals[] = { 500. };
//      const int beamEeVvalsLen = 1;
//      const double beamsizenm = 0.5;
//      const double deepnm = 15;
//
//      const bool trajImg = true;
//      const int trajImgMaxTraj = 50;
//      const double trajImgSize = 200e-9;
//
//      const bool VRML = false;
//      const int VRMLImgMaxTraj = 0;
//
//      SEmaterialT vacuum;
//      vacuum.setName("SE vacuum");
//      ExpQMBarrierSMT vacuumBarrier(&vacuum);
//      ZeroCSDT sZeroCSD;
//
//      MONSEL_MaterialScatterModelT vacuumMSM(&vacuum, &vacuumBarrier, &sZeroCSD);
//
//      const double breakE = ToSI::eV(45.);
//      double density = 1190.;
//      double workfun = 5.5;
//      double bandgap = 5.;
//      double EFermi = -bandgap;
//      double potU = -workfun - EFermi;
//
//      //const ElementT* componentsCOH[] = { &Element::C, &Element::O, &Element::H };
//      //const double compositionCOH[] = { 5. / 15., 2. / 15., 8. / 15. };
//      //SEmaterialT PMMA(componentsCOH, 3, compositionCOH, 3, density, "PMMA");
//      //PMMA.setWorkfunction(ToSI::eV(workfun));
//      //PMMA.setBandgap(ToSI::eV(bandgap));
//      //PMMA.setEnergyCBbottom(ToSI::eV(potU));
//
//      const ElementT& C = Element::C;
//      const ElementT& Ox = Element::O;
//      const ElementT& H = Element::H;
//      CompositionT PMMAcomp;
//      const ElementT* componentsCOH[] = { &C, &Ox, &H };
//      const double compositionCOH[] = { 5, 2, 8 };
//      PMMAcomp.defineByMoleFraction(componentsCOH, 3, compositionCOH, 3);
//      SEmaterialT PMMA(PMMAcomp, density);
//      PMMA.setName("PMMA");
//      PMMA.setWorkfunction(ToSI::eV(workfun));
//      PMMA.setBandgap(ToSI::eV(bandgap));
//      PMMA.setEnergyCBbottom(ToSI::eV(potU));
//
//      SelectableElasticSMT PMMANISTMott(PMMA, NISTMottRS::Factory);
//      JoyLuoNieminenCSDT PMMACSD(PMMA, breakE);
//      FittedInelSMT PMMAfittedInel(PMMA, ToSI::eV(65.4), PMMACSD);
//      GanachaudMokraniPolaronTrapSMT PMMApolaron(2.e7, 1. / ToSI::eV(4.));
//
//      ExpQMBarrierSMT pmmaeqmbsm(&PMMA);
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
//      PMMAMSMDeep.setMinEforTracking(ToSI::eV(50.));
//
//      MONSEL_MaterialScatterModelT& ARCMSM = PMMAMSM;
//
//      density = 1800.;
//      workfun = 5.0;
//      bandgap = 0.;
//      EFermi = 20.4;
//      potU = -workfun - EFermi;
//      const ElementT* glCComponents[] = { &Element::C };
//      const double glCComposition[] = { 1. };
//      SEmaterialT glC(glCComponents, 1, glCComposition, 1, density, "glassy Carbon");
//      double glCWorkfunction = ToSI::eV(workfun);
//      glC.setWorkfunction(glCWorkfunction);
//      glC.setEnergyCBbottom(ToSI::eV(potU));
//      glC.setBandgap(ToSI::eV(bandgap));
//      const double glCCoreEnergy[] = { ToSI::eV(284.2) };
//      glC.setCoreEnergy(glCCoreEnergy, 1);
//
//      StringT tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\";
//      const StringT IIMFPPennInterpglassy = tablePath + "IIMFPPennInterpglassyCSI.csv";
//      const StringT SimReducedDeltaEglassy = tablePath + "interpNUSimReducedDeltaEglassyCSI.csv";
//      const StringT simTableThetaNUglassy = tablePath + "interpsimTableThetaNUglassyCSI.csv";
//      const StringT SimESE0NUglassy = tablePath + "interpSimESE0NUglassyCSI.csv";
//      const char* glCTables[] = {
//         IIMFPPennInterpglassy.c_str(),
//         SimReducedDeltaEglassy.c_str(),
//         simTableThetaNUglassy.c_str(),
//         SimESE0NUglassy.c_str()
//      };
//
//      SelectableElasticSMT glCNISTMott(glC, NISTMottRS::Factory);
//      TabulatedInelasticSMT glCDS(glC, 3, glCTables);
//      printf("%d\n", glCDS.gettableIIMFP()->gettable1d().size());
//      for (auto &i : glCDS.gettableIIMFP()->gettable1d()) {
//         printf("%.10e ", i);
//      }
//      printf("\n");
//
//      printf("%d, %d\n", glCDS.gettableReducedDeltaE()->gettable2d().size(), glCDS.gettableReducedDeltaE()->gettable2d()[0].size());
//      for (auto &i : glCDS.gettableReducedDeltaE()->gettable2d()[0]) {
//         printf("%.10e ", i);
//      }
//      printf("\n");
//
//      printf("%d, %d, %d\n", glCDS.gettableTheta()->gettable3d().size(), glCDS.gettableTheta()->gettable3d()[50].size(), glCDS.gettableTheta()->gettable3d()[50][50].size());
//      for (auto &i : glCDS.gettableTheta()->gettable3d()[50][50]) {
//         printf("%.10e ", i);
//      }
//      printf("\n");
//
//      printf("%d, %d\n", glCDS.gettableSEE0()->gettable2d().size(), glCDS.gettableSEE0()->gettable2d()[0].size());
//      for (auto &i : glCDS.gettableSEE0()->gettable2d()[0]) {
//         printf("%.10e ", i);
//      }
//      printf("\n");
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
//      glCMSMDeep.setMinEforTracking(ToSI::eV(50.));
//
//      const double phononE = 0.063;
//      const double phononStrength = 3.;
//
//      density = 2330.;
//      workfun = 4.85;
//      bandgap = 1.1;
//      EFermi = -bandgap;
//      potU = -workfun - EFermi;
//      const ElementT* SiComponent[] = { &Element::Si };
//      const double SiComposition[] = { 1. };
//      SEmaterialT Si(SiComponent, 1, SiComposition, 1, density, "Silicon");
//      const double SiWorkfunction = ToSI::eV(workfun);
//      Si.setWorkfunction(SiWorkfunction);
//      Si.setEnergyCBbottom(ToSI::eV(potU));
//      Si.setBandgap(ToSI::eV(bandgap));
//      const double SiCoreEnergy[] = { ToSI::eV(99.2), ToSI::eV(99.8), ToSI::eV(149.7), ToSI::eV(1839.) };
//      Si.setCoreEnergy(SiCoreEnergy, 4);
//
//      tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\";
//      const StringT IIMFPFullPennInterpSiSI = tablePath + "IIMFPFullPennInterpSiSI.csv";
//      const StringT interpNUSimReducedDeltaEFullPennSiSI = tablePath + "interpNUSimReducedDeltaEFullPennSiSI.csv";
//      const StringT interpNUThetaFullPennSiBGSI = tablePath + "interpNUThetaFullPennSiBGSI.csv";
//      const StringT interpSimESE0NUSiBGSI = tablePath + "interpSimESE0NUSiBGSI.csv";
//      const char* SiTables[] = {
//         IIMFPFullPennInterpSiSI.c_str(),
//         interpNUSimReducedDeltaEFullPennSiSI.c_str(),
//         interpNUThetaFullPennSiBGSI.c_str(),
//         interpSimESE0NUSiBGSI.c_str()
//      };
//
//      SelectableElasticSMT SiNISTMott(Si, NISTMottRS::Factory);
//      TabulatedInelasticSMT SiDS(Si, 3, SiTables, ToSI::eV(13.54));
//      printf("%d\n", SiDS.gettableIIMFP()->gettable1d().size());
//      for (auto &i : SiDS.gettableIIMFP()->gettable1d()) {
//         printf("%.10e ", i);
//      }
//      printf("\n");
//
//      printf("%d, %d\n", SiDS.gettableReducedDeltaE()->gettable2d().size(), SiDS.gettableReducedDeltaE()->gettable2d()[0].size());
//      for (auto &i : SiDS.gettableReducedDeltaE()->gettable2d()[0]) {
//         printf("%.10e ", i);
//      }
//      printf("\n");
//
//      printf("%d, %d, %d\n", SiDS.gettableTheta()->gettable3d().size(), SiDS.gettableTheta()->gettable3d()[50].size(), SiDS.gettableTheta()->gettable3d()[50][50].size());
//      for (auto &i : SiDS.gettableTheta()->gettable3d()[50][50]) {
//         printf("%.10e ", i);
//      }
//      printf("\n");
//
//      printf("%d, %d\n", SiDS.gettableSEE0()->gettable2d().size(), SiDS.gettableSEE0()->gettable2d()[0].size());
//      for (auto &i : SiDS.gettableSEE0()->gettable2d()[0]) {
//         printf("%.10e ", i);
//      }
//      printf("\n");
//
//      GanachaudMokraniPhononInelasticSMT Siphonon(phononStrength, ToSI::eV(phononE), 300., 11.7, 1.);
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
//      SiMSMDeep.setMinEforTracking(ToSI::eV(50.));
//
//      const double meterspernm = 1.e-9;
//      const double pitch = pitchnm*meterspernm;
//      const double h = hnm*meterspernm;
//      const double w = wnm*meterspernm;
//      const double linelength = linelengthnm*meterspernm;
//
//      const double radperdeg = Math2::PI / 180.;
//      const double thetar = thetardeg*radperdeg;
//      const double thetal = thetaldeg*radperdeg;
//      const double radr = radrnm*meterspernm;
//      const double radl = radlnm*meterspernm;
//      const double layer1thickness = layer1thicknessnm*meterspernm;
//      const double layer2thickness = layer2thicknessnm*meterspernm;
//      const double beamsize = beamsizenm*meterspernm;
//      const double deep = deepnm*meterspernm;
//
//      NullMaterialScatterModelT NULL_MSM;
//      const double center[] = {
//         0.0,
//         0.0,
//         0.0
//      };
//
//      const double beamEeV = beamEeVvals[0];
//      const double beamE = ToSI::eV(beamEeV);
//      SphereT sphere(center, MonteCarloSS::ChamberRadius);
//      GaussianBeamT eg(beamsize, beamE, center);
//
//      RegionT chamber(nullptr, &NULL_MSM, &sphere);
//      chamber.updateMaterial(*chamber.getScatterModel(), vacuumMSM);
//
//      const double normalvector[] = { 0., 0., -1. };
//      const double layer1Pos[] = { 0., 0., 0. };
//      NormalMultiPlaneShapeT layer1;
//      PlaneT pl1(normalvector, layer1Pos);
//      layer1.addPlane(&pl1);
//      RegionT layer1Region(&chamber, &ARCMSM, (NormalShapeT*)&layer1);
//
//      const double layer2Pos[] = { 0., 0., layer1thickness };
//      NormalMultiPlaneShapeT layer2;
//      PlaneT pl2(normalvector, layer2Pos);
//      layer2.addPlane(&pl2);
//      RegionT layer2Region(&layer1Region, &glCMSM, (NormalShapeT*)&layer2);
//
//      const double layer3Pos[] = { 0., 0., layer1thickness + layer2thickness };
//      NormalMultiPlaneShapeT layer3;
//      PlaneT pl3(normalvector, layer2Pos);
//      layer3.addPlane(&pl3);
//      RegionT layer3Region(&layer2Region, &SiMSM, (NormalShapeT*)&layer3);
//
//      const double layer4Pos[] = { 0., 0., layer1thickness + layer2thickness + deep };
//      NormalMultiPlaneShapeT layer4;
//      PlaneT pl4(normalvector, layer4Pos);
//      RegionT layer4Region(&layer3Region, &SiMSM, (NormalShapeT*)&layer4);
//
//      RegionT deepRegion(&layer3Region, &SiMSMDeep, (NormalShapeT*)&layer4);
//
//      const double leftmostLineCenterx = -pitch*(nlines / 2);
//      //for (int i = 0; i < nlines; ++i) {
//      //   double xcenter = leftmostLineCenterx + i*pitch;
//      //   NormalIntersectionShapeT* line = (NormalIntersectionShapeT*)NShapes::createLine(-h, w, linelength, thetal, thetar, radl, radr);
//      //   const double newLinePos[] = { xcenter, 0., 0. };
//      //   line->translate(newLinePos);
//      //   RegionT lineRegion(&chamber, &PMMAMSM, line);
//      //}
//
//      const double xcenter = leftmostLineCenterx + 0 * pitch;
//      NormalIntersectionShapeT* line = (NormalIntersectionShapeT*)NShapes::createLine(-h, w, linelength, thetal, thetar, radl, radr);
//      RegionT lineRegion(&chamber, &PMMAMSM, line);
//
//      //VectorXd yvals = { 0. };
//      VectorXd yvals(128);
//      //for (int i = -100; i < 100; ++i) {
//      //   yvals.push_back(i);
//      //}
//      for (int i = -64; i < 64; i += 1) {
//         yvals.push_back(i);
//      }
//
//      const double xbottom = wnm / 2.;
//      const double xtop = wnm / 2. - hnm * ::tan(thetar);
//      const double xstart = xbottom - 100.5;
//      const double xstop = xbottom + 100.5;
//      const double xfinestart = xtop - 20.5;
//      double xfinestop;
//      if (thetar < 0.) xfinestop = xtop + 20.5;
//      else xfinestop = wnm / 2. + 20.5;
//
//      VectorXd xvals(128);
//      double deltax = 5.;
//      double x = xstart;
//      while (x < xfinestart) {
//         xvals.push_back(x);
//         x += deltax;
//      }
//      x = xfinestart;
//      deltax = 1;
//      while (x < xfinestop) {
//         xvals.push_back(x);
//         x += deltax;
//      }
//      x = xfinestop;
//      deltax = 5.;
//      while (x < xstop) {
//         xvals.push_back(x);
//         x += deltax;
//      }
//      xvals.push_back(xstop);
//
//      const double binSizeEV = 10.;
//
//      MonteCarloSST monte(&eg, &chamber, eg.createElectron());
//
//      output += "\n# Trajectories at each landing position: " + std::to_string(nTrajectories);
//      output += "\n# Pitch of lines (nm): " + std::to_string(pitchnm);
//      output += "\n# lines: " + std::to_string(nlines);
//      output += "\nLine height (nm): " + std::to_string(hnm);
//      output += "\nLine bottom width (nm): " + std::to_string(wnm);
//      output += "\nLine length (nm): " + std::to_string(linelengthnm);
//      output += "\nLeft and right sidewall angles (deg): " + std::to_string(thetaldeg) + " " + std::to_string(thetardeg);
//      output += "\nLeft and right top corner radii (nm): " + std::to_string(radlnm) + " " + std::to_string(radrnm);
//      output += "\nThicknesses of 1st and second layers (nm): " + std::to_string(layer1thicknessnm) + " " + std::to_string(layer2thicknessnm);
//      output += "\nBeam landing energies (eV): ";
//
//      for (int i = 0; i < beamEeVvalsLen; i++) {
//         output += std::to_string(beamEeVvals[i]);
//      }
//      output += "\nBeam size (standard deviation, in nm): " + std::to_string(beamsizenm);
//
//      output += "\n";
//      output += "\nbeamE (eV)	x(nm)	y (nm)	BSE yield	SE yield";
//
//      //amp::vector<ActionListenerT*> v0;
//      //for (int i = 0; i < 100000; ++i) {
//      //   BackscatterStatsT back(monte, (int)(beamEeV / binSizeEV));
//      //   v0.push_back(&back);
//      //   auto itr = amp::find(v0.begin(), v0.end(), (ActionListenerT*)&back);
//      //if (itr != mEventListeners.end()) {
//      //   v0.erase(itr);
//      //}
//      //}
//
//      auto start = std::chrono::system_clock::now();
//      for (auto ynm : yvals) {
//         //double ynm = yvals[0];
//         double y = ynm*meterspernm;
//         for (auto xnm : xvals) {
//            x = xnm*meterspernm;
//            double egCenter[] = { x, y, -h - 20.*meterspernm };
//            eg.setCenter(egCenter);
//
//            static const int nbins = (int)(beamEeV / binSizeEV);
//            BackscatterStatsT back(monte, nbins);
//            monte.addActionListener(back);
//
//            //try {
//               monte.runMultipleTrajectories(nTrajectories);
//
//               const HistogramT& hist = back.backscatterEnergyHistogram();
//
//               double energyperbineV = beamEeV / hist.binCount();
//               double maxSEbin = 50. / energyperbineV;
//               int totalSE = 0;
//               for (int j = 0; j < (int)maxSEbin; ++j) {
//                  totalSE = totalSE + hist.counts(j);
//               }
//
//               double SEf = (float)totalSE / nTrajectories;
//               double bsf = back.backscatterFraction() - SEf;
//               output += "\n" + std::to_string(beamEeV) + " " + std::to_string(xnm) + " " + std::to_string(ynm) + " " + std::to_string(bsf) + " " + std::to_string(SEf);
//
//               std::string tmp("\n" + std::to_string(beamEeV) + " " + std::to_string(xnm) + " " + std::to_string(ynm) + " " + std::to_string(bsf) + " " + std::to_string(SEf));
//               printf("%s", tmp.c_str());
//
//               monte.removeActionListener(back);
//            //}
//            //catch (std::exception&) {
//            //   printf("\nwtfweewew");
//            //}
//         }
//      }
//      auto end = std::chrono::system_clock::now();
//      std::chrono::duration<double> elapsed_seconds = end - start;
//      std::time_t end_time = std::chrono::system_clock::to_time_t(end);
//      std::cout << std::endl << "finished computation at " << std::ctime(&end_time) << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
//      output += "\n" + std::to_string(elapsed_seconds.count());
//
//      std::ofstream myfile;
//      myfile.open("output.txt");
//      myfile << output.c_str();
//      myfile.close();
//   }
//
//   void run1()
//   {
//      const int seed = rand();
//      srand(seed);
//
//      printf("\n Random number seed: ");
//      for (int i = 0; i < 10; ++i) {
//         printf("\n %.10e", Random::random());
//      }
//
//      const int nTrajectories = 100;
//
//      const double pitchnm = 180;
//      const int nlines = 3;
//      const double hnm = 120;
//      const double wnm = 80;
//      const double linelengthnm = 1000;
//      const double thetardeg = 3;
//      const double thetaldeg = 3;
//      const double radrnm = 20;
//      const double radlnm = 20;
//      const double layer1thicknessnm = 80;
//      const double layer2thicknessnm = 200;
//
//      const double beamEeVvals[] = { 500. };
//      const int beamEeVvalsLen = 1;
//      const double beamsizenm = 0.5;
//      const double deepnm = 15;
//
//      const bool trajImg = true;
//      const int trajImgMaxTraj = 50;
//      const double trajImgSize = 200e-9;
//
//      const bool VRML = false;
//      const int VRMLImgMaxTraj = 0;
//
//      SEmaterialT* vacuum = new SEmaterialT();
//      vacuum->setName("SE vacuum");
//      ExpQMBarrierSMT* vacuumBarrier = new ExpQMBarrierSMT(vacuum);
//      ZeroCSDT* sZeroCSD = new ZeroCSDT();
//
//      MONSEL_MaterialScatterModelT* vacuumMSM = new MONSEL_MaterialScatterModelT(vacuum, vacuumBarrier, sZeroCSD);
//
//      const double breakE = ToSI::eV(45.);
//      double density = 1190.;
//      double workfun = 5.5;
//      double bandgap = 5.;
//      double EFermi = -bandgap;
//      double potU = -workfun - EFermi;
//
//      //const ElementT* componentsCOH[] = { &Element::C, &Element::O, &Element::H };
//      //const double compositionCOH[] = { 5. / 15., 2. / 15., 8. / 15. };
//      //SEmaterialT* PMMA = new SEmaterialT(componentsCOH, 3, compositionCOH, 3, density, "PMMA");
//      //PMMA->setWorkfunction(ToSI::eV(workfun));
//      //PMMA->setBandgap(ToSI::eV(bandgap));
//      //PMMA->setEnergyCBbottom(ToSI::eV(potU));
//
//      const ElementT& C = Element::C;
//      const ElementT& Ox = Element::O;
//      const ElementT& H = Element::H;
//      CompositionT PMMAcomp;
//      const ElementT* componentsCOH[] = { &C, &Ox, &H };
//      const double compositionCOH[] = { 5, 2, 8 };
//      PMMAcomp.defineByMoleFraction(componentsCOH, 3, compositionCOH, 3);
//      SEmaterialT* PMMA = new SEmaterialT(PMMAcomp, density);
//      PMMA->setName("PMMA");
//      PMMA->setWorkfunction(ToSI::eV(workfun));
//      PMMA->setBandgap(ToSI::eV(bandgap));
//      PMMA->setEnergyCBbottom(ToSI::eV(potU));
//
//      SelectableElasticSMT* PMMANISTMott = new SelectableElasticSMT(*PMMA, NISTMottRS::Factory);
//
//      JoyLuoNieminenCSDT* PMMACSD = new JoyLuoNieminenCSDT(*PMMA, breakE);
//      FittedInelSMT* PMMAfittedInel = new FittedInelSMT(*PMMA, ToSI::eV(65.4), *PMMACSD);
//      GanachaudMokraniPolaronTrapSMT* PMMApolaron = new GanachaudMokraniPolaronTrapSMT(2.e7, 1. / ToSI::eV(4.));
//
//      ExpQMBarrierSMT* pmmaeqmbsm = new ExpQMBarrierSMT(PMMA);
//
//      MONSEL_MaterialScatterModelT* PMMAMSM = new MONSEL_MaterialScatterModelT(PMMA, pmmaeqmbsm, sZeroCSD);
//      PMMAMSM->addScatterMechanism(PMMANISTMott);
//      PMMAMSM->addScatterMechanism(PMMAfittedInel);
//      PMMAMSM->addScatterMechanism(PMMApolaron);
//
//      PMMAMSM->setCSD(PMMACSD);
//
//      MONSEL_MaterialScatterModelT* PMMAMSMDeep = new MONSEL_MaterialScatterModelT(PMMA, pmmaeqmbsm, sZeroCSD);
//      PMMAMSMDeep->addScatterMechanism(PMMANISTMott);
//      PMMAMSMDeep->addScatterMechanism(PMMAfittedInel);
//      PMMAMSMDeep->addScatterMechanism(PMMApolaron);
//
//      PMMAMSMDeep->setCSD(PMMACSD);
//      PMMAMSMDeep->setMinEforTracking(ToSI::eV(50.));
//
//      MONSEL_MaterialScatterModelT* ARCMSM = PMMAMSM;
//
//      density = 1800.;
//      workfun = 5.0;
//      bandgap = 0.;
//      EFermi = 20.4;
//      potU = -workfun - EFermi;
//      const ElementT* glCComponents[] = { &Element::C };
//      const double glCComposition[] = { 1. };
//      SEmaterialT* glC = new SEmaterialT(glCComponents, 1, glCComposition, 1, density, "glassy Carbon");
//      glC->setWorkfunction(ToSI::eV(workfun));
//      glC->setEnergyCBbottom(ToSI::eV(potU));
//      glC->setBandgap(ToSI::eV(bandgap));
//      const double glCCoreEnergy[] = { ToSI::eV(284.2) };
//      glC->setCoreEnergy(glCCoreEnergy, 1);
//
//      SelectableElasticSMT* glCNISTMott = new SelectableElasticSMT(*glC, NISTMottRS::Factory);
//
//      StringT tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\";
//      const StringT IIMFPPennInterpglassy = tablePath + "IIMFPPennInterpglassyCSI.csv";
//      const StringT SimReducedDeltaEglassy = tablePath + "interpNUSimReducedDeltaEglassyCSI.csv";
//      const StringT simTableThetaNUglassy = tablePath + "interpsimTableThetaNUglassyCSI.csv";
//      const StringT SimESE0NUglassy = tablePath + "interpSimESE0NUglassyCSI.csv";
//      const char* glCTables[] = {
//         IIMFPPennInterpglassy.c_str(),
//         SimReducedDeltaEglassy.c_str(),
//         simTableThetaNUglassy.c_str(),
//         SimESE0NUglassy.c_str()
//      };
//
//      TabulatedInelasticSMT* glCDS = new TabulatedInelasticSMT(*glC, 3, glCTables);
//
//      ExpQMBarrierSMT* glceqmbsm = new ExpQMBarrierSMT(glC);
//
//      MONSEL_MaterialScatterModelT* glCMSM = new MONSEL_MaterialScatterModelT(glC, glceqmbsm, sZeroCSD);
//      glCMSM->addScatterMechanism(glCNISTMott);
//      glCMSM->addScatterMechanism(glCDS);
//
//      MONSEL_MaterialScatterModelT* glCMSMDeep = new MONSEL_MaterialScatterModelT(glC, glceqmbsm, sZeroCSD);
//      glCMSMDeep->addScatterMechanism(glCNISTMott);
//      glCMSMDeep->addScatterMechanism(glCDS);
//
//      glCMSMDeep->setMinEforTracking(ToSI::eV(50.));
//
//      const double phononE = 0.063;
//      const double phononStrength = 3.;
//
//      density = 2330.;
//      workfun = 4.85;
//      bandgap = 1.1;
//      EFermi = -bandgap;
//      potU = -workfun - EFermi;
//      const ElementT* SiComponent[] = { &Element::Si };
//      const double SiComposition[] = { 1. };
//      SEmaterialT* Si = new SEmaterialT(SiComponent, 1, SiComposition, 1, density, "Silicon");
//      const double SiWorkfunction = ToSI::eV(workfun);
//      Si->setWorkfunction(SiWorkfunction);
//      Si->setEnergyCBbottom(ToSI::eV(potU));
//      Si->setBandgap(ToSI::eV(bandgap));
//      const double SiCoreEnergy[] = { ToSI::eV(99.2), ToSI::eV(99.8), ToSI::eV(149.7), ToSI::eV(1839.) };
//      Si->setCoreEnergy(SiCoreEnergy, 4);
//
//      SelectableElasticSMT* SiNISTMott = new SelectableElasticSMT(*Si, NISTMottRS::Factory);
//
//      tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\";
//      const StringT IIMFPFullPennInterpSiSI = tablePath + "IIMFPFullPennInterpSiSI.csv";
//      const StringT interpNUSimReducedDeltaEFullPennSiSI = tablePath + "interpNUSimReducedDeltaEFullPennSiSI.csv";
//      const StringT interpNUThetaFullPennSiBGSI = tablePath + "interpNUThetaFullPennSiBGSI.csv";
//      const StringT interpSimESE0NUSiBGSI = tablePath + "interpSimESE0NUSiBGSI.csv";
//      const char* SiTables[] = {
//         IIMFPFullPennInterpSiSI.c_str(),
//         interpNUSimReducedDeltaEFullPennSiSI.c_str(),
//         interpNUThetaFullPennSiBGSI.c_str(),
//         interpSimESE0NUSiBGSI.c_str()
//      };
//
//      TabulatedInelasticSMT* SiDS = new TabulatedInelasticSMT(*Si, 3, SiTables, ToSI::eV(13.54));
//
//      GanachaudMokraniPhononInelasticSMT* Siphonon = new GanachaudMokraniPhononInelasticSMT(phononStrength, ToSI::eV(phononE), 300., 11.7, 1.);
//
//      ExpQMBarrierSMT* sieqmbsm = new ExpQMBarrierSMT(Si);
//
//      MONSEL_MaterialScatterModelT* SiMSM = new MONSEL_MaterialScatterModelT(Si, sieqmbsm, sZeroCSD);
//      SiMSM->addScatterMechanism(SiNISTMott);
//      SiMSM->addScatterMechanism(SiDS);
//      SiMSM->addScatterMechanism(Siphonon);
//
//      MONSEL_MaterialScatterModelT* SiMSMDeep = new MONSEL_MaterialScatterModelT(Si, sieqmbsm, sZeroCSD);
//      SiMSMDeep->addScatterMechanism(SiNISTMott);
//      SiMSMDeep->addScatterMechanism(SiDS);
//      SiMSMDeep->addScatterMechanism(Siphonon);
//
//      SiMSMDeep->setMinEforTracking(ToSI::eV(50.));
//
//      const double meterspernm = 1.e-9;
//      const double pitch = pitchnm*meterspernm;
//      const double h = hnm*meterspernm;
//      const double w = wnm*meterspernm;
//      const double linelength = linelengthnm*meterspernm;
//
//      const double radperdeg = Math2::PI / 180.;
//      const double thetar = thetardeg*radperdeg;
//      const double thetal = thetaldeg*radperdeg;
//      const double radr = radrnm*meterspernm;
//      const double radl = radlnm*meterspernm;
//      const double layer1thickness = layer1thicknessnm*meterspernm;
//      const double layer2thickness = layer2thicknessnm*meterspernm;
//      const double beamsize = beamsizenm*meterspernm;
//      const double deep = deepnm*meterspernm;
//
//      const double center[] = {
//         0.0,
//         0.0,
//         0.0
//      };
//
//      const double beamEeV = beamEeVvals[0];
//      const double beamE = ToSI::eV(beamEeV);
//      SphereT* sphere = new SphereT(center, MonteCarloSS::ChamberRadius);
//      GaussianBeamT* eg = new GaussianBeamT(beamsize, beamE, center);
//
//      NullMaterialScatterModelT* NULL_MSM = new NullMaterialScatterModelT();
//      RegionT* chamber = new RegionT(nullptr, NULL_MSM, sphere);
//      chamber->updateMaterial(*(chamber->getScatterModel()), *vacuumMSM);
//
//      const double normalvector[] = { 0., 0., -1. };
//      const double layer1Pos[] = { 0., 0., 0. };
//      NormalMultiPlaneShapeT* layer1 = new NormalMultiPlaneShapeT();
//      PlaneT* pl1 = new PlaneT(normalvector, layer1Pos);
//      layer1->addPlane(pl1);
//      RegionT* layer1Region = new RegionT(chamber, ARCMSM, (NormalShapeT*)layer1);
//
//      const double layer2Pos[] = { 0., 0., layer1thickness };
//      NormalMultiPlaneShapeT* layer2 = new NormalMultiPlaneShapeT();
//      PlaneT* pl2 = new PlaneT(normalvector, layer2Pos);
//      layer2->addPlane(pl2);
//      RegionT* layer2Region = new RegionT(layer1Region, glCMSM, (NormalShapeT*)layer2);
//
//      const double layer3Pos[] = { 0., 0., layer1thickness + layer2thickness };
//      NormalMultiPlaneShapeT* layer3 = new NormalMultiPlaneShapeT();
//      PlaneT* pl3 = new PlaneT(normalvector, layer2Pos);
//      layer3->addPlane(pl3);
//      RegionT* layer3Region = new RegionT(layer2Region, SiMSM, (NormalShapeT*)layer3);
//
//      const double layer4Pos[] = { 0., 0., layer1thickness + layer2thickness + deep };
//      NormalMultiPlaneShapeT* layer4 = new NormalMultiPlaneShapeT();
//      PlaneT* pl4 = new PlaneT(normalvector, layer4Pos);
//      RegionT* layer4Region = new RegionT(layer3Region, SiMSM, (NormalShapeT*)layer4);
//
//      RegionT* deepRegion = new RegionT(layer3Region, SiMSMDeep, (NormalShapeT*)layer4);
//
//      const double leftmostLineCenterx = -pitch*(nlines / 2);
//      //for (int i = 0; i < nlines; ++i) {
//      //   double xcenter = leftmostLineCenterx + i*pitch;
//      //   NormalIntersectionShapeT* line = (NormalIntersectionShapeT*)NShapes::createLine(-h, w, linelength, thetal, thetar, radl, radr);
//      //   const double newLinePos[] = { xcenter, 0., 0. };
//      //   line->translate(newLinePos);
//      //   RegionT lineRegion(&chamber, &PMMAMSM, line);
//      //}
//
//      const double xcenter = leftmostLineCenterx + 0 * pitch;
//      NormalIntersectionShapeT* line = (NormalIntersectionShapeT*)NShapes::createLine(-h, w, linelength, thetal, thetar, radl, radr);
//      RegionT* lineRegion = new RegionT(chamber, PMMAMSM, line);
//
//      //VectorXd yvals = { 0. };
//      VectorXd yvals(128);
//      //for (int i = -100; i < 100; ++i) {
//      //   yvals.push_back(i);
//      //}
//      for (int i = -64; i < 64; i += 1) {
//         yvals.push_back(i);
//      }
//
//      const double xbottom = wnm / 2.;
//      const double xtop = wnm / 2. - hnm * ::tan(thetar);
//      const double xstart = xbottom - 100.5;
//      const double xstop = xbottom + 100.5;
//      const double xfinestart = xtop - 20.5;
//      double xfinestop;
//      if (thetar < 0.) xfinestop = xtop + 20.5;
//      else xfinestop = wnm / 2. + 20.5;
//
//      VectorXd xvals(128);
//      double deltax = 5.;
//      double x = xstart;
//      while (x < xfinestart) {
//         xvals.push_back(x);
//         x += deltax;
//      }
//      x = xfinestart;
//      deltax = 1;
//      while (x < xfinestop) {
//         xvals.push_back(x);
//         x += deltax;
//      }
//      x = xfinestop;
//      deltax = 5.;
//      while (x < xstop) {
//         xvals.push_back(x);
//         x += deltax;
//      }
//      xvals.push_back(xstop);
//
//      const double binSizeEV = 10.;
//
//      MonteCarloSST* monte = new MonteCarloSST(eg, chamber, eg->createElectron());
//
//      printf("\n# Trajectories at each landing position: %d", nTrajectories);
//      printf("\n# Pitch of lines (nm): %.10e", pitchnm);
//      printf("\n# lines: %d", nlines);
//      printf("\nLine height (nm): %.10e", hnm);
//      printf("\nLine bottom width (nm): %.10e", wnm);
//      printf("\nLine length (nm): %.10e", linelengthnm);
//      printf("\nLeft and right sidewall angles (deg): %.10e %.10e", thetaldeg, thetardeg);
//      printf("\nLeft and right top corner radii (nm): %.10e %.10e", radlnm, radrnm);
//      printf("\nThicknesses of 1st and second layers (nm): %.10e %.10e", layer1thicknessnm, layer2thicknessnm);
//      printf("\nBeam landing energies (eV): ");
//
//      for (int i = 0; i < beamEeVvalsLen; i++) {
//         printf("\n%.10e", beamEeVvals[i]);
//      }
//      printf("\nBeam size (standard deviation, in nm): %.10e", beamsizenm);
//
//      printf("\n");
//      printf("\nbeamE (eV)\t x(nm)\t y (nm)\t BSE yield\t SE yield");
//
//      std::string output;
//      auto start = std::chrono::system_clock::now();
//      for (auto ynm : yvals) {
//         //double ynm = yvals[0];
//         const double y = ynm*meterspernm;
//         for (auto xnm : xvals) {
//            x = xnm*meterspernm;
//            double egCenter[] = { x, y, -h - 20. * meterspernm };
//            eg->setCenter(egCenter);
//
//            const int nbins = (int)(beamEeV / binSizeEV);
//            BackscatterStatsT* back = new BackscatterStatsT(*monte, nbins);
//            monte->addActionListener(*back);
//
//            monte->runMultipleTrajectories(nTrajectories);
//
//            const HistogramT& hist = back->backscatterEnergyHistogram();
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
//            printf("\n %.10e %.10e %.10e %.10e %.10e", beamEeV, xnm, ynm, bsf, SEf);
//            output += "\n" + std::to_string(beamEeV) + " " + std::to_string(xnm) + " " + std::to_string(ynm) + " " + std::to_string(bsf) + " " + std::to_string(SEf);
//
//            monte->removeActionListener(*back);
//         }
//      }
//      auto end = std::chrono::system_clock::now();
//      std::chrono::duration<double> elapsed_seconds = end - start;
//      std::time_t end_time = std::chrono::system_clock::to_time_t(end);
//      std::cout << std::endl << "finished computation at " << std::ctime(&end_time) << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
//      output += "\n" + std::to_string(elapsed_seconds.count());
//      printf("\n");
//
//      std::ofstream myfile;
//      myfile.open("output.txt");
//      myfile << output.c_str();
//      myfile.close();
//   }
//
//   void initCuda()
//   {
//      NUTableInterpolation::initFactory << <1, 1 >> >();
//      checkCudaErrors(cudaDeviceSynchronize());
//      checkCudaErrors(cudaGetLastError());
//
//      StringT tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\";
//      const StringT IIMFPPennInterpglassy = tablePath + "IIMFPPennInterpglassyCSI.csv";
//      const StringT SimReducedDeltaEglassy = tablePath + "interpNUSimReducedDeltaEglassyCSI.csv";
//      const StringT simTableThetaNUglassy = tablePath + "interpsimTableThetaNUglassyCSI.csv";
//      const StringT SimESE0NUglassy = tablePath + "interpSimESE0NUglassyCSI.csv";
//      NUTableInterpolation::copyDataToCuda(IIMFPPennInterpglassy.c_str());
//      NUTableInterpolation::copyDataToCuda(SimReducedDeltaEglassy.c_str());
//      NUTableInterpolation::copyDataToCuda(simTableThetaNUglassy.c_str());
//      NUTableInterpolation::copyDataToCuda(SimESE0NUglassy.c_str());
//
//      tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\";
//      const StringT IIMFPFullPennInterpSiSI = tablePath + "IIMFPFullPennInterpSiSI.csv";
//      const StringT interpNUSimReducedDeltaEFullPennSiSI = tablePath + "interpNUSimReducedDeltaEFullPennSiSI.csv";
//      const StringT interpNUThetaFullPennSiBGSI = tablePath + "interpNUThetaFullPennSiBGSI.csv";
//      const StringT interpSimESE0NUSiBGSI = tablePath + "interpSimESE0NUSiBGSI.csv";
//      NUTableInterpolation::copyDataToCuda(IIMFPFullPennInterpSiSI.c_str());
//      NUTableInterpolation::copyDataToCuda(interpNUSimReducedDeltaEFullPennSiSI.c_str());
//      NUTableInterpolation::copyDataToCuda(interpNUThetaFullPennSiBGSI.c_str());
//      NUTableInterpolation::copyDataToCuda(interpSimESE0NUSiBGSI.c_str());
//   }
//
//   __global__ void runCuda()
//   {
//      //String dest = DefaultOutput;
//      //new File(dest).mkdirs();
//      //String filename = DefaultOutput + PathSep + "output.txt";
//
//      //const int seed = rand();
//      //const int seed = 0;
//      //srand(seed);
//      printf("\n Random number seed: ");
//      //output += amp::to_string(seed);
//      for (int i = 0; i < 10; ++i) {
//         printf("\n %.10e", Random::random());
//      }
//
//      const int nTrajectories = 100;
//
//      const double pitchnm = 180;
//      const int nlines = 3;
//      const double hnm = 120;
//      const double wnm = 80;
//      const double linelengthnm = 1000;
//      const double thetardeg = 3;
//      const double thetaldeg = 3;
//      const double radrnm = 20;
//      const double radlnm = 20;
//      const double layer1thicknessnm = 80;
//      const double layer2thicknessnm = 200;
//
//      const double beamEeVvals[] = { 500. };
//      const int beamEeVvalsLen = 1;
//      const double beamsizenm = 0.5;
//      const double deepnm = 15;
//
//      const bool trajImg = true;
//      const int trajImgMaxTraj = 50;
//      const double trajImgSize = 200e-9;
//
//      const bool VRML = false;
//      const int VRMLImgMaxTraj = 0;
//
//      SEmaterialT* vacuum = new SEmaterialT(); printf("0");
//      vacuum->setName("SE vacuum");
//      ExpQMBarrierSMT* vacuumBarrier = new ExpQMBarrierSMT(vacuum); printf("1");
//      ZeroCSDT* sZeroCSD = new ZeroCSDT(); printf("2");
//
//      MONSEL_MaterialScatterModelT* vacuumMSM = new MONSEL_MaterialScatterModelT(vacuum, vacuumBarrier, sZeroCSD); printf("3");
//
//      const double breakE = ToSI::eV(45.);
//      double density = 1190.;
//      double workfun = 5.5;
//      double bandgap = 5.;
//      double EFermi = -bandgap;
//      double potU = -workfun - EFermi;
//
//      const ElementT* componentsCOH[] = { Element::dC, Element::dO, Element::dH };
//      const double compositionCOH[] = { 5. / 15., 2. / 15., 8. / 15. };
//      SEmaterialT* PMMA = new SEmaterialT(componentsCOH, 3, compositionCOH, 3, density, "PMMA"); printf("4");
//      PMMA->setWorkfunction(ToSI::eV(workfun));
//      PMMA->setBandgap(ToSI::eV(bandgap));
//      PMMA->setEnergyCBbottom(ToSI::eV(potU));
//
//      SelectableElasticSMT* PMMANISTMott = new SelectableElasticSMT(*PMMA, *NISTMottRS::d_Factory); printf("5");
//
//      JoyLuoNieminenCSDT* PMMACSD = new JoyLuoNieminenCSDT(*PMMA, breakE); printf("6");
//      FittedInelSMT* PMMAfittedInel = new FittedInelSMT(*PMMA, ToSI::eV(65.4), *PMMACSD); printf("7");
//      GanachaudMokraniPolaronTrapSMT* PMMApolaron = new GanachaudMokraniPolaronTrapSMT(2.e7, 1. / ToSI::eV(4.)); printf("8");
//
//      ExpQMBarrierSMT* pmmaeqmbsm = new ExpQMBarrierSMT(PMMA); printf("9");
//
//      MONSEL_MaterialScatterModelT* PMMAMSM = new MONSEL_MaterialScatterModelT(PMMA, pmmaeqmbsm, sZeroCSD); printf("10");
//      PMMAMSM->addScatterMechanism(PMMANISTMott);
//      PMMAMSM->addScatterMechanism(PMMAfittedInel);
//      PMMAMSM->addScatterMechanism(PMMApolaron);
//
//      PMMAMSM->setCSD(PMMACSD);
//
//      MONSEL_MaterialScatterModelT* PMMAMSMDeep = new MONSEL_MaterialScatterModelT(PMMA, pmmaeqmbsm, sZeroCSD); printf("11");
//      PMMAMSMDeep->addScatterMechanism(PMMANISTMott);
//      PMMAMSMDeep->addScatterMechanism(PMMAfittedInel);
//      PMMAMSMDeep->addScatterMechanism(PMMApolaron);
//
//      PMMAMSMDeep->setCSD(PMMACSD);
//      PMMAMSMDeep->setMinEforTracking(ToSI::eV(50.));
//
//      MONSEL_MaterialScatterModelT* ARCMSM = PMMAMSM; printf("12");
//
//      density = 1800.;
//      workfun = 5.0;
//      bandgap = 0.;
//      EFermi = 20.4;
//      potU = -workfun - EFermi;
//      const ElementT* glCComponents[] = { Element::dC };
//      const double glCComposition[] = { 1. };
//      SEmaterialT* glC = new SEmaterialT(glCComponents, 1, glCComposition, 1, density, "glassy Carbon"); printf("13");
//      glC->setWorkfunction(ToSI::eV(workfun));
//      glC->setEnergyCBbottom(ToSI::eV(potU));
//      glC->setBandgap(ToSI::eV(bandgap));
//      const double glCCoreEnergy[] = { ToSI::eV(284.2) };
//      glC->setCoreEnergy(glCCoreEnergy, 1);
//
//      SelectableElasticSMT* glCNISTMott = new SelectableElasticSMT(*glC, *NISTMottRS::d_Factory); printf("14");
//
//      StringT tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\glassyCTables\\";
//      const StringT IIMFPPennInterpglassy = tablePath + "IIMFPPennInterpglassyCSI.csv";
//      const StringT SimReducedDeltaEglassy = tablePath + "interpNUSimReducedDeltaEglassyCSI.csv";
//      const StringT simTableThetaNUglassy = tablePath + "interpsimTableThetaNUglassyCSI.csv";
//      const StringT SimESE0NUglassy = tablePath + "interpSimESE0NUglassyCSI.csv";
//      const char* glCTables[] = {
//         IIMFPPennInterpglassy.c_str(),
//         SimReducedDeltaEglassy.c_str(),
//         simTableThetaNUglassy.c_str(),
//         SimESE0NUglassy.c_str()
//      };
//
//      TabulatedInelasticSMT* glCDS = new TabulatedInelasticSMT(*glC, 3, glCTables); printf("15");
//
//      ExpQMBarrierSMT* glceqmbsm = new ExpQMBarrierSMT(glC); printf("16");
//
//      MONSEL_MaterialScatterModelT* glCMSM = new MONSEL_MaterialScatterModelT(glC, glceqmbsm, sZeroCSD); printf("17");
//      glCMSM->addScatterMechanism(glCNISTMott);
//      glCMSM->addScatterMechanism(glCDS);
//
//      MONSEL_MaterialScatterModelT* glCMSMDeep = new MONSEL_MaterialScatterModelT(glC, glceqmbsm, sZeroCSD); printf("18");
//      glCMSMDeep->addScatterMechanism(glCNISTMott);
//      glCMSMDeep->addScatterMechanism(glCDS);
//
//      glCMSMDeep->setMinEforTracking(ToSI::eV(50.));
//
//      const double phononE = 0.063;
//      const double phononStrength = 3.;
//
//      density = 2330.;
//      workfun = 4.85;
//      bandgap = 1.1;
//      EFermi = -bandgap;
//      potU = -workfun - EFermi;
//      const ElementT* SiComponent[] = { Element::dSi };
//      const double SiComposition[] = { 1. };
//      SEmaterialT* Si = new SEmaterialT(SiComponent, 1, SiComposition, 1, density, "Silicon"); printf("19");
//      const double SiWorkfunction = ToSI::eV(workfun);
//      Si->setWorkfunction(SiWorkfunction);
//      Si->setEnergyCBbottom(ToSI::eV(potU));
//      Si->setBandgap(ToSI::eV(bandgap));
//      const double SiCoreEnergy[] = { ToSI::eV(99.2), ToSI::eV(99.8), ToSI::eV(149.7), ToSI::eV(1839.) };
//      Si->setCoreEnergy(SiCoreEnergy, 4);
//
//      SelectableElasticSMT* SiNISTMott = new SelectableElasticSMT(*Si, *NISTMottRS::d_Factory); printf("20");
//
//      tablePath = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\";
//      const StringT IIMFPFullPennInterpSiSI = tablePath + "IIMFPFullPennInterpSiSI.csv";
//      const StringT interpNUSimReducedDeltaEFullPennSiSI = tablePath + "interpNUSimReducedDeltaEFullPennSiSI.csv";
//      const StringT interpNUThetaFullPennSiBGSI = tablePath + "interpNUThetaFullPennSiBGSI.csv";
//      const StringT interpSimESE0NUSiBGSI = tablePath + "interpSimESE0NUSiBGSI.csv";
//      const char* SiTables[] = {
//         IIMFPFullPennInterpSiSI.c_str(),
//         interpNUSimReducedDeltaEFullPennSiSI.c_str(),
//         interpNUThetaFullPennSiBGSI.c_str(),
//         interpSimESE0NUSiBGSI.c_str()
//      };
//
//      TabulatedInelasticSMT* SiDS = new TabulatedInelasticSMT(*Si, 3, SiTables, ToSI::eV(13.54)); printf("21");
//
//      GanachaudMokraniPhononInelasticSMT* Siphonon = new GanachaudMokraniPhononInelasticSMT(phononStrength, ToSI::eV(phononE), 300., 11.7, 1.); printf("22");
//
//      ExpQMBarrierSMT* sieqmbsm = new ExpQMBarrierSMT(Si); printf("23");
//
//      MONSEL_MaterialScatterModelT* SiMSM = new MONSEL_MaterialScatterModelT(Si, sieqmbsm, sZeroCSD); printf("24");
//      SiMSM->addScatterMechanism(SiNISTMott);
//      SiMSM->addScatterMechanism(SiDS);
//      SiMSM->addScatterMechanism(Siphonon);
//
//      MONSEL_MaterialScatterModelT* SiMSMDeep = new MONSEL_MaterialScatterModelT(Si, sieqmbsm, sZeroCSD); printf("25");
//      SiMSMDeep->addScatterMechanism(SiNISTMott);
//      SiMSMDeep->addScatterMechanism(SiDS);
//      SiMSMDeep->addScatterMechanism(Siphonon);
//
//      SiMSMDeep->setMinEforTracking(ToSI::eV(50.));
//
//      const double meterspernm = 1.e-9;
//      const double pitch = pitchnm*meterspernm;
//      const double h = hnm*meterspernm;
//      const double w = wnm*meterspernm;
//      const double linelength = linelengthnm*meterspernm;
//
//      const double radperdeg = Math2::PI / 180.;
//      const double thetar = thetardeg*radperdeg;
//      const double thetal = thetaldeg*radperdeg;
//      const double radr = radrnm*meterspernm;
//      const double radl = radlnm*meterspernm;
//      const double layer1thickness = layer1thicknessnm*meterspernm;
//      const double layer2thickness = layer2thicknessnm*meterspernm;
//      const double beamsize = beamsizenm*meterspernm;
//      const double deep = deepnm*meterspernm;
//
//      const double center[] = {
//         0.0,
//         0.0,
//         0.0
//      };
//
//      const double beamEeV = beamEeVvals[0];
//      const double beamE = ToSI::eV(beamEeV);
//      SphereT* sphere = new SphereT(center, MonteCarloSS::ChamberRadius); printf("26");
//      GaussianBeamT* eg = new GaussianBeamT(beamsize, beamE, center); printf("27");
//
//      NullMaterialScatterModelT* NULL_MSM = new NullMaterialScatterModelT(); printf("28");
//      RegionT* chamber = new RegionT(nullptr, NULL_MSM, sphere); printf("29");
//      chamber->updateMaterial(*(chamber->getScatterModel()), *vacuumMSM);
//
//      const double normalvector[] = { 0., 0., -1. };
//      const double layer1Pos[] = { 0., 0., 0. };
//      NormalMultiPlaneShapeT* layer1 = new NormalMultiPlaneShapeT(); printf("30");
//      PlaneT* pl1 = new PlaneT(normalvector, layer1Pos); printf("31");
//      layer1->addPlane(pl1);
//      RegionT* layer1Region = new RegionT(chamber, ARCMSM, (NormalShapeT*)layer1); printf("32");
//
//      const double layer2Pos[] = { 0., 0., layer1thickness };
//      NormalMultiPlaneShapeT* layer2 = new NormalMultiPlaneShapeT(); printf("33");
//      PlaneT* pl2 = new PlaneT(normalvector, layer2Pos); printf("34");
//      layer2->addPlane(pl2);
//      RegionT* layer2Region = new RegionT(layer1Region, glCMSM, (NormalShapeT*)layer2); printf("35");
//
//      const double layer3Pos[] = { 0., 0., layer1thickness + layer2thickness };
//      NormalMultiPlaneShapeT* layer3 = new NormalMultiPlaneShapeT(); printf("36");
//      PlaneT* pl3 = new PlaneT(normalvector, layer2Pos); printf("37");
//      layer3->addPlane(pl3);
//      RegionT* layer3Region = new RegionT(layer2Region, SiMSM, (NormalShapeT*)layer3); printf("38");
//
//      const double layer4Pos[] = { 0., 0., layer1thickness + layer2thickness + deep };
//      NormalMultiPlaneShapeT* layer4 = new NormalMultiPlaneShapeT(); printf("39");
//      PlaneT* pl4 = new PlaneT(normalvector, layer4Pos); printf("40");
//      RegionT* layer4Region = new RegionT(layer3Region, SiMSM, (NormalShapeT*)layer4); printf("41");
//
//      RegionT* deepRegion = new RegionT(layer3Region, SiMSMDeep, (NormalShapeT*)layer4); printf("42");
//
//      const double leftmostLineCenterx = -pitch*(nlines / 2);
//      //for (int i = 0; i < nlines; ++i) {
//      //   double xcenter = leftmostLineCenterx + i*pitch;
//      //   NormalIntersectionShapeT* line = (NormalIntersectionShapeT*)NShapes::createLine(-h, w, linelength, thetal, thetar, radl, radr);
//      //   const double newLinePos[] = { xcenter, 0., 0. };
//      //   line->translate(newLinePos);
//      //   RegionT lineRegion(&chamber, &PMMAMSM, line);
//      //}
//
//      const double xcenter = leftmostLineCenterx + 0 * pitch;
//      NormalIntersectionShapeT* line = (NormalIntersectionShapeT*)NShapes::createLine(-h, w, linelength, thetal, thetar, radl, radr); printf("43");
//      RegionT* lineRegion = new RegionT(chamber, PMMAMSM, line); printf("44");
//
//      //VectorXd yvals = { 0. };
//      VectorXd yvals(128); printf("45");
//      //for (int i = -100; i < 100; ++i) {
//      //   yvals.push_back(i);
//      //}
//      for (int i = -64; i < 64; i += 1) {
//         yvals.push_back(i);
//      }
//
//      const double xbottom = wnm / 2.;
//      const double xtop = wnm / 2. - hnm * ::tan(thetar);
//      const double xstart = xbottom - 100.5;
//      const double xstop = xbottom + 100.5;
//      const double xfinestart = xtop - 20.5;
//      double xfinestop;
//      if (thetar < 0.) xfinestop = xtop + 20.5;
//      else xfinestop = wnm / 2. + 20.5;
//
//      VectorXd xvals(128); printf("46");
//      double deltax = 5.;
//      double x = xstart;
//      while (x < xfinestart) {
//         xvals.push_back(x);
//         x += deltax;
//      }
//      x = xfinestart;
//      deltax = 1;
//      while (x < xfinestop) {
//         xvals.push_back(x);
//         x += deltax;
//      }
//      x = xfinestop;
//      deltax = 5.;
//      while (x < xstop) {
//         xvals.push_back(x);
//         x += deltax;
//      }
//      xvals.push_back(xstop);
//
//      const double binSizeEV = 10.;
//
//      MonteCarloSST* monte = new MonteCarloSST(eg, chamber, eg->createElectron()); printf("47");
//
//      printf("\n# Trajectories at each landing position: %d", nTrajectories);
//      printf("\n# Pitch of lines (nm): %.10e", pitchnm);
//      printf("\n# lines: %d", nlines);
//      printf("\nLine height (nm): %.10e", hnm);
//      printf("\nLine bottom width (nm): %.10e", wnm);
//      printf("\nLine length (nm): %.10e", linelengthnm);
//      printf("\nLeft and right sidewall angles (deg): %.10e %.10e", thetaldeg, thetardeg);
//      printf("\nLeft and right top corner radii (nm): %.10e %.10e", radlnm, radrnm);
//      printf("\nThicknesses of 1st and second layers (nm): %.10e %.10e", layer1thicknessnm, layer2thicknessnm);
//      printf("\nBeam landing energies (eV): ");
//
//      for (int i = 0; i < beamEeVvalsLen; i++) {
//         printf("\n%.10e", beamEeVvals[i]);
//      }
//      printf("\nBeam size (standard deviation, in nm): %.10e", beamsizenm);
//
//      printf("\n");
//      printf("\nbeamE (eV)\t x(nm)\t y (nm)\t BSE yield\t SE yield");
//
//      for (auto ynm : yvals) {
//         //double ynm = yvals[0];
//         const double y = ynm*meterspernm;
//         for (auto xnm : xvals) {
//            x = xnm*meterspernm;
//            const double egCenter[] = { x, y, -h - 20.*meterspernm };
//            eg->setCenter(egCenter);
//
//            const int nbins = (int)(beamEeV / binSizeEV);
//            BackscatterStatsT* back = new BackscatterStatsT(*monte, nbins); printf("48");
//            monte->addActionListener(*back);
//
//            monte->runMultipleTrajectories(nTrajectories);
//
//            const HistogramT& hist = back->backscatterEnergyHistogram(); printf("49");
//
//            const double energyperbineV = beamEeV / hist.binCount();
//            const double maxSEbin = 50. / energyperbineV;
//            int totalSE = 0;
//            for (int j = 0; j < (int)maxSEbin; ++j) {
//               totalSE = totalSE + hist.counts(j);
//            }
//
//            double SEf = (float)totalSE / nTrajectories;
//            double bsf = back->backscatterFraction() - SEf;
//            printf("\n %.10e %.10e %.10e %.10e %.10e", beamEeV, xnm, ynm, bsf, SEf);
//
//            monte->removeActionListener(*back);
//         }
//      }
//      printf("\n");
//   }
//}
