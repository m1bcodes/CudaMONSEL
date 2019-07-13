#include "gov\nist\microanalysis\EPQLibrary\CzyzewskiMottCrossSection.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "Amphibian\Algorithm.cuh"

#include <sstream>

namespace CzyzewskiMottCrossSection
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const int SpecialEnergyCount = 26;
   __constant__ const int SpecialAngleCount = 96;
   __constant__ const double MaxEnergy = 4.8065296e-15;

   __constant__ static const int kTotalIndex = 96;
   __constant__ static const int kMeanFreePathIndex = 97;

   __constant__ static const double kEnergy[] = {
      3.2043530600e-018,
      8.0108826500e-018,
      1.2016323975e-017,
      1.6021765300e-017,
      3.2043530600e-017,
      4.8065295900e-017,
      6.4087061200e-017,
      8.0108826500e-017,
      9.6130591800e-017,
      1.1215235710e-016,
      1.2817412240e-016,
      1.4419588770e-016,
      1.6021765300e-016,
      3.2043530600e-016,
      4.8065295900e-016,
      6.4087061200e-016,
      8.0108826500e-016,
      9.6130591800e-016,
      1.1215235710e-015,
      1.2817412240e-015,
      1.4419588770e-015,
      1.6021765300e-015,
      2.4032647950e-015,
      3.2043530600e-015,
      4.0054413250e-015,
      4.8065295900e-015
   };

   __constant__ static const double kAngle[] = {
      0.0000000000e+000,
      1.7453292520e-002,
      3.4906585040e-002,
      5.2359877560e-002,
      6.9813170080e-002,
      8.7266462600e-002,
      1.0471975512e-001,
      1.2217304764e-001,
      1.3962634016e-001,
      1.5707963268e-001,
      1.7453292520e-001,
      2.0943951024e-001,
      2.4434609528e-001,
      2.7925268032e-001,
      3.1415926536e-001,
      3.4906585040e-001,
      3.8397243544e-001,
      4.1887902048e-001,
      4.5378560552e-001,
      4.8869219056e-001,
      5.2359877560e-001,
      5.5850536064e-001,
      5.9341194568e-001,
      6.2831853072e-001,
      6.6322511576e-001,
      6.9813170080e-001,
      7.3303828584e-001,
      7.6794487088e-001,
      8.0285145592e-001,
      8.3775804096e-001,
      8.7266462600e-001,
      9.0757121104e-001,
      9.4247779608e-001,
      9.7738438112e-001,
      1.0122909662e+000,
      1.0471975512e+000,
      1.0821041362e+000,
      1.1170107213e+000,
      1.1519173063e+000,
      1.1868238914e+000,
      1.2217304764e+000,
      1.2566370614e+000,
      1.2915436465e+000,
      1.3264502315e+000,
      1.3613568166e+000,
      1.3962634016e+000,
      1.4311699866e+000,
      1.4660765717e+000,
      1.5009831567e+000,
      1.5358897418e+000,
      1.5707963268e+000,
      1.6057029118e+000,
      1.6406094969e+000,
      1.6755160819e+000,
      1.7104226670e+000,
      1.7453292520e+000,
      1.7802358370e+000,
      1.8151424221e+000,
      1.8500490071e+000,
      1.8849555922e+000,
      1.9198621772e+000,
      1.9547687622e+000,
      1.9896753473e+000,
      2.0245819323e+000,
      2.0594885174e+000,
      2.0943951024e+000,
      2.1293016874e+000,
      2.1642082725e+000,
      2.1991148575e+000,
      2.2340214426e+000,
      2.2689280276e+000,
      2.3038346126e+000,
      2.3387411977e+000,
      2.3736477827e+000,
      2.4085543678e+000,
      2.4434609528e+000,
      2.4783675378e+000,
      2.5132741229e+000,
      2.5481807079e+000,
      2.5830872930e+000,
      2.6179938780e+000,
      2.6529004630e+000,
      2.6878070481e+000,
      2.7227136331e+000,
      2.7576202182e+000,
      2.7925268032e+000,
      2.8274333882e+000,
      2.8623399733e+000,
      2.8972465583e+000,
      2.9321531434e+000,
      2.9670597284e+000,
      3.0019663134e+000,
      3.0368728985e+000,
      3.0717794835e+000,
      3.1066860685e+000,
      3.1415926536e+000
   };
#else
   const int SpecialEnergyCount = 26;
   const int SpecialAngleCount = 96;
   const double MaxEnergy = ToSI::keV(30.0);

   static const int kTotalIndex = 96;
   static const int kMeanFreePathIndex = 97;

   static const double kEnergy[] = {
      ToSI::eV(20),
      ToSI::eV(50),
      ToSI::eV(75),
      ToSI::eV(100),
      ToSI::eV(200),
      ToSI::eV(300),
      ToSI::eV(400),
      ToSI::eV(500),
      ToSI::eV(600),
      ToSI::eV(700),
      ToSI::eV(800),
      ToSI::eV(900),
      ToSI::keV(1),
      ToSI::keV(2),
      ToSI::keV(3),
      ToSI::keV(4),
      ToSI::keV(5),
      ToSI::keV(6),
      ToSI::keV(7),
      ToSI::keV(8),
      ToSI::keV(9),
      ToSI::keV(10),
      ToSI::keV(15),
      ToSI::keV(20),
      ToSI::keV(25),
      ToSI::keV(30)
   };

   static const double kAngle[] = {
      Math2::toRadians(0.0),
      Math2::toRadians(1),
      Math2::toRadians(2),
      Math2::toRadians(3),
      Math2::toRadians(4),
      Math2::toRadians(5),
      Math2::toRadians(6),
      Math2::toRadians(7),
      Math2::toRadians(8),
      Math2::toRadians(9),
      Math2::toRadians(10),
      Math2::toRadians(12),
      Math2::toRadians(14),
      Math2::toRadians(16),
      Math2::toRadians(18),
      Math2::toRadians(20),
      Math2::toRadians(22),
      Math2::toRadians(24),
      Math2::toRadians(26),
      Math2::toRadians(28),
      Math2::toRadians(30),
      Math2::toRadians(32),
      Math2::toRadians(34),
      Math2::toRadians(36),
      Math2::toRadians(38),
      Math2::toRadians(40),
      Math2::toRadians(42),
      Math2::toRadians(44),
      Math2::toRadians(46),
      Math2::toRadians(48),
      Math2::toRadians(50),
      Math2::toRadians(52),
      Math2::toRadians(54),
      Math2::toRadians(56),
      Math2::toRadians(58),
      Math2::toRadians(60),
      Math2::toRadians(62),
      Math2::toRadians(64),
      Math2::toRadians(66),
      Math2::toRadians(68),
      Math2::toRadians(70),
      Math2::toRadians(72),
      Math2::toRadians(74),
      Math2::toRadians(76),
      Math2::toRadians(78),
      Math2::toRadians(80),
      Math2::toRadians(82),
      Math2::toRadians(84),
      Math2::toRadians(86),
      Math2::toRadians(88),
      Math2::toRadians(90),
      Math2::toRadians(92),
      Math2::toRadians(94),
      Math2::toRadians(96),
      Math2::toRadians(98),
      Math2::toRadians(100),
      Math2::toRadians(102),
      Math2::toRadians(104),
      Math2::toRadians(106),
      Math2::toRadians(108),
      Math2::toRadians(110),
      Math2::toRadians(112),
      Math2::toRadians(114),
      Math2::toRadians(116),
      Math2::toRadians(118),
      Math2::toRadians(120),
      Math2::toRadians(122),
      Math2::toRadians(124),
      Math2::toRadians(126),
      Math2::toRadians(128),
      Math2::toRadians(130),
      Math2::toRadians(132),
      Math2::toRadians(134),
      Math2::toRadians(136),
      Math2::toRadians(138),
      Math2::toRadians(140),
      Math2::toRadians(142),
      Math2::toRadians(144),
      Math2::toRadians(146),
      Math2::toRadians(148),
      Math2::toRadians(150),
      Math2::toRadians(152),
      Math2::toRadians(154),
      Math2::toRadians(156),
      Math2::toRadians(158),
      Math2::toRadians(160),
      Math2::toRadians(162),
      Math2::toRadians(164),
      Math2::toRadians(166),
      Math2::toRadians(168),
      Math2::toRadians(170),
      Math2::toRadians(172),
      Math2::toRadians(174),
      Math2::toRadians(176),
      Math2::toRadians(178),
      Math2::toRadians(180)
   };
#endif

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

   void CzyzewskiMottCrossSection::loadTables(int atomicNo)
   {
      std::string name(".\\gov\\nist\\microanalysis\\EPQLibrary\\CzyzewskiXSec/" + (atomicNo < 10 ? "0" + std::to_string(atomicNo) : std::to_string(atomicNo)) + ".dat");
      printf("Reading: %s\n", name.c_str());
      try {
         std::ifstream t(name);
         if (!t.good()) throw 0;
         std::string line;
         const char delim[4] = "  ";
         char* tok = NULL;
         char* next_token = NULL;

         mValues.resize(SpecialEnergyCount, VectorXd(kMeanFreePathIndex + 1, 0));
         std::getline(t, line);
         tok = strtok_s((char*)line.c_str(), delim, &next_token);
         // ensure
         // that
         // numbers are
         // parsed
         // correctly
         for (int r = 0; r < SpecialEnergyCount; ++r) {
            for (int c = 0; c <= kMeanFreePathIndex; ++c) {
               if (!tok) {
                  std::getline(t, line);
                  tok = strtok_s((char*)line.c_str(), delim, &next_token);
               }
               mValues[r][c] = (float)sciToDub(tok);
               tok = strtok_s(NULL, delim, &next_token);
            }
         }
         t.close();
      }
      catch (std::exception ex) {
         //throw new EPQFatalException("Fatal error loading a Mott cross section data file. - " + name);
         printf("Fatal error loading a Mott cross section data file. - %s\n", name.c_str());
      }
   }

   CzyzewskiMottCrossSection::CzyzewskiMottCrossSection(const ElementT& el) : mElement(el)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
#else
      loadTables(el.getAtomicNumber());
#endif
   }

   CzyzewskiMottCrossSection::CzyzewskiMottCrossSection(int an) : mElement(Element::byAtomicNumber(an))
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
#else
      loadTables(an);
#endif

   }

   StringT CzyzewskiMottCrossSection::toString() const
   {
      return "CzyzewskiMott[" + StringT(mElement.toAbbrev()) + "]";
   }

   __host__ __device__ const ElementT& CzyzewskiMottCrossSection::getElement() const
   {
      return mElement;
   }

   __host__ __device__ double getSpecialEnergy(const int index)
   {
      return kEnergy[index];
   }

   __host__ __device__ double getSpecialAngle(const int index)
   {
      return kAngle[index];
   }

   __host__ __device__ int getEnergyIndex(const double energy)
   {
      const int ei = Algorithm::binarySearch(kEnergy, 0, sizeof(kEnergy)/sizeof(double), energy);
      return ei >= 0 ? ei : -(ei + 1);
   }

   double CzyzewskiMottCrossSection::totalCrossSection(const double energy) const
   {
      int ei = Algorithm::binarySearch(kEnergy, 0, sizeof(kEnergy) / sizeof(double), energy);
      if (ei >= 0)
         return ToSI::sqrAngstrom(mValues[ei][kTotalIndex]);
      else {
         ei = -(ei + 1);
         return ToSI::sqrAngstrom(mValues[ei - 1][kTotalIndex] + (energy - kEnergy[ei - 1]) / (kEnergy[ei] - kEnergy[ei - 1]) * (mValues[ei - 1][kTotalIndex] - mValues[ei][kTotalIndex]));
      }
   }

   double CzyzewskiMottCrossSection::partialCrossSection(double elevation, double azimuth, double energy) const
   {
      int ei = Algorithm::binarySearch(kEnergy, 0, sizeof(kEnergy) / sizeof(double), energy);
      if (ei < 0)
         ei = -(ei + 1);
      if (!(energy <= getSpecialEnergy(ei))) printf("CzyzewskiMottCrossSection::partialCrossSection 1: %lf, %lf", energy, getSpecialEnergy(ei));
      if (!((ei == 0) || (energy >= getSpecialEnergy(ei - 1)))) printf("CzyzewskiMottCrossSection::partialCrossSection 2: ", energy, getSpecialEnergy(ei - 1));
      if (ei == 0)
         ei = 1; // Extrapolate across the boundary...

      int ai = Algorithm::binarySearch(kAngle, 0, sizeof(kAngle) / sizeof(double), elevation);
      if (ai < 0)
         ai = -(ai + 1);
      if (!(elevation <= getSpecialAngle(ai))) printf("CzyzewskiMottCrossSection::partialCrossSection 3: %lf, %lf\n", elevation, getSpecialAngle(ai));
      if (!((ai == 0) || (elevation >= getSpecialAngle(ai - 1)))) printf("CzyzewskiMottCrossSection::partialCrossSection 4: %lf, %lf\n", elevation, getSpecialAngle(ai-1));
      if (ai == 0)
         ai = 1;
      const double t = (energy - kEnergy[ei - 1]) / (kEnergy[ei] - kEnergy[ei - 1]);
      const double u = (elevation - kAngle[ai - 1]) / (kAngle[ai] - kAngle[ai - 1]);
      return ToSI::sqrAngstrom((1.0 - t) * (1.0 - u) * mValues[ei - 1][ai - 1] + t * (1.0 - u) * mValues[ei][ai - 1] + t * u * mValues[ei][ai] + (1.0 - t) * u * mValues[ei - 1][ai]);
   }

   double CzyzewskiMottCrossSection::partialCrossSection(double elevation, double energy) const
   {
      return partialCrossSection(elevation, 0, energy) * 2.0 * Math2::PI * ::sin(elevation);
   }

   double CzyzewskiMottCrossSection::meanFreePath(double energy) const
   {
      int ei = Algorithm::binarySearch(kEnergy, 0, sizeof(kEnergy) / sizeof(double), energy);
      if (ei >= 0)
         return ToSI::angstrom(mValues[ei][kTotalIndex]);
      else {
         ei = -(ei + 1);
         return ToSI::angstrom(mValues[ei - 1][kMeanFreePathIndex] + (energy - kEnergy[ei - 1]) / (kEnergy[ei] - kEnergy[ei - 1]) * (mValues[ei - 1][kMeanFreePathIndex] - mValues[ei][kMeanFreePathIndex]));
      }
   }
}
