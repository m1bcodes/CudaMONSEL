#include "gov\nist\microanalysis\EPQLibrary\MeanIonizationPotential.cuh"
#include "gov\nist\microanalysis\EPQLibrary\CaveatBase.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Composition.cuh"
#include "gov\nist\microanalysis\EPQLibrary\FromSI.cuh"

namespace MeanIonizationPotential
{
   MeanIonizationPotential::MeanIonizationPotential(StringT name, const ReferenceT& reference) : AlgorithmClass("Mean Ionization Potential", name, reference) {
   }

   void MeanIonizationPotential::initializeDefaultStrategy() {}

   StringT caveat(const ElementT& el)
   {
      return CaveatBase::None;
   }

   StringT caveat(const CompositionT& comp)
   {
      StringT res(CaveatBase::None);
      for (auto el : comp.getElementSet())
         res = CaveatBase::append(res, caveat(*el));
      return res;
   }

   double MeanIonizationPotential::computeLn(const CompositionT& comp) const
   {
      double m = 0.0;
      double lnJ = 0.0;
      for (auto el : comp.getElementSet()) {
         double cz_a = comp.weightFraction(*el, true) * el->getAtomicNumber() / el->getAtomicWeight();
         m += cz_a;
         lnJ += cz_a * ::log(FromSI::keV(compute(*el)));
      }
      return ToSI::keV(::exp(lnJ / m));
   }

   Reference::CrudeReference SternheimerCR("Sternheimer quoted in Berger MJ, Seltzer S. NASA Technical Publication SP-4012 (1964)");
   Sternheimer64MeanIonizationPotential::Sternheimer64MeanIonizationPotential() : MeanIonizationPotential("Sternheimer 1964", SternheimerCR)
   {
   }

   double Sternheimer64MeanIonizationPotential::compute(const ElementT& el) const
   {
      double z = el.getAtomicNumber();
      return ToSI::eV(9.76 * z + 58.8 * ::pow(z, -0.19));
   }

   const Sternheimer64MeanIonizationPotential Sternheimer64Ref;
   const MeanIonizationPotential& Sternheimer64 = Sternheimer64Ref;

   Reference::CrudeReference BergerSeltzerCR("Berger and Seltzer as implemented by CITZAF 3.06");
   BergerAndSeltzerCITZAFMeanIonizationPotential::BergerAndSeltzerCITZAFMeanIonizationPotential() : MeanIonizationPotential("Berger & Seltzer as per JTA", BergerSeltzerCR)
   {
   }

   double BergerAndSeltzerCITZAFMeanIonizationPotential::compute(const ElementT& el) const
   {
      double z = el.getAtomicNumber();
      return ToSI::eV(9.76 * z + 58.5 * ::pow(z, -0.19));
   }

   const BergerAndSeltzerCITZAFMeanIonizationPotential BergerAndSeltzerCITZAFRef;
   const MeanIonizationPotential& BergerAndSeltzerCITZAF = BergerAndSeltzerCITZAFRef;

   Reference::CrudeReference Bloch33CR("Bloch F, F. Z. Phys. 81, 363 (1933)");
   Bloch33MeanIonizationPotential::Bloch33MeanIonizationPotential() : MeanIonizationPotential("Bloch 1933", Bloch33CR)
   {
   }

   double Bloch33MeanIonizationPotential::compute(const ElementT& el) const
   {
      double z = el.getAtomicNumber();
      return ToSI::eV(13.5 * z);
   }

   const Bloch33MeanIonizationPotential Bloch33Ref;
   const MeanIonizationPotential& Bloch33 = Bloch33Ref;

   Reference::CrudeReference Wilson41CR("Wilson RR. Phys Rev. 60. 749 (1941)");
   Wilson41MeanIonizationPotential::Wilson41MeanIonizationPotential() : MeanIonizationPotential("Wilson 1941", Wilson41CR)
   {
   }

   double Wilson41MeanIonizationPotential::compute(const ElementT& el) const
   {
      double z = el.getAtomicNumber();
      return ToSI::eV(11.5 * z);
   }

   const Wilson41MeanIonizationPotential Wilson41Ref;
   const MeanIonizationPotential& Wilson41 = Wilson41Ref;

   Reference::CrudeReference Springer67CR("Springer G. Meues Jahrbuch Fuer Mineralogie, Monatshefte (1967) 9/10, p. 304");
   Springer67MeanIonizationPotential::Springer67MeanIonizationPotential() : MeanIonizationPotential("Springer 1967", Springer67CR)
   {
   }

   double Springer67MeanIonizationPotential::compute(const ElementT& el) const
   {
      double z = el.getAtomicNumber();
      return ToSI::eV(z * (9.0 * (1.0 + ::pow(z, -0.67)) + 0.03 * z));
   }

   const Springer67MeanIonizationPotential Springer67Ref;
   const MeanIonizationPotential& Springer67 = Springer67Ref;

   Reference::CrudeReference Heinrich70CR("Heinrich KFJ, Yakowitz H. Mikrochim Acta (1970) p 123");
   Heinrich70MeanIonizationPotential::Heinrich70MeanIonizationPotential() : MeanIonizationPotential("Heinrich & Yakowitz 1970", Heinrich70CR)
   {
   }

   double Heinrich70MeanIonizationPotential::compute(const ElementT& el) const
   {
      double z = el.getAtomicNumber();
      return ToSI::eV(z * (12.4 + 0.027 * z));
   }

   const Heinrich70MeanIonizationPotential Heinrich70Ref;
   const MeanIonizationPotential& Heinrich70 = Heinrich70Ref;

   Reference::CrudeReference Duncumb69CR("Duncumb P, Shields-Mason PK, DeCasa C. Proc. 5th Int. Congr. on X-ray Optics and Microanalysis, Springer, Berlin, 1969 p. 146");
   Duncumb69MeanIonizationPotential::Duncumb69MeanIonizationPotential() : MeanIonizationPotential("Duncumb & DeCasa 1969", Duncumb69CR)
   {
   }

   double Duncumb69MeanIonizationPotential::compute(const ElementT& el) const
   {
      double z = el.getAtomicNumber();
      return ToSI::eV((14.0 * (1.0 - ::exp(-0.1 * z)) + 75.5 / ::pow(z, z / 7.5) - z / (100 + z)) * z);
   }

   const Duncumb69MeanIonizationPotential Duncumb69Ref;
   const MeanIonizationPotential& Duncumb69 = Duncumb69Ref;

   Reference::CrudeReference Zeller75CR("Zeller C in Ruste J, Gantois M, J. Phys. D. Appl. Phys 8, 872 (1975)");
   Zeller75MeanIonizationPotential::Zeller75MeanIonizationPotential() : MeanIonizationPotential("Zeller 1975", Zeller75CR)
   {
   }

   double Zeller75MeanIonizationPotential::compute(const ElementT& el) const
   {
      double z = el.getAtomicNumber();
      return ToSI::eV((10.04 + 8.25 * ::exp(-z / 11.22)) * z);
   }

   const Zeller75MeanIonizationPotential Zeller75Ref;
   const MeanIonizationPotential& Zeller75 = Zeller75Ref;

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
   Berger64MeanIonizationPotential::Berger64MeanIonizationPotential() : MeanIonizationPotential("Berger & Seltzer 1964", Berger64CR)
   {
   }

   VectorXd Berger64MeanIonizationPotential::mMeasured; // nominal, in Joules
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

   double Berger64MeanIonizationPotential::compute(const ElementT& el) const
   {
      return mMeasured[el.getAtomicNumber() - 1];
   }

   const Berger64MeanIonizationPotential Berger64Ref;
   const MeanIonizationPotential& Berger64 = Berger64Ref;

   Reference::CrudeReference Berger83CR("Berger MJ, Seltzer S. NBSIR 82-2550-A - US Dept of Commerce, Washington DC (1983)");
   Berger83MeanIonizationPotential::Berger83MeanIonizationPotential() : MeanIonizationPotential("Berger & Seltzer 1983", Berger83CR)
   {
   }

   VectorXd Berger83MeanIonizationPotential::mMeasured; // nominal, in Joules
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

   double Berger83MeanIonizationPotential::compute(const ElementT& el) const
   {
      return mMeasured[el.getAtomicNumber() - 1];
   }

   const Berger83MeanIonizationPotential Berger83Ref;
   const MeanIonizationPotential& Berger83 = Berger83Ref;

   const AlgorithmClassT *  mAllImplementations[] = {
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

   AlgorithmClassT const * const * MeanIonizationPotential::getAllImplementations() const
   {
      return mAllImplementations;
   }
}
