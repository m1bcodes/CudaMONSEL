#include "gov\nist\microanalysis\EPQLibrary\EdgeEnergy.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Reference.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AtomicShell.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

#include "gov\nist\microanalysis\Utility\CSVReader.h"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

namespace EdgeEnergy
{
   static const Reference::CrudeReference sCrudeRef("Edge Energy");

   __host__ __device__ void EdgeEnergy::initializeDefaultStrategy()
   {
   }

   EdgeEnergy::EdgeEnergy(StringT name, StringT ref) : AlgorithmClassT(name, ref, sCrudeRef)
   {
   }

   EdgeEnergy::EdgeEnergy(StringT name, const ReferenceT& ref) : AlgorithmClassT("Edge Energy", name, ref)
   {
   }

   //float compute(XRayTransition xrt)
   //{
   //   return compute(xrt.getDestination());
   //}

   // https://www.oreilly.com/library/view/c-cookbook/0596007612/ch03s06.html
   static float sciToDub(const std::string& str)
   {
      std::stringstream ss(str);
      float d = 0;
      ss >> d;

      if (ss.fail()) {
         std::string s = "Unable to format ";
         s += str;
         s += " as a number!";
         throw (s);
      }

      return d;
   }

   MatrixXd DiracHartreeSlaterIonizationEnergies::Uis;
   void DiracHartreeSlaterIonizationEnergies::loadxionUis()
   {
      if (!Uis.empty()) return;

      std::string name(".\\gov\\nist\\microanalysis\\EPQLibrary\\SalvatXion/xionUis.csv");
      printf("Reading: %s\n", name.c_str());
      try {
         std::ifstream file(name);
         if (!file.good()) throw 0;
         Uis.resize(100);
         int idx = 0;
         for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
            int z = std::stoi((*loop)[0]);
            if (z != idx + 1) printf("loadxionUis: wrong line %d\n", z);
            for (int k = 1; k < (*loop).size(); ++k) {
               float cur = sciToDub((*loop)[k]);
               Uis[z].push_back(cur);
            }
            ++idx;
         }
         file.close();
      }
      catch (const std::exception&) {
         printf("cannot load loadxionUis: %s\n", name.c_str());
      }
   }

   static const Reference::Author* al[] = { &Reference::Bote, &Reference::FSalvat };
   static const Reference::JournalArticle ja(Reference::PhysRevA, "77", "042701-1 to 24", 2008, al, 2);

   DiracHartreeSlaterIonizationEnergies::DiracHartreeSlaterIonizationEnergies() : EdgeEnergy("Bote-Salvat 2008", ja)
   {
   }

   float DiracHartreeSlaterIonizationEnergies::compute(const AtomicShellT& shell) const
   {
      if (Uis.empty()) printf("DTSAEdgeEnergy::compute: DTSAEdgeEnergies.empty()\n");
      return ToSI::eV(Uis[shell.getElement().getAtomicNumber()][shell.getShell()]);
   }

   bool DiracHartreeSlaterIonizationEnergies::isSupported(const AtomicShellT& shell) const
   {
      if (Uis.empty()) printf("DTSAEdgeEnergy::isSupported: DTSAEdgeEnergies.empty()\n");
      int z = shell.getElement().getAtomicNumber();
      int sh = shell.getShell();
      return (z >= 1) && (z <= 99) && (sh < Uis[z].size());
   };

   static const DiracHartreeSlaterIonizationEnergies DHSIonizationEnergyRef;
   const EdgeEnergy& DHSIonizationEnergy = DHSIonizationEnergyRef;

   MatrixXd NISTEdgeEnergy::NISTEdgeEnergies;
   void NISTEdgeEnergy::loadNISTxrtdb()
   {
      if (!NISTEdgeEnergies.empty()) return;

      std::string name(".\\gov\\nist\\microanalysis\\EPQLibrary\\NISTxrtdb.csv");
      printf("Reading: %s\n", name.c_str());
      try {
         std::ifstream file(name);
         if (!file.good()) throw 0;
         NISTEdgeEnergies.resize(91, VectorXd(4, 0));
         int idx = 0;
         for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
            for (int k = 0; k < (*loop).size(); ++k) {
               float cur = sciToDub((*loop)[k]);
               NISTEdgeEnergies[idx][k] = cur;
            }
            ++idx;
         }
         file.close();
      }
      catch (const std::exception&) {
         printf("cannot load loadNISTxrtdb: %s\n", name.c_str());
      }
   }

   NISTEdgeEnergy::NISTEdgeEnergy() : EdgeEnergy("NIST X-ray transition database", "http://physics.nist.gov/PhysRefData/XrayTrans/")
   {
   }

   float NISTEdgeEnergy::compute(const AtomicShellT& shell) const
   {
      if (NISTEdgeEnergies.empty()) printf("DTSAEdgeEnergy::compute: DTSAEdgeEnergies.empty()\n");
      int an = shell.getElement().getAtomicNumber();
      if ((an < Element::elmNe) || (an > Element::elmFm)) printf("NISTEdgeEnergy::compute: %d - only supports elements Ne (10) to Fm (100)\n", an);
      int sh = shell.getShell();
      if ((sh < AtomicShell::K) || (sh > AtomicShell::LIII)) printf("NISTEdgeEnergy::compute: %d - only supports shells K, L1, L2 and L3\n", sh);
      return ToSI::eV(NISTEdgeEnergies[an - Element::elmNe][sh - AtomicShell::K]);
   }

   bool NISTEdgeEnergy::isSupported(const AtomicShellT& shell) const
   {
      if (NISTEdgeEnergies.empty()) printf("DTSAEdgeEnergy::isSupported: DTSAEdgeEnergies.empty()\n");
      int an = shell.getElement().getAtomicNumber();
      int sh = shell.getShell();
      return (an >= Element::elmNe) && (an <= Element::elmFm) && (sh >= AtomicShell::K) && (sh <= AtomicShell::LIII) && (NISTEdgeEnergies[an - Element::elmNe][sh - AtomicShell::K] > 0.0);
   }

   static const NISTEdgeEnergy NISTxrtdbRef;
   const EdgeEnergy& NISTxrtdb = NISTxrtdbRef;

   MatrixXd ChantlerEdgeEnergy::ChantlerEdgeEnergies;
   void ChantlerEdgeEnergy::loadFFastEdgeDB()
   {
      if (!ChantlerEdgeEnergies.empty()) return;

      std::string name(".\\gov\\nist\\microanalysis\\EPQLibrary\\FFastEdgeDB.csv");
      printf("Reading: %s\n", name.c_str());
      try {
         std::ifstream file(name);
         if (!file.good()) throw 0;
         ChantlerEdgeEnergies.resize(92);
         int idx = 0;
         for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
            ChantlerEdgeEnergies[idx].reserve((*loop).size());
            for (int k = 0; k < (*loop).size(); ++k) {
               float cur = sciToDub((*loop)[k]);
               if (cur) ChantlerEdgeEnergies[idx].push_back(ToSI::eV(cur));
            }
            ++idx;
         }
         file.close();
      }
      catch (const std::exception&) {
         printf("cannot load loadFFastEdgeDB: %s\n", name.c_str());
      }
   }

   ChantlerEdgeEnergy::ChantlerEdgeEnergy() : EdgeEnergy("NIST-Chantler 2005", "http://physics.nist.gov/ffast")
   {
   }

   static int index(int sh)
   {
      if ((sh < AtomicShell::K) || (sh > AtomicShell::PIII))
         return -1;
      if (sh <= AtomicShell::OV)
         return sh;
      if (sh >= AtomicShell::PI)
         return sh - 4; // (AtomicShell::PI-AtomicShell::OV + 1);
      return -1;
   }

   float ChantlerEdgeEnergy::compute(const AtomicShellT& shell) const
   {
      if (ChantlerEdgeEnergies.empty()) printf("ChantlerEdgeEnergy::compute ChantlerEdgeEnergies empty!\n");
      int an = shell.getElement().getAtomicNumber();
      if ((an < Element::elmH) || (an > Element::elmU)) printf("ChantlerEdgeEnergy::compute: %d - only supports elements H (1) to U (92)\n", an);
      int sh = shell.getShell();
      int i = index(sh);
      if (i == -1) printf("ChantlerEdgeEnergy::compute: %d - only supports shells K to O5, P1 to P3\n", i);
      if (!(ChantlerEdgeEnergies.size() > an - Element::elmH)) printf("Too few elements in EdgeEnergy database.\n");
      return ChantlerEdgeEnergies[an - Element::elmH].size() > i ? ChantlerEdgeEnergies[an - Element::elmH][i] : 0.0;
   }

   bool ChantlerEdgeEnergy::isSupported(const AtomicShellT& shell) const
   {
      if (ChantlerEdgeEnergies.empty()) printf("ChantlerEdgeEnergy::isSupported ChantlerEdgeEnergies empty!\n");
      int an = shell.getElement().getAtomicNumber();
      int i = index(shell.getShell());
      return (an >= Element::elmLi) && (an <= Element::elmU) && (i >= 0) && (ChantlerEdgeEnergies[an - Element::elmH].size() > i) && (ChantlerEdgeEnergies[an - Element::elmH][i] > 0.0);
   }

   static const ChantlerEdgeEnergy Chantler2005Ref;
   const EdgeEnergy& Chantler2005 = Chantler2005Ref;

   WernishEdgeEnergy::WernishEdgeEnergy() : EdgeEnergy("Wernisch et al., 1984", "Wernisch et al., 1984 - Taken from Markowitz in the Handbook of X-ray Spectroscopy") {
   }

   float WernishEdgeEnergy::compute(const AtomicShellT& shell) const
   {
      float z = shell.getElement().getAtomicNumber();
      if (!isSupported(shell)) printf(StringT("Unsupported shell " + StringT(shell.toString()) + " in " + toString() + "\n").c_str());
      if (shell.getShell() == AtomicShell::K) return ToSI::keV(-1.304e-1 + z * (-2.633e-3 + z * (9.718e-3 + z * 4.144e-5)));
      if (shell.getShell() == AtomicShell::LI) return ToSI::keV(-4.506e-1 + z * (1.566e-2 + z * (7.599e-4 + z * 1.792e-5)));
      if (shell.getShell() == AtomicShell::LII) return ToSI::keV(-6.018e-1 + z * (1.964e-2 + z * (5.935e-4 + z * 1.843e-5)));
      if (shell.getShell() == AtomicShell::LIII) return ToSI::keV(3.390e-1 + z * (-4.931e-2 + z * (2.336e-3 + z * 1.836e-6)));
      if (shell.getShell() == AtomicShell::MI) return ToSI::keV(-8.645 + z * (3.977e-1 + z * (-5.963e-3 + z * 3.624e-5)));
      if (shell.getShell() == AtomicShell::MII) return ToSI::keV(-7.499 + z * (3.459e-1 + z * (-5.250e-3 + z * 3.263e-5)));
      if (shell.getShell() == AtomicShell::MIII) return ToSI::keV(-6.280 + z * (2.831e-1 + z * (-4.117e-3 + z * 2.505e-5)));
      if (shell.getShell() == AtomicShell::MIV) return ToSI::keV(-4.778 + z * (2.184e-1 + z * (-3.303e-3 + z * 2.115e-5)));
      if (shell.getShell() == AtomicShell::MV) return ToSI::keV(-2.421 + z * (1.172e-1 + z * (-1.845e-3 + z * 1.397e-5)));
      printf(StringT("Unsupported shell in " + toString()).c_str());
      return -1;
   }

   bool WernishEdgeEnergy::isSupported(const AtomicShellT& shell) const
   {
      int z = shell.getElement().getAtomicNumber();
      if (shell.getShell() == AtomicShell::K) return ((z >= 11) && (z <= 63));
      if (shell.getShell() == AtomicShell::LI) return ((z >= 28) && (z <= 83));
      if (shell.getShell() == AtomicShell::LII) return ((z >= 30) && (z <= 83));
      if (shell.getShell() == AtomicShell::LIII) return ((z >= 30) && (z <= 83));
      if (shell.getShell() == AtomicShell::MI) return ((z >= 52) && (z <= 83));
      if (shell.getShell() == AtomicShell::MII) return ((z >= 55) && (z <= 83));
      if (shell.getShell() == AtomicShell::MIII) return ((z >= 55) && (z <= 83));
      if (shell.getShell() == AtomicShell::MIV) return ((z >= 60) && (z <= 83));
      if (shell.getShell() == AtomicShell::MV) return ((z >= 61) && (z <= 83));
      return false;
   }

   static const WernishEdgeEnergy Wernish84Ref;
   const EdgeEnergy& Wernish84 = Wernish84Ref;

   MatrixXd DTSAEdgeEnergy::DTSAEdgeEnergies; // TODO: fix initialization
   void DTSAEdgeEnergy::loadEdgeEnergies()
   {
      if (!DTSAEdgeEnergies.empty()) return;

      std::string name(".\\gov\\nist\\microanalysis\\EPQLibrary\\EdgeEnergies.csv");
      printf("Reading: %s\n", name.c_str());
      try {
         std::ifstream file(name);
         if (!file.good()) throw 0;
         DTSAEdgeEnergies.resize(99);
         int idx = 0;
         for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
            if ((*loop).size()) {
               DTSAEdgeEnergies[idx].reserve((*loop).size());
               for (int k = 0; k < (*loop).size(); ++k) {
                  float cur = sciToDub((*loop)[k]);
                  if (cur) DTSAEdgeEnergies[idx].push_back(cur);
               }
               ++idx;
            }
         }
         file.close();
      }
      catch (const std::exception&) {
         printf("cannot load loadEdgeEnergies: %s\n", name.c_str());
      }
   }

   DTSAEdgeEnergy::DTSAEdgeEnergy() : EdgeEnergy("DTSA", "From DTSA at http://www.cstl.nist.gov/div837/Division/outputs/DTSA/DTSA.htm")
   {
   }

   float DTSAEdgeEnergy::compute(const AtomicShellT& shell) const
   {
      StringT shellinfo(shell.toString());
      if (DTSAEdgeEnergies.empty()) printf("DTSAEdgeEnergy::compute: DTSAEdgeEnergies.empty()\n");
      int sh = shell.getShell();
      if ((sh < AtomicShell::K) || (sh > AtomicShell::NI)) printf(StringT("Unsupported shell " + shellinfo + " in " + toString() + "\n").c_str());
      int i = shell.getElement().getAtomicNumber() - 1;
      if ((i < 0) || (i >= DTSAEdgeEnergies.size())) printf(StringT("Unsupported element " + shellinfo + " in " + toString() + "\n").c_str());
      return (!DTSAEdgeEnergies[i].empty() && DTSAEdgeEnergies[i].size() > sh) ? ToSI::eV(DTSAEdgeEnergies[i][sh]) : 0.0;
   }

   bool DTSAEdgeEnergy::isSupported(const AtomicShellT& shell) const
   {
      if (DTSAEdgeEnergies.empty()) printf("DTSAEdgeEnergy::isSupported: DTSAEdgeEnergies.empty()\n");
      int sh = shell.getShell();
      int zp = shell.getElement().getAtomicNumber() - 1;
      return (zp < DTSAEdgeEnergies.size() && !DTSAEdgeEnergies[zp].empty()) && (DTSAEdgeEnergies[zp].size() > sh);
   }

   static const DTSAEdgeEnergy DTSARef;
   const EdgeEnergy& DTSA = DTSARef;

   SuperSetEdgeEnergy::SuperSetEdgeEnergy() : EdgeEnergy("Superset", "Chantler then NIST then DTSA")
   {
   }

   float SuperSetEdgeEnergy::compute(const AtomicShellT& shell) const
   {
      try {
         return Chantler2005.compute(shell);
      }
      catch (const std::exception&) {
         try {
            return NISTxrtdb.compute(shell);
         }
         catch (const std::exception&) {
            try {
               return DTSA.compute(shell);
            }
            catch (const std::exception&) {
               return -1.0;
            }
         }
      }
   }

   bool SuperSetEdgeEnergy::isSupported(const AtomicShellT& shell) const
   {
      return Chantler2005.isSupported(shell) || NISTxrtdb.isSupported(shell) || DTSA.isSupported(shell);
   }

   static const SuperSetEdgeEnergy SuperSetRef;
   const EdgeEnergy& SuperSet = SuperSetRef;

   const AlgorithmClassT* mAllImplementations[] = {
      &DTSA,
      &Chantler2005,
      &NISTxrtdb,
      &Wernish84,
      &DHSIonizationEnergy,
      &SuperSet
   };

   AlgorithmClassT const * const * EdgeEnergy::getAllImplementations() const
   {
      return mAllImplementations;
   }
};
