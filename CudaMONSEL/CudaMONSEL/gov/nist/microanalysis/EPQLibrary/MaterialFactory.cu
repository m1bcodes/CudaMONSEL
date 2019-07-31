#include "gov\nist\microanalysis\EPQLibrary\MaterialFactory.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Composition.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"

namespace MaterialFactory
{
   const char* K3189 = "K3189";
   const char* RynasAlTiAlloy = "Ryna's Al-Ti alloy";
   const char* Mylar = "Mylar";
   const char* VanadiumPentoxide = "Vanadium Pentoxide";
   const char* SiliconDioxide = "Silicon Dioxide";
   const char* Ice = "Ice";
   const char* PerfectVacuum = "Perfect vacuum";
   const char* CaCO3 = "Calcium carbonate";
   const char* Al2O3 = "Alumina";
   const char* SS316 = "Stainless Steel 316";
   const char* UraniumOxide = "Uranium oxide";
   const char* K227 = "K227";
   const char* K309 = "K309";
   const char* K411 = "K411";
   const char* K412 = "K412";
   const char* K961 = "K961";
   const char* K1080 = "K1080";
   const char* K2450 = "K2450";
   const char* K2451 = "K2451";
   const char* K2466 = "K2466";
   const char* K2469 = "K2469";
   const char* K2472 = "K2472";
   const char* K2496 = "K2496";
   const char* ParaleneC = "Paralene C";
   const char* MagnesiumOxide = "Magnesium Oxide";
   const char* Chloroapatite = "Chloroapatite";
   const char* CalciumFluoride = "Calcium Fluoride";
   const char* GalliumPhosphate = "Gallium Phosphate";
   const char* Nothing = "None";

   typedef std::unordered_map<StringT, CompositionT*> CompositionMap;

   static const char* PBaseMaterials[] = {
      K3189,
      RynasAlTiAlloy,
      Mylar,
      VanadiumPentoxide,
      SiliconDioxide,
      Ice,
      CaCO3,
      Al2O3,
      SS316,
      UraniumOxide,
      K227,
      K309,
      K411,
      K412,
      K961,
      K1080,
      K2450,
      K2451,
      K2466,
      K2469,
      K2472,
      K2496,
      ParaleneC,
      MagnesiumOxide,
      Chloroapatite,
      CalciumFluoride
   };

   //public static List<String> BaseMaterials = Collections.unmodifiableList(Arrays.asList(PBaseMaterials));

   const ::Material::data_type mElementalDensities[] = {
      0.0f,
      0.0f,
      0.534f,
      1.85f,
      2.53f,
      2.25f,
      0.0f,
      0.0f,
      0.0f,
      0.0f,
      0.97f,
      1.74f,
      2.70f, // H to Ne
      2.42f,
      1.83f, /* Yellow */
      1.92f,
      0.0f,
      0.0f,
      0.86f,
      1.55f,
      3.02f,
      4.5f,
      5.98f,
      7.14f,
      7.41f,
      7.88f,
      8.71f,
      8.88f,
      8.96f,
      7.1f,
      5.93f, // Ga
      5.46f,
      5.73f,
      4.82f,
      0.0f,
      0.0f,
      1.53f,
      2.56f,
      3.8f,
      6.4f,
      8.57f,
      10.22f,
      11.5f,
      12.1f,
      12.44f,
      12.16f,
      10.49f,
      8.65f,
      7.28f,
      7.3f, // Sn
      6.62f,
      6.25f,
      4.94f,
      0.0f,
      1.87f,
      3.5f,
      6.15f,
      6.90f,
      6.48f,
      6.96f,
      -1.0f,
      7.75f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f, // Cs to Dy
      -1.0f,
      13.3f,
      16.6f,
      19.3f,
      -1.0f,
      22.5f,
      22.42f,
      21.45f,
      19.3f,
      14.19f, /* solid at -39 C */
      11.86f,
      11.34f,
      9.78f,
      -1.0f,
      -1.0f, // Ho to Pt
      0.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      11.3f,
      -1.0f,
      18.7f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f, // Fr to Es
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f,
      -1.0f
      // Fm to Uub
   };

   bool canCreate(const ElementT& el)
   {
      return mElementalDensities[el.getAtomicNumber() - 1] > 0.0;
   }

   typedef std::unordered_map <const ElementT*, int> ElementMap; 

   static bool isUpperCase(char c)
   {
      return c >= 65 && c <= 90;
   }

   static bool isDigit(char c)
   {
      return c >= 48&& c <= 57;
   }

   ElementMap parseCompound(std::string str)
   {
      ElementMap elMap;

      std::string elStr;
      char c = (str.length() > 0 ? str.at(0) : INT_MIN);
      for (int i = 0; c != NULL; ++i) {
         c = (i < str.length() ? str.at(i) : NULL);
         if (isUpperCase(c) || isDigit(c) || (c == '(') || (c == NULL)) {
            if (elStr.length() > 0) {
               const ElementT& el = Element::byName(elStr.c_str());
               int count = 1;
               if (!el.isValid()) printf("Unrecognized element %s in %s.", elStr.c_str(), str.c_str());
               elStr = "";
               if (isDigit(c)) {
                  while (isDigit(c)) {
                     elStr += c;
                     ++i;
                     c = (i < str.length() ? str.at(i) : NULL);
                  }
                  count = std::atoi(elStr.c_str());
               }
               int ii = elMap[&el];
               elMap[&el] = ii + count;
               elStr = "";
            }
            while (c == '(') {
               int br = 1;
               for (++i; br > 0; ++i) {
                  c = (i < str.length() ? str.at(i) : NULL);
                  if (c == NULL)
                     printf("Unmatched bracket in %s\n", str.c_str());
                  if (c == ')')
                     --br;
                  if (c == '(')
                     ++br;
                  if (br > 0)
                     elStr += c;
               }
               // Parse bracketed items recursively
               ElementMap subComp = parseCompound(elStr);
               elStr = "";
               c = (i < str.length() ? str.at(i) : NULL);
               int count = 1;
               // Get the count of the bracketed items
               if (isDigit(c)) {
                  while (isDigit(c)) {
                     elStr += c;
                     ++i;
                     c = (i < str.length() ? str.at(i) : NULL);
                  }
                  count = std::atoi(elStr.c_str());
                  elStr = "";
               }
               for (auto e : subComp) {
                  const ElementT* elm = e.first;
                  int ic = elMap[elm];
                  elMap[elm] = ic + count * e.second;
               }
            }
            elStr += c;
         }
         else
            elStr += c;
      }
      // Check for K411 and other anomalies
      if (elMap.size() == 1)
         if (elMap.begin()->second > 1)
            printf("Unusually formated pure element in, ", str.c_str());
      return elMap;
   }

   CompositionT createCompound1(const char* str)
   {
      ElementMap elMap = parseCompound(str);
      CompositionT comp;
      comp.setName(str);
      for (auto e : elMap) comp.addElementByStoiciometry(*e.first, e.second);
      return comp;
   }

   static MaterialT createCompound(StringT str, ::Material::data_type density)
   {
      return MaterialT(createCompound(str), density);
   }

   CompositionT createCompound(StringT str)
   {
      {
         int p = str.find("+");
         if (p != std::string::npos) {
            CompositionT cl = createCompound(str.substr(0, p));
            CompositionT cr = createCompound(str.substr(p + 1));
            const Element::UnorderedSetT& cles = cl.getElementSet();
            const Element::UnorderedSetT& cres = cr.getElementSet();
            std::unordered_set<const ElementT*> elms;
            //elms.insert(cles.begin(), cles.end());
            //elms.insert(cres.begin(), cres.end());
            for (auto e : cles) {
               elms.insert(e);
            }for (auto e : cres) {
               elms.insert(e);
            }

            ::Composition::data_type* massFracs = new ::Composition::data_type[elms.size()];
            std::vector<const ElementT*> elmA(elms.begin(), elms.end());
            for (int i = 0; i < elmA.size(); ++i)
               massFracs[i] = cl.weightFraction(*elmA[i], false) + cr.weightFraction(*elmA[i], false);
            //double* massFracs = new double[elms.size()];
            //for (auto i : elms) {
            //}
            
            CompositionT ret(elmA.data(), elmA.size(), massFracs, elms.size());
            delete[] massFracs;
            return ret;
         }
      }
      {
         int p = str.find("*");
         if (p != std::string::npos) {
            ::Composition::data_type k;
            try {
               k = std::atof(str.substr(0, p).c_str());
            }
            catch (std::exception&) {
               printf("Error parsing number: %s\n", str.substr(0, p).c_str());
            }
            CompositionT cr = createCompound(str.substr(p + 1));
            const Element::UnorderedSetT& elms = cr.getElementSet();
            //std::vector<const ElementT*> elmA(elms.begin(), elms.end());
            std::vector<const ElementT*> elmA;
            for (auto e : elms) {
               elmA.push_back(e);
            }
            std::vector<::Composition::data_type> massFracs(elms.size());
            for (int i = 0; i < elmA.size(); ++i)
               massFracs[i] = k * cr.weightFraction(*elmA[i], false);
            CompositionT ret(elmA.data(), elmA.size(), massFracs.data(), massFracs.size());
            return ret;
         }
      }
      return createCompound1(str.c_str());
   }

   //static const CompositionT compoundSiO2 = createCompound("SiO2");
   //static const CompositionT compoundAl2O3 = createCompound("Al2O3");
   //static const CompositionT compoundCaO = createCompound("CaO");
   //static const CompositionT compoundMgO = createCompound("MgO");
   //static const CompositionT compoundTiO2 = createCompound("TiO2");
   //static const CompositionT compoundFe2O3 = createCompound("Fe2O3");
   //static const CompositionT compoundPbO = createCompound("PbO");
   //static const CompositionT compoundBaO = createCompound("BaO");
   //static const CompositionT compoundFeO = createCompound("FeO");
   //static const CompositionT compoundV2O5 = createCompound("V2O5");

   static CompositionT compoundSiO2;
   static CompositionT compoundAl2O3;
   static CompositionT compoundCaO;
   static CompositionT compoundMgO;
   static CompositionT compoundTiO2;
   static CompositionT compoundFe2O3;
   static CompositionT compoundPbO;
   static CompositionT compoundBaO;
   static CompositionT compoundFeO;
   static CompositionT compoundV2O5;

   void init()
   {
      //StringT compoundSiO2Name("SiO2");
      compoundSiO2 = createCompound("SiO2");
      //StringT compoundAl2O3Name("Al2O3");
      compoundAl2O3 = createCompound("Al2O3");
      //StringT compoundCaOName("CaO");
      compoundCaO = createCompound("CaO");
      //StringT compoundMgOName("MgO");
      compoundMgO = createCompound("MgO");
      //StringT compoundTiO2Name("TiO2");
      compoundTiO2 = createCompound("TiO2");
      //StringT compoundFe2O3Name("Fe2O3");
      compoundFe2O3 = createCompound("Fe2O3");
      //StringT compoundPbOName("PbO");
      compoundPbO = createCompound("PbO");
      //StringT compoundBaOName("BaO");
      compoundBaO = createCompound("BaO");
      //StringT compoundFeOName("FeO");
      compoundFeO = createCompound("FeO");
      //StringT compoundV2O5Name("V2O5");
      compoundV2O5 = createCompound("V2O5");
   }

   CompositionT createMaterial(StringT name)
   {
      try {
         if (name == K3189) {
            const CompositionT* comp[] = { &compoundSiO2, &compoundAl2O3, &compoundCaO, &compoundMgO, &compoundTiO2, &compoundFe2O3 };
            ::Composition::data_type massFracs[] = { 0.400000, 0.140000, 0.140000, 0.100000, 0.020000, 0.200000 };
            CompositionT tmp;
            tmp.defineByMaterialFraction(comp, 6, massFracs, 6);
            MaterialT mat(tmp, ToSI::gPerCC(3.23));
            mat.setName(K3189);
            return mat;
         }
         else if (name == (RynasAlTiAlloy)) {
            const ElementT* elms[] = { &Element::Ti, &Element::Al, &Element::Nb, &Element::W };
            ::Composition::data_type massFracs[] = { 0.54, 0.31, 0.11, 0.04 };
            ::Material::data_type density = ToSI::gPerCC(8.57);
            StringT name = RynasAlTiAlloy;
            return MaterialT(elms, 4, massFracs, 4, density, name.c_str());;
         } 
         else if (name == (Mylar)) {
            MaterialT mat = createCompound("C10H8O4", ToSI::gPerCC(1.39));
            mat.setName(Mylar);
            return mat;
         }
         else if (name == (VanadiumPentoxide)) {
            MaterialT mat = createCompound("V2O5", ToSI::gPerCC(3.357));
            mat.setName(VanadiumPentoxide);
            return mat;
         }
         else if (name == (SiliconDioxide)) {
            MaterialT mat = createCompound("SiO2", ToSI::gPerCC(2.65));
            mat.setName(SiliconDioxide);
            return mat;
         }
         else if (name == (Ice)) {
            MaterialT mat = createCompound("H2O", ToSI::gPerCC(0.917));
            mat.setName(Ice);
            return mat;
         }
         else if (name == (PerfectVacuum)) {
            MaterialT mat(0.0);
            mat.setName(PerfectVacuum);
            return mat;
         }
         else if (name == (CaCO3)) {
            MaterialT mat = createCompound("CaCO3", ToSI::gPerCC(2.7));
            mat.setName(CaCO3);
            return mat;
         }
         else if (name == (Al2O3)) {
            MaterialT mat = createCompound("Al2O3", ToSI::gPerCC(3.97));
            mat.setName(Al2O3);
            return mat;
         }
         else if (name == (SS316)) {
            const ElementT* elms[] = { &Element::Fe, &Element::Ni, &Element::Cr, &Element::Mn, &Element::Si };
            ::Composition::data_type massFracs[] = { 0.50, 0.205, 0.245, 0.02, 0.03 };
            ::Material::data_type density = ToSI::gPerCC(7.8);
            StringT name("SS316");
            MaterialT mat(elms, 5, massFracs, 5, density, name.c_str());
            return mat;
         } else if (name == (UraniumOxide)) {
            MaterialT mat = createCompound("UO2", ToSI::gPerCC(10.0));
            mat.setName(UraniumOxide);
         }
         else if (name == (K227)) {
            const CompositionT* cons[] = { &compoundSiO2, &compoundPbO };
            ::Composition::data_type massFracs[] = { 0.20000, 0.80000, };
            CompositionT comp;
            comp.defineByMaterialFraction(cons, 2, massFracs, 2);
            comp.setName(K227);
            return comp;
         }
         else if (name == (K309)) {
            CompositionT comp;
            const CompositionT* cons[] = { &compoundAl2O3, &compoundSiO2, &compoundCaO, &compoundFe2O3, &compoundBaO };
            ::Composition::data_type massFracs[] = { 0.15000, 0.40000, 0.15000, 0.15000, 0.15000, };
            comp.defineByMaterialFraction(cons, 5, massFracs, 5);
            comp.setName(K309);
            return comp;
         }
         else if (name == (K411)) {
            const CompositionT* cons[] = { &compoundMgO, &compoundSiO2, &compoundCaO, &compoundFeO };
            ::Composition::data_type massFracs[] = { 0.146700, 0.543000, 0.154700, 0.144200 };
            CompositionT comp;
            comp.defineByMaterialFraction(cons, 4, massFracs, 4);
            comp.setName(K411);
            return MaterialT(comp, ToSI::gPerCC(5.0));
         }
         else if (name == (K412)) {
            CompositionT comp;
            const CompositionT* cons[] = { &compoundMgO, &compoundAl2O3, &compoundSiO2, &compoundCaO, &compoundFeO };
            ::Composition::data_type massFracs[] = { 0.193300, 0.092700, 0.453500, 0.152500, 0.099600 };
            comp.defineByMaterialFraction(cons, 5, massFracs, 5);
            comp.setName(K412);
            return comp;
         }
         else if (name == (K961)) {
            const ElementT* elms[] = { &Element::Na, &Element::Mg, &Element::Al, &Element::Si, &Element::P, &Element::K, &Element::Ca, &Element::Ti, &Element::Mn, &Element::Fe, &Element::O };
            ::Composition::data_type massFracs[] = { 0.029674, 0.030154, 0.058215, 0.299178, 0.002182, 0.024904, 0.035735, 0.011990, 0.003160, 0.034972, 0.469837 };
            return MaterialT(elms, 11, massFracs, 11, ToSI::gPerCC(6.0), K961);
         }
         else if (name == (K1080)) {
            const ElementT* elms[] = { &Element::Li, &Element::B, &Element::Mg, &Element::Al, &Element::Si, &Element::Ca, &Element::Ti, &Element::Sr, &Element::Zr, &Element::Lu, &Element::O };
            ::Composition::data_type massFracs[] = { 0.027871, 0.006215, 0.008634, 0.079384, 0.186986, 0.107204, 0.011990, 0.126838, 0.007403, 0.017588, 0.416459 };
            return MaterialT(elms, 11, massFracs, 11, ToSI::gPerCC(6.0), K1080);
         }
         else if (name == (K2450)) {
            const CompositionT* cons[] = { &compoundSiO2, &compoundAl2O3, &compoundCaO, &compoundTiO2 };
            ::Composition::data_type massFracs[] = { 0.30000, 0.30000, 0.30000, 0.10000 };
            CompositionT comp;
            comp.defineByMaterialFraction(cons, 4, massFracs, 4);
            comp.setName(name.c_str());
            return comp;
         }
         else if (name == (K2451)) {
            const CompositionT* cons[] = { &compoundSiO2, &compoundAl2O3, &compoundCaO, &compoundV2O5 };
            ::Composition::data_type massFracs[] = { 0.300000, 0.300000, 0.300000, 0.100000 };
            CompositionT comp;
            comp.defineByMaterialFraction(cons, 4, massFracs, 4);
            comp.setName(name.c_str());
            return comp;
         }
         else if (name == (K2466)) {
            CompositionT comp;
            const CompositionT* cons[] = { &compoundSiO2, &compoundBaO, &compoundTiO2 };
            ::Composition::data_type massFracs[] = { 0.44, 0.48, 0.08 };
            comp.defineByMaterialFraction(cons, 3, massFracs, 3);
            comp.setName(name.c_str());
            return comp;
         }
         else if (name == (K2469)) {
            const CompositionT* cons[] = { &compoundSiO2, &compoundBaO, &compoundTiO2 };
            ::Composition::data_type massFracs[] = { 0.36, 0.48, 0.16 };
            CompositionT comp;
            comp.defineByMaterialFraction(cons, 3, massFracs, 3);
            comp.setName(name.c_str());
            return comp;
         }
         else if (name == (K2472)) {
            const CompositionT* cons[] = { &compoundSiO2, &compoundBaO, &compoundTiO2, &compoundV2O5 };
            ::Composition::data_type massFracs[] = { 0.36, 0.48, 0.10, 0.06 };
            CompositionT comp;
            comp.defineByMaterialFraction(cons, 4, massFracs, 4);
            comp.setName(name.c_str());
            return comp;
         }
         else if (name == (K2496)) {
            CompositionT comp;
            const CompositionT* cons[] = { &compoundSiO2, &compoundBaO, &compoundTiO2 };
            ::Composition::data_type massFracs[] = { 0.49, 0.48, 0.03 };
            comp.defineByMaterialFraction(cons, 3, massFracs, 3);
            comp.setName(name.c_str());
            return comp;
         }
         else if (name == (ParaleneC)) {
            MaterialT mat = createCompound("C8H7Cl", ToSI::gPerCC(1.2));
            mat.setName(ParaleneC);
            return mat;
         }
         else if (name == (MagnesiumOxide)) {
            MaterialT mat = createCompound("MgO", ToSI::gPerCC(3.55));
            mat.setName(MagnesiumOxide);
            return mat;
         }
         else if (name == (CalciumFluoride)) {
            MaterialT mat = createCompound("CaF2", ToSI::gPerCC(3.18));
            mat.setName(CalciumFluoride);
            return mat;
         }
         else if (name == (Chloroapatite)) {
            MaterialT mat = createCompound("Ca5(PO4)3Cl", ToSI::gPerCC(3.15));
            mat.setName(Chloroapatite);
            return mat;
         }
         else if (name == (GalliumPhosphate)) {
            MaterialT mat = createCompound("GaPO4", ToSI::gPerCC(3.570));
            mat.setName(GalliumPhosphate);
            return mat;
         }
         else if (name == (Nothing)) {
            MaterialT mat(0.0);
            mat.setName("None");
            return mat;
         }
      }
      catch (std::exception&) {
         printf("failed creating material");
      }
   }
}