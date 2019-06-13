#ifndef _MATERIAL_FACTORY_CUH_
#define _MATERIAL_FACTORY_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace MaterialFactory
{
   extern const char* K3189;
   extern const char* RynasAlTiAlloy;
   extern const char* Mylar;
   extern const char* VanadiumPentoxide;
   extern const char* SiliconDioxide;
   extern const char* Ice;
   extern const char* PerfectVacuum;
   extern const char* CaCO3;
   extern const char* Al2O3;
   extern const char* SS316;
   extern const char* UraniumOxide;
   extern const char* K227;
   extern const char* K309;
   extern const char* K411;
   extern const char* K412;
   extern const char* K961;
   extern const char* K1080;
   extern const char* K2450;
   extern const char* K2451;
   extern const char* K2466;
   extern const char* K2469;
   extern const char* K2472;
   extern const char* K2496;
   extern const char* ParaleneC;
   extern const char* MagnesiumOxide;
   extern const char* Chloroapatite;
   extern const char* CalciumFluoride;
   extern const char* GalliumPhosphate;
   extern const char* Nothing;

   extern bool canCreate(const ElementT& el);
   // extern ElementMap parseCompound(StringT str);
   //extern CompositionT createCompound1(StringT str);
   extern CompositionT createCompound(StringT str);
   extern CompositionT createMaterial(StringT name);

   void init();
}

#endif
