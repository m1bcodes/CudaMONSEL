#ifndef _MATERIAL_FACTORY_CUH_
#define _MATERIAL_FACTORY_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace MaterialFactory
{
   extern const StringT K3189;
   extern const StringT RynasAlTiAlloy;
   extern const StringT Mylar;
   extern const StringT VanadiumPentoxide;
   extern const StringT SiliconDioxide;
   extern const StringT Ice;
   extern const StringT PerfectVacuum;
   extern const StringT CaCO3;
   extern const StringT Al2O3;
   extern const StringT SS316;
   extern const StringT UraniumOxide;
   extern const StringT K227;
   extern const StringT K309;
   extern const StringT K411;
   extern const StringT K412;
   extern const StringT K961;
   extern const StringT K1080;
   extern const StringT K2450;
   extern const StringT K2451;
   extern const StringT K2466;
   extern const StringT K2469;
   extern const StringT K2472;
   extern const StringT K2496;
   extern const StringT ParaleneC;
   extern const StringT MagnesiumOxide;
   extern const StringT Chloroapatite;
   extern const StringT CalciumFluoride;
   extern const StringT GalliumPhosphate;
   extern const StringT Nothing;

   extern bool canCreate(const ElementT& el);
   // extern ElementMap parseCompound(StringT str);
   //extern CompositionT createCompound1(StringT str);
   extern CompositionT createCompound(StringT str);
   extern CompositionT createMaterial(StringT name);

   void init();
}

#endif
