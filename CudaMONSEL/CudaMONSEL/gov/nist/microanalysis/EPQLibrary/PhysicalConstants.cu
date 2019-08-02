#include "PhysicalConstants.cuh"

namespace PhysicalConstants
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const float AvagadroNumber = 6.0221415e23f; // dimensionless
   __constant__ const float SpeedOfLight = 299792458.f; // m/s

   __constant__ const float ElectronCharge = 1.60217653e-19f; // c
   __constant__ const float PlanckConstant = 6.6260693e-34f; // j s
   __constant__ const float PlanckReduced = 1.05457266e-34f; // j s
   __constant__ const float ElectronMass = 9.1093826e-31f; // kg

   __constant__ const float ElectronRestMass = 8.1871047e-14f; // j
   __constant__ const float ProtonMass = 1.66053886e-27f; // kg
   __constant__ const float NeutronMass = 1.67492728e-27f; // kg
   __constant__ const float UnifiedAtomicMass = 1.66053886e-27f; // kg
   __constant__ const float PermittivityOfFreeSpace = 8.8541878176203898505365630317108e-12f; // f/m
   __constant__ const float PermeabilityOfFreeSpace = 12.566370614359172953850573533118e-7f; // n/(a^2)
   __constant__ const float BoltzmannConstant = 1.3806505e-23f; // j/k

   __constant__ const float GravitationConstant = 6.6742e-11f; // kg m^2
   __constant__ const float PlanckLength = 1.61624e-35f; // m
   __constant__ const float PlanckMass = 2.17645e-8f; // kg
   __constant__ const float PlanckTemperature = 1.41679e32f; // dimensionless
   __constant__ const double PlanckTime = 5.39121e-44f; // s

   __constant__ const float RydbergEnergy = 2.17987209e-18f; // joules

   __constant__ const float BohrRadius = 5.291772083e-11f; // meter

   __constant__ const float FineStructure = 7.297352568e-3f;

   __constant__ const float ClassicalElectronRadius = 2.817940325e-15f;

   __constant__ const float IcePoint = 273.15f; // kelvin
   __constant__ const float StandardAtmosphere = 101325.0f; // pascal

   __constant__ const float PI = 3.14159265358979323846f;
#else
   const float AvagadroNumber = 6.0221415e23f; // dimensionless
   const float SpeedOfLight = 299792458.f; // m/s

   const float ElectronCharge = 1.60217653e-19f; // c
   const float PlanckConstant = 6.6260693e-34f; // j s
   const float PlanckReduced = 1.05457266e-34f; // j s
   const float ElectronMass = 9.1093826e-31f; // kg

   const float ElectronRestMass = 8.1871047e-14f; // j
   const float ProtonMass = 1.66053886e-27f; // kg
   const float NeutronMass = 1.67492728e-27f; // kg
   const float UnifiedAtomicMass = 1.66053886e-27f; // kg
   const float PermittivityOfFreeSpace = 8.8541878176203898505365630317108e-12f; // f/m
   const float PermeabilityOfFreeSpace = 12.566370614359172953850573533118e-7f; // n/(a^2)
   const float BoltzmannConstant = 1.3806505e-23f; // j/k

   const float GravitationConstant = 6.6742e-11f; // kg m^2
   const float PlanckLength = 1.61624e-35f; // m
   const float PlanckMass = 2.17645e-8f; // kg
   const float PlanckTemperature = 1.41679e32f; // dimensionless
   const double PlanckTime = 5.39121e-44f; // s

   const float RydbergEnergy = 2.17987209e-18f; // joules

   const float BohrRadius = 5.291772083e-11f; // meter

   const float FineStructure = 7.297352568e-3f;

   const float ClassicalElectronRadius = 2.817940325e-15f;

   const float IcePoint = 273.15f; // kelvin
   const float StandardAtmosphere = 101325.0f; // pascal

   const float PI = 3.14159265358979323846f;
#endif
}