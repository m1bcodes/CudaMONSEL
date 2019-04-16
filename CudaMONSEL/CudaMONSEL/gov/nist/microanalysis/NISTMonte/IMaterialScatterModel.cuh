//#ifndef _I_MATERIAL_SCATTER_MODEL_CUH_
//#define _I_MATERIAL_SCATTER_MODEL_CUH_
//
//class RegionBase;
//class Electron::Electron;
//class Material::Material;
//
//class IMaterialScatterModel
//{
//public:
//   virtual bool operator==(const IMaterialScatterModel& other) { return this == &other; };
//
//   virtual Material::Material& getMaterial() = 0;
//   virtual double getMinEforTracking() = 0;
//   virtual void setMinEforTracking(double minEforTracking) = 0;
//   virtual double randomMeanPathLength(Electron::Electron& pe) = 0;
//   virtual Electron::Electron& scatter(Electron::Electron& pe) = 0;
//   virtual Electron::Electron& barrierScatter(Electron::Electron& pe, RegionBase& nextRegion) = 0;
//   virtual double calculateEnergyLoss(double len, Electron::Electron& pe) = 0;
//};
//
//#endif