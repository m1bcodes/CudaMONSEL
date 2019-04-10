//#include "MonteCarloSS.cu"
//
//const float MonteCarloSS::ChamberRadius = 0.1f;
//const float MonteCarloSS::SMALL_DISP = 1.0e-15f;
//
//class MonteCarloSS::RegionBase
//{
//protected:
//   int id;
//public:
//   __host__ __device__ RegionBase(int k) { this->id = id; }
//   //__host__ __device__ ~RegionBase() { };
//
//   __host__ __device__ int GetId() { return id; }
//};
//
//class MonteCarloSS::ElectronGun
//{
//public:
//   virtual void setBeamEnergy(double beamEnergy) = 0;
//   virtual double getBeamEnergy() = 0;
//   virtual void setCenter(double center[]) = 0;
//   virtual double* getCenter() = 0;
//   virtual Electron createElectron() = 0;
//};
//
//class MonteCarloSS::Shape
//{
//   virtual bool contains(double pos[]) = 0;
//   virtual double getFirstIntersection(double pos0[], double pos1[]) = 0;
//};
//
