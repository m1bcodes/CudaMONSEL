#ifndef _LINES_ON_LAYERS_0_CUH_
#define _LINES_ON_LAYERS_0_CUH_

#include "gov\nist\nanoscalemetrology\JMONSEL\NShapes.cuh"

#include <cuda_runtime.h>

namespace LinesOnLayers
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   extern __device__ unsigned int ysize;
   extern __device__ unsigned int xsize;

   extern __device__ unsigned int nTrajectories;
   extern __device__ float beamsizenm;
   extern __device__ float beamsize;
   extern __device__ float beamEeV;
   extern __device__ float beamE;
   extern __device__ NShapes::LineParams** lineParams;
   extern __device__ unsigned int nlines;
   extern __device__ unsigned int linemat;
#else
   extern unsigned int ysize;
   extern unsigned int xsize;

   extern unsigned int nTrajectories;
   extern float beamsizenm;
   extern float beamsize;
   extern float beamEeV;
   extern float beamE;
   extern NShapes::LineParams** lineParams;
   extern unsigned int nlines;
   extern unsigned int linemat;
#endif
   extern void loadNUTable();
   extern void transferNUTableToCuda();

   extern __global__ void initCuda();
   extern __global__ void runCuda(float*, float*, float*);

   extern __host__ __device__ void runSinglePixel(const unsigned int r, const unsigned int c, float*, float*, float*);
   extern __host__ __device__ void initImgRange();

   extern void runSinglePixelThread(int id, const unsigned int r, const unsigned int c, float*, float*, float*);

   //extern void setSimParamsFromCSVFile(const char* fname);
   extern __host__ __device__ void setSimParams();
   extern __host__ __device__ void destroySimParams();

   extern void lineProjection(const unsigned int n, char* gt);

   extern void writeSerializedParams(const char res[]);
   extern void readSerializedParams(const char fname[]);
}

#endif