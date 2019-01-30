#include "Math.cuh"

__device__ __host__ double Math::sqrt(double number, double error)
{
   double s = number;

   while ((s - number / s) > error) {
      s = (s + number / s) / 2;
   }
   return s;
}

__device__ __host__ double Math::abs(double number)
{
   return number > 0 ? number : -number;
}
