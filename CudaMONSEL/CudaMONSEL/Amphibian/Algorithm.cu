#include "Amphibian\Algorithm.cuh"

namespace Algorithm
{
   __host__ __device__ int binarySearch(const double a[], int l, int r, const double t)
   {
      while (l <= r) {
         const int mid = l + (r - l) / 2;
         if (a[mid] == t) return mid;
         if (a[mid] < t) l = mid + 1;
         else r = mid - 1;
      }
      return l;
   }

   __host__ __device__ int binarySearch(const float a[], int l, int r, const float t)
   {
      while (l <= r) {
         const int mid = l + (r - l) / 2;
         if (a[mid] == t) return mid;
         if (a[mid] < t) l = mid + 1;
         else r = mid - 1;
      }
      return l;
   }

   __host__ __device__ static int partition(double a[], int l, int r)
   {
      int j = l - 1; // points to the last element less than pivot
      for (int i = l; i < r; ++i) { // going from l to r
         if (a[i] < a[r]) { // when sees an element less than pivot
            ++j; // points to the first element greater than or equal to pivot
            // swap i and j, so j points to last element less than pivot
            double t = a[i];
            a[i] = a[j];
            a[j] = t;
         }
      }
      ++j; // points to first element greater than or equal to pivot
      // swap
      double t = a[r];
      a[r] = a[j];
      a[j] = t;

      return j;
   }

   __host__ __device__ void quicksort(double a[], int l, int r)
   {
      if (l >= r) return;
      int pv = partition(a, l, r);

      quicksort(a, l, pv - 1);
      quicksort(a, pv + 1, r);
   }
}