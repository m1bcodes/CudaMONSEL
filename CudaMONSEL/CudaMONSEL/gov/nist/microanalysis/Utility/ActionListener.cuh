#ifndef _ACTION_LISTENER_CUH_
#define _ACTION_LISTENER_CUH_

namespace ActionListener
{
   class ActionListener
   {
   public:
      __host__ __device__ virtual void actionPerformed(const int) = 0;
   };
}

#endif