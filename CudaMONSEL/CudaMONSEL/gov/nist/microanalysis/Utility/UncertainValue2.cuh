#ifndef _UNCERTAIN_VALUE_2_CUH_
#define _UNCERTAIN_VALUE_2_CUH_

#include "..\..\..\..\Amphibian\LinkedList.cuh"
#include "..\..\..\..\Amphibian\String.cuh"

class UncertainValue2
{
public:
   static const char DEFAULT[];

   static const UncertainValue2 ONE;
   static const UncertainValue2 ZERO;
   //static const UncertainValue2 NaN;
   //static const UncertainValue2 POSITIVE_INFINITY;
   //static const UncertainValue2 NEGATIVE_INFINITY;

   UncertainValue2(double v, char source[], double dv);
   UncertainValue2(double v);
   UncertainValue2(double v, double dv);
   UncertainValue2(double v, Node<String, double>* sigmas);

   void assignComponent(String name, double sigma);
   double fractionalUncertainty();
   double doubleValue();
   bool isUncertain();
   double variance();
   double uncertainty();
   bool equals(UncertainValue2 * const uv);

   bool isNaN();

private:
   static const int MAX_LEN = 11;

   static int sDefIndex; // transient
   static const long serialVersionUID;

   const double mValue;
   Node<String, double>* mSigmas;

   bool mNotANumber;
   bool mPosInfinity, mNegInfinity;
};

#endif