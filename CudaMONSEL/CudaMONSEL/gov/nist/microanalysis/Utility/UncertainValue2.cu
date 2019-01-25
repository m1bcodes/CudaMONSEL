#ifndef UNCERTAIN_VALUE_2_H
#define UNCERTAIN_VALUE_2_H

class UncertainValue2
{
public:
   static const int MAX_LEN = 15;
   static const char DEFAULT[];

   static const UncertainValue2 ONE;
   static const UncertainValue2 ZERO;
   static const UncertainValue2 NaN;
   static const UncertainValue2 POSITIVE_INFINITY;
   static const UncertainValue2 NEGATIVE_INFINITY;

   UncertainValue2(double v, char source[], double dv);
   UncertainValue2(double v);
   UncertainValue2(double v, double dv);
   UncertainValue2(double v, char* sigmaName[], float sigmaVal);

private:
   struct Sigma
   {
      const char Name[MAX_LEN];
      const float Value;
      Sigma* Next;
   };

   static int sDefIndex; // transient
   static const long serialVersionUID;

   const double mValue;
   const char* mSigmaNames[];
   const float mSigmaVals;
};

#endif