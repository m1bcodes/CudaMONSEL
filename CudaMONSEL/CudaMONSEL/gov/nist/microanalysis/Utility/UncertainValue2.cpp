//#include "UncertainValue2.cu"
//
//#include <stdlib.h>
//
//const char UncertainValue2::DEFAULT[] = "Default";
//
//int UncertainValue2::sDefIndex = 0;
//const long UncertainValue2::serialVersionUID = 119495064970078787L;
//
//const UncertainValue2 UncertainValue2::ONE(1.0);
//const UncertainValue2 UncertainValue2::ZERO(0.0);
////const UncertainValue2 UncertainValue2::NaN(Double.NaN);
////const UncertainValue2 UncertainValue2::POSITIVE_INFINITY(Double.POSITIVE_INFINITY);
////const UncertainValue2 UncertainValue2::NEGATIVE_INFINITY(Double.NEGATIVE_INFINITY);
//
//UncertainValue2::UncertainValue2(double v, double dv) : mValue(v)
//{
//   char tmpName[MAX_LEN];
//   itoa(++sDefIndex, tmpName, 10);
//   UncertainValue2::UncertainValue2(v, tmpName, dv);
//}
//
//UncertainValue2::UncertainValue2(double v) : mValue(v)
//{
//   UncertainValue2::UncertainValue2(v, 0.0);
//}
//
//UncertainValue2::UncertainValue2(double v, char source[], double dv) : mValue(v)
//{
//   assignComponent(source, dv);
//}
//
//UncertainValue2::UncertainValue2(double v, char* sigmaName[], float sigmaVal) : mValue(v)
//{
//   for (Map.Entry<String, Double> me : sigmas.entrySet())
//      assignComponent(me.getKey(), me.getValue());
//}
