
#ifndef _ITO33_TESTS_ERROR_QUALITY_H_
#define _ITO33_TESTS_ERROR_QUALITY_H_

#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/common.h"

namespace ito33
{


enum ErrorQuality
{
  ErrorQuality_pass,
  ErrorQuality_warn,
  ErrorQuality_fail,
  ErrorQuality_new,

  ErrorQuality_Max
};

class DifferenceQuality
{
public:
  DifferenceQuality() :
      lowLevel(ErrorQuality_pass), strictLevel(ErrorQuality_pass)
  {
  }

  ErrorQuality lowLevel;
  ErrorQuality strictLevel;
};

inline DifferenceQuality ComparePrice(double dPrice1, double dPrice2)
{
  DifferenceQuality result;
  
  double dDiff = dPrice1 - dPrice2;
  if( fabs(dPrice1) > 100 )
    dDiff = (dDiff / dPrice1) * 100;
  dDiff = fabs(dDiff);
  
  result.lowLevel = (dDiff < 0.05) ? ErrorQuality_pass : ErrorQuality_fail;
  result.strictLevel = (dDiff < 0.02) ? ErrorQuality_pass : ErrorQuality_fail;

  return result;
}


}

#endif // #define _ITO33_TESTS_ERROR_QUALITY_H_
