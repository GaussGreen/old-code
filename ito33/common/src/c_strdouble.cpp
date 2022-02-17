/////////////////////////////////////////////////////////////////////////////
// Name:        utils/c_strdouble.cpp
// Purpose:     ito33::String::From/ToCDouble functions
// Author:      Vaclav Slavik
// Created:     2005-05-30
// RCS-ID:      $Id: c_strdouble.cpp,v 1.3 2006/08/02 13:29:27 cosmin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <stdexcept>
#include <algorithm>

#ifdef _MSC_VER
  #include <float.h>
#else
  #include <math.h>
#endif

#ifdef _MSC_VER
  #define isfinite(x)   _finite(x)
#endif

namespace ito33
{

namespace String
{

// functions & variables from c_dtoa.cpp
extern double c_strtod(const char *s00, char **se);
extern char *c_dtoa(double d, int mode, int ndigits,
                    int *decpt, int *sign, char **rve);
extern void c_freedtoa(char *s);
extern int c_dtoa_NO_DECIMAL_POINT;

std::string FromCDouble(double x, int precision)
{
  int decimalPoint, sign;

  char *str;

  if ( precision == -1 )
  {
    str = c_dtoa(x, 0, 0, &decimalPoint, &sign, NULL);
  }
  else
  {
    if ( isfinite(x) )
    {
      // determine # of digits after decimal point, by counting how many
      // of them we need to express trunc(x):
      int digs = 0;
      for ( double x2 = x; x2 >= 1.0; x2 /= 10 )
        digs++;
      precision = std::max(0, precision - digs);
    }
    else
    {
      precision = 0;
    }

    // and use it for formatting:
    str = c_dtoa(x, 3, precision, &decimalPoint, &sign, NULL);
  }

  std::string result(str);
  c_freedtoa(str);

  if ( decimalPoint == c_dtoa_NO_DECIMAL_POINT )
  {
    // NaN or inifinity don't have decimal point in str representation
    return result;
  }

  int len = (int)result.length();

  if ( decimalPoint != len )
  {
    if ( decimalPoint < 0 )
    {
      result.insert((std::string::size_type)0, -decimalPoint, '0');
      result.insert(0, "0.");
    }
    else if ( decimalPoint == 0 )
    {
      result.insert(0, "0.");
    }
    else if ( decimalPoint == len )
    {
      // do nothing
    }
    else if ( decimalPoint > len )
    {
      result.append(decimalPoint - len, '0');
    }
    else
    {
      result.insert(decimalPoint, ".");
    }
  }

  if ( sign )
    result.insert(0, "-");

  return result;
}

bool ToCDouble(const std::string& s, double *val)
{
  char *end;
  const char *start = s.c_str();
  *val = c_strtod(start, &end);

  return !*end && (end != start);
}

} // namespace String

} // namespace ito33
