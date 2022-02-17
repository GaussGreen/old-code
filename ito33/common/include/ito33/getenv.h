/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/getenv.h
// Purpose:     defines a wrapper of getenv()
// Created:     2005/02/01
// RCS-ID:      $Id: getenv.h,v 1.2 2006/05/31 15:04:29 willy Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/getenv.h
    @brief  defines wrapper of getenv() because getenv() was declared
            deprecated in VC8.
 */

#ifndef _ITO33_GETENV_H_
#define _ITO33_GETENV_H_

#include <cstdlib>

namespace ito33
{

/**
    Gets a value from a environment variable.

    @param envName the name of the environement variable.
    @return the value of the environement variable.
 */
inline std::string GetEnv(const std::string& envName)
{
  #if ( defined(_MSC_VER) && _MSC_VER >= 1400 )
  {
    size_t requiredSize;
    getenv_s( &requiredSize, NULL, 0, envName.c_str());

    std::string ret(requiredSize, ' ');
    if ( requiredSize )
      getenv_s( &requiredSize, &ret[0], requiredSize, envName.c_str());

    return ret;
  }
  #else
  {
    std::string ret;
	const char *instDir = getenv(envName.c_str());
	if (  instDir )
		ret = instDir;

    return ret;
  }
  #endif
}

} // namespace ito33

#endif // #ifndef _ITO33_GETENV_H_
