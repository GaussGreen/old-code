/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/showversion.h
// Purpose:     file containing version information for HG product
// Created:     2005/01/17
// RCS-ID:      $Id: showversion.h,v 1.1 2005/01/17 17:24:44 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_HG_SHOWVERSION_H_
#define _ITO33_HG_SHOWVERSION_H_

#include "ito33/string.h"
#include "ito33/hg/version.h"

namespace ito33
{

namespace hg
{

/**
   Gets the version for HG product.

   @return version for HG product
 */
inline std::string ShowVersion()
{
  return ITO33_HG_VERSION_DOT_STRING;
}


} // namespace hg

} // namespace ito33

#endif // #ifndef _ITO33_HG_SHOWVERSION_H_

