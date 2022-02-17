/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/license.cpp
// Purpose:     Implementation of HG license management
// Created:     2005/01/13
// RCS-ID:      $Id: license.cpp,v 1.2 2006/06/12 17:03:49 zhang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/hg/version.h"

#include "ito33/hg/error.h"
#include "hg/license.h"

extern const ito33::hg::Error ITO33_HG_LICENSE; 

namespace ito33
{

namespace hg
{

License& GetHGLicense()
{
  static License 
    licHG("hg_base",
          ITO33_MAKE_SERIES_DOT_STRING(ITO33_HG_VERSION_MAJOR, 
                                       ITO33_HG_VERSION_MINOR),
          ITO33_HG_LICENSE,
          ITO33_HG_LICENSE.GetMessage());

  return licHG;
}

License& GetHGCBLicense()
{
  static License 
    licHG("hg_cb",
          ITO33_MAKE_SERIES_DOT_STRING(ITO33_HG_VERSION_MAJOR, 
                                       ITO33_HG_VERSION_MINOR),
          ITO33_HG_LICENSE,
          ITO33_HG_LICENSE.GetMessage());
  
  return licHG;
}


} // namespace hg

} // namespace ito33
