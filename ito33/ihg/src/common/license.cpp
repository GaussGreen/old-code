/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/license.cpp
// Purpose:     Implementation of license class and functions
// Author:      ZHANG Yunzhi
// Created:     2004/06/30
// RCS-ID:      $Id: license.cpp,v 1.9 2006/06/12 17:03:49 zhang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/ihg/version.h"
#include "ihg/license.h"
#include "ito33/ihg/error.h"

extern const ito33::ihg::Error
  ITO33_IHG_LICENSE,
  ITO33_IHG_LICENSE_CB; 

namespace ito33
{

namespace ihg
{

License& GetIHGLicense()
{
  static License 
            licIHG("ihg_base",
                   ITO33_MAKE_SERIES_DOT_STRING(ITO33_IHG_VERSION_MAJOR, 
                                                ITO33_IHG_VERSION_MINOR),
                   ITO33_IHG_LICENSE,
                   ITO33_IHG_LICENSE.GetMessage()
                  );
  return licIHG;
}

License& GetIHGCBLicense()
{
  static License 
            licIHG("ihg_cb",
                   ITO33_MAKE_SERIES_DOT_STRING(ITO33_IHG_VERSION_MAJOR, 
                                                ITO33_IHG_VERSION_MINOR),
                   ITO33_IHG_LICENSE_CB,
                   ITO33_IHG_LICENSE_CB.GetMessage()
                  );
  return licIHG;
}



} // namespace ihg

} // namespace ito33

