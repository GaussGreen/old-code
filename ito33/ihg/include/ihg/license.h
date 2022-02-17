/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/license.h
// Purpose:     license utils for ihg project
// Author:      ZHANG Yunzhi
// RCS-ID:      $Id: license.h,v 1.6 2005/03/30 12:32:22 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ihg/license.h
  @brief declaration of class and functions for ihg license issue

 */
#ifndef _IHG_LICENSE_H_
#define _IHG_LICENSE_H_

#include "ito33/beforestd.h"
#include "ito33/afterstd.h"

#include "ito33/license.h"

namespace ito33
{

namespace ihg
{

License& GetIHGLicense();

License& GetIHGCBLicense();

}
}

#endif // _IHG_LICENSE_H_ 
