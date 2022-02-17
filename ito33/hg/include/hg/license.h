/////////////////////////////////////////////////////////////////////////////
// Name:        hg/license.h
// Purpose:     license utils for HG project
// Created:     2005/01/13
// RCS-ID:      $Id: license.h,v 1.1.1.1 2005/01/17 15:52:17 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/license.h
   @brief declaration of class and functions for HG license issue
 */
#ifndef _HG_LICENSE_H_
#define _HG_LICENSE_H_

#include "ito33/beforestd.h"
#include "ito33/afterstd.h"

#include "ito33/license.h"

namespace ito33
{

namespace hg
{


License& GetHGLicense();

License& GetHGCBLicense();


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_LICENSE_H_ 
