/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/common.h
// Purpose:     common definitions for IHG specific code
// Created:     2005/03/31
// RCS-ID:      $Id: common.h,v 1.7 2006/08/20 09:38:12 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/ihg/common.h
    @brief  This file contains miscellaneous stuff which is used for IHG
            specific code.

    This file should be always included, directly or indirectly, by all IHG
    source code.
 */

#ifndef _ITO33_IHG_COMMON_H_
#define _ITO33_IHG_COMMON_H_

// do this to disable some warnings (beforestd.h does it and afterstd.h is
// needed to restore the warning level -- but it doesn't reenable these
// particular warnings)
#include "ito33/beforestd.h"
#include "ito33/afterstd.h"

#include "ito33/dlldecl.h"
#include "ito33/ihg/dlldecl.h"

#endif // _ITO33_IHG_COMMON_H_
