/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/common.h
// Purpose:     common definitions for HG specific code
// Created:     2005/04/06
// RCS-ID:      $Id: common.h,v 1.2 2006/08/19 23:46:21 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/hg/common.h
    @brief  This file contains miscellaneous stuff which is used for HG
            specific code.

    This file should be always included, directly or indirectly, by all HG
    source code.
 */

#ifndef _ITO33_HG_COMMON_H_
#define _ITO33_HG_COMMON_H_

// do this to disable some warnings (beforestd.h does it and afterstd.h is
// needed to restore the warning level -- but it doesn't reenable these
// particular warnings)
#include "ito33/beforestd.h"
#include "ito33/afterstd.h"

#include "ito33/dlldecl.h"
#include "ito33/hg/dlldecl.h"

#endif // _ITO33_HG_COMMON_H_
