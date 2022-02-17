/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/list.h
// Purpose:     wrapper for the standard list header
// Author:      Vadim Zeitlin
// Created:     13.02.03
// RCS-ID:      $Id: list.h,v 1.4 2004/10/05 09:13:35 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/list.h
    @brief  Wrapper for the standard list header.

    Including the standard list header results in too many warnings at maximum
    warning level with MS VC++ so we use this file to disable them.
 */

#ifndef _ITO33_LIST_H_
#define _ITO33_LIST_H_

#include "ito33/beforestd.h"
#include <list>
#include "ito33/afterstd.h"

#endif // _ITO33_LIST_H_


