/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/vector.h
// Purpose:     wrapper for the standard vector header
// Author:      Vadim Zeitlin
// Created:     16.04.03
// RCS-ID:      $Id: vector.h,v 1.3 2004/10/05 09:13:35 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/vector.h
    @brief  Just a wrapper for the standard vector header.

    Including the standard vector header results in too many warnings at
    maximum warning level with MS VC++ so we use this file to disable them.
 */

#ifndef _ITO33_VECTOR_H_
#define _ITO33_VECTOR_H_

#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"

#endif // _ITO33_VECTOR_H_


