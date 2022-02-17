/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/erroroffsets.h
// Purpose:     Manages the error code range for all ITO modules.
// Author:      ZHANG Yunzhi
// Created:     16.12.02
// RCS-ID:      $Id: erroroffsets.h,v 1.3 2006/06/15 18:57:04 zhang Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/erroroffsets.h
    @brief This file manages the error code ranges for all ITO modules.

    That means, this file is the only place to define and so to give a quick
    view of all ITO33 project modules and their error code ranges.

    Note that only "big" modules can reserve code range here. Each module
    itself should manage the code ranges for its sub-modules, if any.
 */
#ifndef _ITO33_ERROROFFSETS_H_
#define _ITO33_ERROROFFSETS_H_

/**
    Ito33's error codes start from 10000.
    In particular, basic error codes range is from 10000 to 19999
 */
#define ITO33_BASE_ERROR_START 10000

/// finance error codes range is from 20000 to 29999
#define ITO33_FINANCE_ERROR_START 20000

/// ihg error codes range is from 30000 to 39999
#define ITO33_IHG_ERROR_START 30000

/// ihg error codes range is from 40000 to 49999
#define ITO33_HG_ERROR_START 40000

/// opsession common (library etc.) error codes range is from 50000 to 59999
#define ITO33_OPSESSION_ERROR_START 50000

/// opsession inblue error codes range is from 60000 to 69999
#define ITO33_INBLUE_ERROR_START 60000

/// opssesion opsgui error codes range is from 70000 to 79999
#define ITO33_OPSGUI_ERROR_START 70000

/// datastore codes start from 80000.
#define ITO33_DATASTORE_ERROR_START 80000

#endif /* _ITO33_ERROROFFSETS_H_ */
