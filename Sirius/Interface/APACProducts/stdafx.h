// stdafx.h : include file for standard system include files,
//			  or project specific include files that are used frequently,
//			  but are changed infrequently
//
//////////////////////////////////////////////////////////////////////

#ifndef _STDAFX_H
#define _STDAFX_H

#pragma once

#define STRICT
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0400
#endif
#define _ATL_APARTMENT_THREADED

#include "warningsoff.h"

// ----------------------------------------------------------------------------
// STL
// ----------------------------------------------------------------------------
#include <math.h>
#include <vector>
#include <map>
#include <hash_map>
#include <sstream>
#include <algorithm>
#include <fstream.h>
#include <memory>
#include <set>

// ----------------------------------------------------------------------------
// ATL + Sirius
// ----------------------------------------------------------------------------
#include <atlbase_noassert.h>
#include <sirius.h>
#include "apacproductscommodule.h"
extern CAPACProductsComModule		_Module;
#include <atlcom.h>
#include <comdef.h>

// ----------------------------------------------------------------------------
// Additional Sirius includes
// ----------------------------------------------------------------------------
#include <siriusinterfacebase.h>

// ----------------------------------------------------------------------------
// GDA
// ----------------------------------------------------------------------------
#ifdef _DEBUG
	#define SIRIUS_DEBUG_TEMP_COPY _DEBUG
	#undef _DEBUG
	#include <GdaCommon.h>			// public parts of GdaCommon
	#include <GdaInterface.h>		// public parts of GdaInterface
	#define _DEBUG SIRIUS_DEBUG_TEMP_COPY 
	#undef SIRIUS_DEBUG_TEMP_COPY
#else
	#include <GdaCommon.h>			// public parts of GdaCommon
	#include <GdaInterface.h>		// public parts of GdaInterface
#endif

// ----------------------------------------------------------------------------
// ModelLib includes
// ----------------------------------------------------------------------------
#include <mleqnormal.h>

//{{AFX_INSERT_LOCATION}}

#endif
