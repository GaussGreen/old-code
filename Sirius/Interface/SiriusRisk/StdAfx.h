// stdafx.h : Include file for standard system include files,
//			  or project specific include files that are used frequently,
//			  but are changed infrequently.
//
// Author:	  David Cuin.
//
// Note:      Don't mess about with this file!
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
// MFC Classes
// ----------------------------------------------------------------------------
#include <afxwin.h>
#include <afxdisp.h>
#include <afxdao.h>
#ifdef _MSC_VER
	#pragma warning(disable:4275) // non dll-interface class 'exception' used as base for dll-interface class
#endif


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
#include "siriusriskcommodule.h"
extern CSiriusRiskComModule		_Module;
#include <atlcom.h>
#include <comdef.h>

// ----------------------------------------------------------------------------
// Additional Sirius includes
// ----------------------------------------------------------------------------
#include <siriusinterfacebase.h>

// ----------------------------------------------------------------------------
// GDA
// ----------------------------------------------------------------------------
// #include <GdaCommon.h>			// public parts of GdaCommon
// #include <GdaInterface.h>		// public parts of GdaInterface

// ----------------------------------------------------------------------------
// ModelLib includes
// ----------------------------------------------------------------------------
// #include <mleqnormal.h>

//{{AFX_INSERT_LOCATION}}

#endif
