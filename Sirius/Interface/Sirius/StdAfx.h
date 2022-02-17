// stdafx.h : include file for standard system include files,
//			  or project specific include files that are used frequently,
//			  but are changed infrequently.
//
//////////////////////////////////////////////////////////////////////

#ifndef _STDAFX_H
#define _STDAFX_H

#pragma once

#define STRICT
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0501
#endif
#define _ATL_APARTMENT_THREADED

#include "warningsoff.h"

// ----------------------------------------------------------------------------
// STL
// ----------------------------------------------------------------------------
#include <math.h>
#include <string>
#include <vector>
#include <stack>
#include <map>
#include <hash_map>
#include <sstream>
#include <algorithm>
#include <limits>
#include <iostream>
#include <locale>
#include <fstream>
#include <memory>
#include <set>

// ----------------------------------------------------------------------------
// ATL + Windows Stuff
// ----------------------------------------------------------------------------
#include <atlbase_noassert.h>
#include <commctrl.h>
#ifdef ATLTRACE2
	#undef ATLTRACE2													// suppress trace output messages from ATL library
#endif
#define ATLTRACE2
#import "C:\Program Files\Common Files\System\ADO\msado15.dll"	\
		no_namespace											\
		rename("EOF", "adoEOF")									\
		rename("Parameters", "adoParameters")					\
		rename("Parameter", "adoParameter")
#include <sirius.h>

#include "siriuscommodule.h"
extern CSiriusComModule					_Module;						// do not change the name _Module
#ifndef max
	#define max( a, b ) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
	#define min( a, b ) ((a) < (b) ? (a) : (b))
#endif
#include <atlcom.h>
#undef max
#undef min
#include <atlwin.h>
#include <comdef.h>
#include <atlhost.h>
#include <sqltypes.h>
#include <sqlext.h>
#include <shlobj.h>						// for SQLRETURN etc.
#pragma comment(lib, "mpr.lib")			// for WNetGetUniversalName function etc.

// ----------------------------------------------------------------------------
// Additional Sirius includes
// ----------------------------------------------------------------------------
#include <siriusinterfacebase.h>
// For IN-PROCESS Bug whereby  creation of sirius application sub objects  hangs 
// due the use of DECLARE_CLASSFACTORY_SINGLETON which makes this creation occur
// when ole32.dll calls DllGetClassObject and then CreateInstance
#include <AltSingleton_SafeProc.h>

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
	//#include <GdaCommon.h>			// public parts of GdaCommon
	//#include <GdaInterface.h>		// public parts of GdaInterface
#endif

// ----------------------------------------------------------------------------
// Macros etc
// ----------------------------------------------------------------------------
#define check_publishing_enabled														{																										\
																							if (!_Module.GetEnablePublishing()) return CParameterMap::ReturnErrorR(IDS_PUBLISHING_NOT_ENABLED);	\
																						}

#define check_load_enabled																{																										\
																							if (!_Module.GetEnableLoad()) return CParameterMap::ReturnErrorR(IDS_LOAD_NOT_ENABLED);				\
																						}

																						
//	This handles the declaration and retrieval function for an object in the CApplication class.
#define	DECLARE_COLLECTION(CollectionName)												protected:																\
																							CComPtr<I##CollectionName>			m_sp##CollectionName;			\
																						public:																	\
																							STDMETHOD(get_##CollectionName)(I##CollectionName** pVal)			\
																							{																	\
																								if (!m_sp##CollectionName) return E_POINTER;					\
																								return m_sp##CollectionName.CopyTo(pVal);						\
																							}
//	These macros implement CApplication::InstallStandardObjects.
#define BEGIN_IMPLEMENT_COLLECTION()													HRESULT CSiriusApplication::InstallStandardObjects(void)						\
																						{																				\
																							HRESULT								hr;
#define IMPLEMENT_COLLECTION(LibraryName, CollectionName, ObjectName, ProductLike)			if (hr = m_ObjectManager.InstallStandardObject(LibraryName, CLSID_##ObjectName, IID_I##ObjectName, ProductLike, m_sp##CollectionName, CLSID_##CollectionName)) return hr;
#define END_IMPLEMENT_COLLECTION()															return S_OK;																\
																						}
//{{AFX_INSERT_LOCATION}}

#endif
