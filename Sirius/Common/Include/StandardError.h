///////////////////////////////////////////////////////////////////////////////
// Filename: StandardError.h 
// Author: Adam Gladstone 
// Date: Saturday, July 7, 2003 
// 
// Description
//    Error handling for the server.
//
///////////////////////////////////////////////////////////////////////////////

#if !defined(CLSGEN_CSTANDARDERROR_H__63902D4E_5E09_11D7_BC27_00805F9B19DC__INCLUDED_)
#define CLSGEN_CSTANDARDERROR_H__63902D4E_5E09_11D7_BC27_00805F9B19DC__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000


/////////////////////////////////////////////////////////////////////////////////
// Definitions

#define STDTRY { \
			try 


#define STDCATCH() \
    catch (const _com_error &e) { \
		CStdString strError; \
		if(e.ErrorInfo() != NULL) \
		{ \
			strError = e.Description(); \
		} \
		ATLTRACE(_T("%s(%d) : %s\n"), __FILE__, __LINE__, strError.c_str()); \
		return AtlReportError(GetObjectCLSID(), (LPTSTR)(LPCTSTR)strError, guid, e.Error()); \
    } \
	catch (const std::bad_alloc&) { \
		CStdString strError; \
		strError.Format(_T("An allocation error has occurred.")); \
		ATLTRACE(_T("%s(%d) : %s\n"), __FILE__, __LINE__, strError.c_str()); \
		return AtlReportError(GUID_NULL, (LPTSTR)(LPCTSTR)strError, GUID_NULL, E_OUTOFMEMORY); \
	} \
	catch (const std::exception& rcException) { \
		CStdString strError; \
		strError.Format(_T("An unexpected STD error has occurred. %s."), rcException.what()); \
		ATLTRACE(_T("%s(%d) : %s\n"), __FILE__, __LINE__, strError.c_str()); \
		return AtlReportError(GUID_NULL, (LPTSTR)(LPCTSTR)strError, GUID_NULL, E_UNEXPECTED); \
	} \
	catch (...)	{ \
		CStdString strError; \
		strError.Format(_T("An unexpected error has occurred.")); \
		ATLTRACE(_T("%s(%d) : %s\n"), __FILE__, __LINE__, strError.c_str()); \
		return AtlReportError(GUID_NULL, (LPTSTR)(LPCTSTR)strError, GUID_NULL, E_UNEXPECTED); \
	} }




#define COMTRY_EX { \
			try 




#define COMCATCH_EX(guid) \
    catch (const _com_error &e) { \
		CStdString strError; \
		if(e.ErrorInfo() != NULL) \
		{ \
			strError = e.Description(); \
		} \
		ATLTRACE(_T("%s(%d) : %s\n"), __FILE__, __LINE__, strError.c_str()); \
		return AtlReportError(GetObjectCLSID(), (LPTSTR)(LPCTSTR)strError, guid, e.Error()); \
    } \
	catch (const std::bad_alloc&) { \
		CStdString strError; \
		strError.Format(_T("An allocation error has occurred.")); \
		ATLTRACE(_T("%s(%d) : %s\n"), __FILE__, __LINE__, strError.c_str()); \
		return AtlReportError(GetObjectCLSID(), (LPTSTR)(LPCTSTR)strError, guid, E_OUTOFMEMORY); \
	} \
	catch (const std::exception& rcException) { \
		CStdString strError; \
		strError.Format(_T("An unexpected STD error has occurred.\n%s."), rcException.what()); \
		ATLTRACE(_T("%s(%d) : %s\n"), __FILE__, __LINE__, strError.c_str()); \
		return AtlReportError(GetObjectCLSID(), (LPTSTR)(LPCTSTR)strError, guid, E_UNEXPECTED); \
	} \
	catch (...)	{ \
		CStdString strError; \
		strError.Format(_T("An unexpected error has occurred.")); \
		ATLTRACE(_T("%s(%d) : %s\n"), __FILE__, __LINE__, strError.c_str()); \
		return AtlReportError(GetObjectCLSID(), (LPTSTR)(LPCTSTR)strError, guid, E_UNEXPECTED); \
	} }



/////////////////////////////////////////////////////////////////////////////////
// Includes

/////////////////////////////////////////////////////////////////////////////////
// Forward Declarations

/////////////////////////////////////////////////////////////////////////////////
// Class Definition

#endif		//CLSGEN_CSTANDARDERROR_H__63902D4E_5E09_11D7_BC27_00805F9B19DC__INCLUDED_
