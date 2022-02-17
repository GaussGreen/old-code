//	comdispatchdriverex.h: Interface for the CComDispatchDriverEx class
//
//////////////////////////////////////////////////////////////////////

#ifndef _COMDISPATCHDRIVER_H
#define _COMDISPATCHDRIVER_H

#pragma once

#include "atlcom.h"

class CComDispatchDriverEx : public CComDispatchDriver
{
protected:
	static HRESULT GetPropertyWithParameters(OLECHAR FAR* szMember, const CComDispatchDriverEx& dd, VARIANT* pparams, unsigned int nParams, VARIANT* result)
	{
		DISPID								dispid;
		HRESULT								hr;
		DISPPARAMS							dispparams;

		if (hr = dd.p->GetIDsOfNames(IID_NULL, &szMember, 1, LOCALE_SYSTEM_DEFAULT, &dispid)) return hr;
		dispparams.rgdispidNamedArgs = NULL;
		dispparams.cArgs = nParams;
		dispparams.cNamedArgs = 0;
		dispparams.rgvarg = pparams;
		hr = dd.p->Invoke(dispid, IID_NULL, LOCALE_SYSTEM_DEFAULT, DISPATCH_PROPERTYGET, &dispparams, result, NULL, NULL);
		// ToDo - do something with pExcepInfo and puArgErr
		return hr;
	}
	static HRESULT PutPropertyWithParameters(OLECHAR FAR* szMember, const CComDispatchDriverEx& dd, VARIANT* pparams, unsigned int nParams)
	{
		DISPID								dispid;
		HRESULT								hr;
		DISPPARAMS							dispparams;
		DISPID								dispidPut = DISPID_PROPERTYPUT;

		if (hr = dd.p->GetIDsOfNames(IID_NULL, &szMember, 1, LOCALE_SYSTEM_DEFAULT, &dispid)) return hr;		
		dispparams.rgdispidNamedArgs = &dispidPut;		
		dispparams.cArgs = nParams;
		dispparams.cNamedArgs = 1;
		dispparams.rgvarg = pparams;
		hr = dd.p->Invoke(dispid, IID_NULL, LOCALE_SYSTEM_DEFAULT, DISPATCH_PROPERTYPUT, &dispparams, NULL, NULL, NULL);
		return hr;
	}

public:	
	CComDispatchDriverEx(){}
	CComDispatchDriverEx(IDispatch* lp) : CComDispatchDriver(lp){}
	CComDispatchDriverEx(IUnknown* lp) : CComDispatchDriver(lp){}
	void Release()
	{
		if (p) p->Release(); 
		p=NULL;
	}
	operator IDispatch*()
	{
		return p;
	}
	IDispatch& operator*()
	{
		ATLASSERT(p!=NULL); 
		return *p; 
	}
	IDispatch** operator&()
	{
		ATLASSERT(p==NULL);
		return &p;
	}
	IDispatch* operator->()
	{
		ATLASSERT(p!=NULL);
		return p;
	}
	IDispatch* operator=(IDispatch* lp)
	{
		return (IDispatch*)AtlComPtrAssign((IUnknown**)&p, lp);
	}
	IDispatch* operator=(IUnknown* lp)
	{
		return (IDispatch*)AtlComQIPtrAssign((IUnknown**)&p, lp, IID_IDispatch);
	}
	BOOL operator!(){
		return (p == NULL) ? TRUE : FALSE;
	}
		
	HRESULT GetPropertyByName(LPCOLESTR lpsz, CComDispatchDriverEx& dd)
	{
		HRESULT							hr;
		CComVariant						v;

		if (hr = GetPropertyByName(lpsz, &v)) return hr;
		if (v.vt == VT_DISPATCH){
			dd = v.pdispVal;
			return S_OK;
		} else {			
			ATLASSERT(false);
			return E_FAIL;
		}
	}

	HRESULT GetPropertyByName(LPCOLESTR lpsz, VARIANT* pVar)
	{
		ATLASSERT(p);
		ATLASSERT(pVar);
		DISPID dwDispID;
		HRESULT hr = GetIDOfName(lpsz, &dwDispID);
		if (SUCCEEDED(hr))
			hr = GetProperty(p, dwDispID, pVar);
		return hr;
	}

	HRESULT GetPropertyByName(OLECHAR FAR* szMember, VARIANT* pparams, unsigned int nParams, CComDispatchDriverEx& dd)
	{
		HRESULT							hr;
		CComVariant						v;

		if (hr = GetPropertyByName(szMember, pparams, nParams, &v)) return hr;
		if (v.vt == VT_DISPATCH){
			dd = v.pdispVal;
			return S_OK;
		} else {			
			ATLASSERT(false);
			return E_FAIL;
		}
	}

	HRESULT GetPropertyByName(OLECHAR FAR* szMember, VARIANT* pparams, unsigned int nParams, VARIANT* result) const
	{
		return GetPropertyWithParameters(szMember, *this, pparams, nParams, result);
	}
	
	HRESULT PutPropertyByName(const std::string& szMember, VARIANT* pVar)
	{			
		HRESULT							hr;
		OLECHAR*						sMember;		
		ULONG							ulCnt;
		
		ulCnt = MultiByteToWideChar(CP_ACP, MB_PRECOMPOSED, szMember.c_str(), -1, NULL, 0) + 1;	// see how big the unicode string will be	   
		sMember = new OLECHAR [ulCnt];	// allocate a buffer for the unicode string
		MultiByteToWideChar(CP_ACP, 0, szMember.c_str() , -1, sMember, ulCnt);	
		hr = PutPropertyByName(sMember, pVar);
		delete sMember;
		return hr;
	}

	HRESULT PutPropertyByName(LPCOLESTR lpsz, VARIANT* pVar)
	{
		return CComDispatchDriver::PutPropertyByName(lpsz, pVar);		
	}

	HRESULT PutPropertyByName(OLECHAR FAR* szMember, VARIANT* pparams, unsigned int nParams) const
	{
		return PutPropertyWithParameters(szMember, *this, pparams, nParams);
	}

	// Invoke a method with no parameters.
	HRESULT Invoke0(DISPID dispid, VARIANT* pvarRet = NULL)
	{
		return InvokeN(dispid, NULL, 0, pvarRet);
	}	
	HRESULT Invoke0(LPCOLESTR lpszName, VARIANT* pvarRet = NULL)
	{
		return InvokeN(lpszName, NULL, 0, pvarRet);		
	}

	// Invoke a method by name with a single parameter
	HRESULT Invoke1(LPCOLESTR lpszName, VARIANT* pvarParam1, CComDispatchDriverEx& dd)
	{
		HRESULT							hr;
		CComVariant						v;
		if (hr = Invoke1(lpszName,  pvarParam1, &v)) return hr;
		if (v.vt == VT_DISPATCH){
			dd = v.pdispVal;
			return S_OK;
		} else {			
			ATLASSERT(false);
			return E_FAIL;
		}
	}		
	HRESULT Invoke1(LPCOLESTR lpszName, VARIANT* pvarParam1, VARIANT* pvarRet = NULL)
	{
		HRESULT hr;
		DISPID dispid;
		hr = GetIDOfName(lpszName, &dispid);
		if (SUCCEEDED(hr)) hr = Invoke1(dispid, pvarParam1, pvarRet);
		return hr;
	}

	// Invoke a method by DISPID with a single parameter
	HRESULT Invoke1(DISPID dispid, VARIANT* pvarParam1, CComDispatchDriverEx& dd)
	{
		HRESULT							hr;
		CComVariant						v;
		if (hr = Invoke1(dispid, pvarParam1, &v)) return hr;
		if (v.vt == VT_DISPATCH){
			dd = v.pdispVal;
			return S_OK;
		} else {			
			ATLASSERT(false);
			return E_FAIL;
		}
	}
	HRESULT Invoke1(DISPID dispid, VARIANT* pvarParam1, VARIANT* pvarRet = NULL)
	{
		DISPPARAMS dispparams = { pvarParam1, NULL, 1, 0};
		return p->Invoke(dispid, IID_NULL, LOCALE_USER_DEFAULT, DISPATCH_METHOD, &dispparams, pvarRet, NULL, NULL);
	}
				
	// This implementation of InvokeN throws a string error for the DISP_E_EXCEPTION return case.
	// This facilitates easy capture of, say, a Visual Basic err.raise() call.
	HRESULT InvokeN(LPCOLESTR lpszName, VARIANT* pvarParams, int nParams, VARIANT* pvarRet = NULL)
	{
		HRESULT hr;
		DISPID dispid;
		hr = GetIDOfName(lpszName, &dispid);
		if (SUCCEEDED(hr))
			hr = InvokeN(dispid, pvarParams, nParams, pvarRet);
		return hr;
	}	
	HRESULT InvokeN(DISPID dispid, VARIANT* pvarParams, int nParams, VARIANT* pvarRet = NULL)
	{
		HRESULT			hr;
		EXCEPINFO		excepinfo;
		
		::memset(&excepinfo, NULL, sizeof(EXCEPINFO));
		DISPPARAMS dispparams = {pvarParams, NULL, nParams, 0};
		if ((hr =  p->Invoke(dispid, IID_NULL, LOCALE_USER_DEFAULT, DISPATCH_METHOD, &dispparams, pvarRet, &excepinfo, NULL)) == DISP_E_EXCEPTION){
			// This means that in this case an error raised inside a successfully invoked function.
			// This mechanism allows us to use "err.raise" in VBA class modules. The error text
			// propagates through here.
			// Microsoft state that excepinfo is populated for only this hr value.
			CComBSTR sDescription(excepinfo.bstrDescription);			
			if (!sDescription.Length()) sDescription = "Unspecified error occurred";
			::SysFreeString(excepinfo.bstrSource);		
			::SysFreeString(excepinfo.bstrDescription);
			::SysFreeString(excepinfo.bstrHelpFile);
			throw sDescription;
		}
		return hr;
	}
};

#endif
