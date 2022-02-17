// apacproductscommodule.h: interface for the CAPACProductsComModule class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _APACPRODUCTSCOMMODULE_H
#define _APACPRODUCTSCOMMODULE_H

#pragma once
class estring;

class CAPACProductsComModule : public CComModule  
{
public:
	CAPACProductsComModule();
	virtual									~CAPACProductsComModule(){}
	HRESULT									AddError(const estring& sz) const;
	CComPtr<ISiriusApplication>				GetSiriusApplication(void) const;
	void									Term(void);

private:
	mutable CComPtr<ISiriusApplication>		m_spApplication;
};

#endif
