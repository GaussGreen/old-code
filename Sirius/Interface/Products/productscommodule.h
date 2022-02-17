// productscommodule.h: interface for the CProductsComModule class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _PRODUCTSCOMMODULE_H
#define _PRODUCTSCOMMODULE_H

#pragma once
class estring;

class CProductsComModule : public CComModule  
{
public:
	CProductsComModule();
	virtual									~CProductsComModule(){}
	CComPtr<ISiriusApplication>				GetSiriusApplication(void) const;
	HRESULT									AddError(const estring& sz) const;
	void									Term(void);
		
private:
	mutable CComPtr<ISiriusApplication>		m_spApplication;
};

#endif
