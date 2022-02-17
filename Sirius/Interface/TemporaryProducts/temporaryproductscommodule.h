// temporaryproductscommodule.h: interface for the CTemporaryProductsComModule class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _TEMPORARYPRODUCTSCOMMODULE_H
#define _TEMPORARYPRODUCTSCOMMODULE_H

#pragma once
class estring;

class CTemporaryProductsComModule : public CComModule  
{
public:
	CTemporaryProductsComModule();
	virtual									~CTemporaryProductsComModule(){}
	HRESULT									AddError(const estring& sz) const;
	CComPtr<ISiriusApplication>				GetSiriusApplication(void) const;
	void									Term(void);

private:
	mutable CComPtr<ISiriusApplication>		m_spApplication;
};

#endif
