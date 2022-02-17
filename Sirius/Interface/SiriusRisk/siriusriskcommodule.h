// siriusriskcommodule.h: interface for the CSiriusRiskComModule class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _SIRIUSRISKCOMMODULE_H
#define _SIRIUSRISKCOMMODULE_H

#pragma once
class estring;

class CSiriusRiskComModule : public CComModule  
{
public:
	CSiriusRiskComModule();
	virtual									~CSiriusRiskComModule(){}
	HRESULT									AddError(const estring& sz) const;
	CComPtr<ISiriusApplication>				GetSiriusApplication(void) const;
	void									Term(void);

private:
	mutable CComPtr<ISiriusApplication>		m_spApplication;
};

#endif
