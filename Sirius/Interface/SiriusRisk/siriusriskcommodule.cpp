// siriusriskcommodule.cpp: implementation of the CSiriusRiskComModule class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "siriusriskcommodule.h"

CSiriusRiskComModule::CSiriusRiskComModule()
{
}

HRESULT CSiriusRiskComModule::AddError(const estring& sz) const
{
	ATLASSERT(GetSiriusApplication());	
	return GetSiriusApplication()->AddError((CComBSTR)sz);
}

// get a pointer to the singleton sirius application object; creating it if necessary
CComPtr<ISiriusApplication> CSiriusRiskComModule::GetSiriusApplication(void) const
{
	if (!m_spApplication){
		if (m_spApplication.CoCreateInstance(_CLSID_SiriusApplication)){
			CParameterMap::DisplayError(L"Fatal error encountered when creating the global Sirius Application object!\n\nYou cannot use Sirius in this state.", MB_ICONSTOP);
		}
	}
	return m_spApplication;
}

void CSiriusRiskComModule::Term(void)
{
	AfxDaoTerm();
	m_spApplication = NULL;
	CComModule::Term();
}