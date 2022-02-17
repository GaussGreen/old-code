// temporaryproductscommodule.cpp: implementation of the CTemporaryProductsComModule class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "temporaryproductscommodule.h"

CTemporaryProductsComModule::CTemporaryProductsComModule()
{
}

HRESULT CTemporaryProductsComModule::AddError(const estring& sz) const
{
	ATLASSERT(GetSiriusApplication());	
	return GetSiriusApplication()->AddError((CComBSTR)sz);
}

// get a pointer to the singleton sirius application object; creating it if necessary
CComPtr<ISiriusApplication> CTemporaryProductsComModule::GetSiriusApplication(void) const
{
	if (!m_spApplication){
		if (m_spApplication.CoCreateInstance(_CLSID_SiriusApplication)){
			CParameterMap::DisplayError(L"Fatal error encountered when creating the global Sirius Application object!\n\nYou cannot use Sirius in this state.", MB_ICONSTOP);
		}
	}
	return m_spApplication;
}

void CTemporaryProductsComModule::Term(void)
{
	m_spApplication = NULL;
	CComModule::Term();
}