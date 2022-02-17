//	progress.cpp : Implementation of CProgress
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "progress.h"
#include "excelinterface.h"

HRESULT CProgress::GetWindowText(std::string* psz) const
{
	if (!m_hWnd) return E_FAIL;
	int									nLength = GetWindowTextLength(m_hWnd);
	LPTSTR								lpString = new TCHAR[nLength + 1];	
	
	if (::GetWindowText(m_hWnd, lpString, nLength + 1)){
		psz->assign(lpString);
	}	
	delete lpString;	
	return S_OK;
}

STDMETHODIMP CProgress::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IProgress };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

HRESULT CProgress::FinalConstruct()
{	
	if (!(m_hWnd = CExcelInterface::GetWindow())){	
		m_hWnd = ::GetActiveWindow();
	}
	m_nPercent = 0;
	return GetWindowText(&m_szInitial);	
}

HRESULT CProgress::FinalRelease()
{
	if (m_szInitial.size() && m_hWnd) ::SetWindowText(m_hWnd, m_szInitial.c_str());
	return S_OK;
}

STDMETHODIMP CProgress::get_PercentCompleted(short *pVal)
{
	*pVal = m_nPercent;
	return S_OK;
}

STDMETHODIMP CProgress::put_PercentCompleted(short newVal)
{
	short int							nPercent = __max(__min(newVal, 100), 0);	
	if (m_nPercent == nPercent) return S_OK;
	m_nPercent = nPercent;
	estring sz(m_szInitial + " (" + estring(m_nPercent) + "%)");
	
	if (m_hWnd) ::SetWindowText(m_hWnd, sz.c_str());
	return S_OK;
}