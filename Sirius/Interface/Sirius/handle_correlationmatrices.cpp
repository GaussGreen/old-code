//	handle_correlationmatrices.cpp : Implementation of CCorrelationMatrices
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_correlationmatrices.h"
#include "handle_correlationmatrix.h"

void CCorrelationMatrices::CheckAddHandle(std::string* pszHandle)
{
	bool								bIsMaintainedCollection = g_pApplication->GetObjectManager().IsMaintainedCollection(CComPtr<ICorrelationMatrices>(this));
	if (bIsMaintainedCollection) pszHandle->erase();
}

/*static*/ HRESULT CCorrelationMatrices::Load(const std::string& szDummy, DataSourceEnum ds, long nDate, CComPtr<ICorrelationMatrices>& spCorrelationMatrices)
{	
	// This function loads the blank correlation matrix (which will contain
	// all the correlations) and adds this single element to a new collection.
	HRESULT								hr;
	CComPtr<ICorrelationMatrix>			spCorrelationMatrix;
	
	if (szDummy.size()) throw "The identifier value '" + szDummy + "' is invalid for the load correlation matrices request";
	if (hr = CCorrelationMatrix::Load("", ds, nDate, spCorrelationMatrix)) return hr;
	spCorrelationMatrices.CoCreateInstance(CLSID_CorrelationMatrices);
	if (hr = spCorrelationMatrices->Add(CComVariant(CCorrelationMatrix::s_szCorrelationMatrixName), spCorrelationMatrix)) return hr;
	return S_OK;
}

// See comments in "comobjectcollectionserialisable.h" for return value conventions.
// For the maintained correlation matrix collection, we are reluctant to remove existing
// correlation pairs when adding a new object. What we do is merge in the data in the input 
// spSingular into any existing item at that point (i.e. date and data source in the maintained
// collection.)
// 
HRESULT CCorrelationMatrices::SpecialAdd(CComPtr<ICorrelationMatrix> spSingular, const CComBSTR& sKey)
{
	bool								bIsMaintainedCollection = g_pApplication->GetObjectManager().IsMaintainedCollection(CComPtr<ICorrelationMatrices>(this));
	KeyToHandleType::const_iterator		itKeyToHandle = m_KeyToHandle.find(sKey);

	if (!bIsMaintainedCollection) return S_FALSE;
	if (itKeyToHandle == m_KeyToHandle.end()) return S_FALSE;

	// We need to merge in the correlation matrix spSingular with the correlation matrix.
	CComPtr<ICorrelationMatrix> spGlobal = dynamic_cast<ICorrelationMatrix*>(m_coll[itKeyToHandle->second].pdispVal);
	ATLASSERT(spGlobal);
	if (spGlobal->Add(spSingular, VARIANT_FALSE)) ATLASSERT(false);
	return S_OK;
}