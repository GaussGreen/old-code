//	handle_products.cpp : Implementation of CProducts
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqObjects.h"
#include "comobjectcollectionfunctions.h"
#include "handle_products.h"
#include "handle_position.h"
#include "handle_assets.h"
#include "handle_correlationmatrix.h"
#include "handle_deals.h"
#include "handle_dividendschedule.h"
#include "handle_spotschedule.h"
#include "handle_zerocurve.h"
#include "handle_zerocurves.h"


STDMETHODIMP CProducts::Evaluate(BSTR Calculate, IResult** pVal)
{
	return CComObjectCollectionFunctions<IProduct>(&m_coll, &m_ObjectToNotional).Evaluate(Calculate, pVal);
}

STDMETHODIMP CProducts::get_Date(DATE* pVal)
{
	begin_function
	*pVal = CComObjectCollectionFunctions<IProduct>(&m_coll).GetDate();	
	end_function
}

STDMETHODIMP CProducts::GetFxUnderlyings(IAssets** pVal)
{	
	return CComObjectCollectionFunctions<IProduct>(&m_coll).GetFxUnderlyings(pVal);
}

STDMETHODIMP CProducts::GetUnderlyings(IAssets** pVal)
{	
	return CComObjectCollectionFunctions<IProduct>(&m_coll).GetUnderlyings(pVal);
}

STDMETHODIMP CProducts::GetVolatilityStructures(IVolatilityStructures** pVal)
{
	return CComObjectCollectionFunctions<IProduct>(&m_coll).GetVolatilityStructures(pVal);
}

STDMETHODIMP CProducts::GetZeroCurves(IZeroCurves** pVal)
{	
	return CComObjectCollectionFunctions<IProduct>(CComQIPtr<IEvaluatable>(this)).GetZeroCurves(pVal);	
}

STDMETHODIMP CProducts::PutDataSource(DataSourceEnum newVal)
{
	begin_function
	CComObjectCollectionFunctions<IProduct>(&m_coll).PutDataSource(newVal);
	end_function
}

STDMETHODIMP CProducts::PutDate(DATE newVal)
{
	begin_function
	CComObjectCollectionFunctions<IProduct>(&m_coll).PutDate(newVal);
	end_function
}
