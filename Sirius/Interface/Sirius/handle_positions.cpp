//	handle_positions.cpp : Implementation of CPositions
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqObjects.h"
#include "comobjectcollectionfunctions.h"
#include "siriusapplication.h"
#include "xmlstreamer.h"
#include "handle_assets.h"
#include "handle_correlationmatrix.h"
#include "handle_dividendschedule.h"
#include "handle_spotschedule.h"
#include "handle_position.h"
#include "handle_positions.h"
#include "handle_zerocurve.h"
#include "handle_zerocurves.h"

STDMETHODIMP CPositions::Evaluate(BSTR Calculate, IResult** pVal)
{			
	return CComObjectCollectionFunctions<IPosition>(&m_coll, &m_ObjectToNotional).Evaluate(Calculate, pVal);
}

STDMETHODIMP CPositions::get_Date(DATE* pVal)
{
	begin_function
	*pVal = CComObjectCollectionFunctions<IPosition>(&m_coll).GetDate();	
	end_function
}

STDMETHODIMP CPositions::GetFxUnderlyings(IAssets** pVal)
{
	return CComObjectCollectionFunctions<IPosition>(&m_coll).GetFxUnderlyings(pVal);
}

STDMETHODIMP CPositions::GetUnderlyings(IAssets** pVal)
{
	return CComObjectCollectionFunctions<IPosition>(&m_coll).GetUnderlyings(pVal);
}
	
STDMETHODIMP CPositions::GetVolatilityStructures(IVolatilityStructures** pVal)
{
	return CComObjectCollectionFunctions<IPosition>(&m_coll).GetVolatilityStructures(pVal);
}

STDMETHODIMP CPositions::GetZeroCurves(IZeroCurves** pVal)
{	
	return CComObjectCollectionFunctions<IPosition>(CComQIPtr<IEvaluatable>(this)).GetZeroCurves(pVal);	
}

/*static*/ HRESULT CPositions::Load(const std::string& szBook, DataSourceEnum ds, long nDate, CComPtr<IPositions>& spPositions)
{		
	HRESULT								hr;		
	FileSystemEnum						fs;
			
	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();
	if (hr = spPositions.CoCreateInstance(CLSID_Positions)) return hr;
							
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;
		std::string						szError;
		if (szBook.size()){
			szError = "No positions found in book '" + szBook + "' on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
		} else {
			szError = "No positions found on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
		}
		ssQuery << "sp_user_get_positions '" << szBook << "', " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "'";
		if (hr = CComObjectCollectionFunctions<IPosition>(&szError).ImplementLoadCollection_SQL(ssQuery, ds, spPositions)) return hr;
	} else if (fs == fsNTFS){		
		if (hr = CComObjectCollectionFunctions<IPosition>().ImplementLoadCollection_NTFS(spPositions, IDS_HEADER_BOOK, szBook, _Module.GetLocation(), "", ds, nDate, CPosition::Load)) return hr;
	} else {
		return E_FAIL;
	}
	
	// Set the data source to ds for every product in every position in spPositions
	long nCount = 0L;
	spPositions->get_Count(&nCount);
	for (CComVariant v = 1L; v.lVal < nCount; v.lVal++){
		CComPtr<IPosition> spPosition;
		spPositions->get_Item(v, &spPosition);
		CComPtr<IProducts> spProducts;
		spPosition->get_Products(&spProducts);
		CComObjectCollectionFunctions<IProduct>().PutDataSource(spProducts, ds);
	}
	return S_OK;
}

STDMETHODIMP CPositions::PutDataSource(DataSourceEnum newVal)
{
	begin_function
	CComObjectCollectionFunctions<IPosition>(&m_coll).PutDataSource(newVal);
	end_function
}

STDMETHODIMP CPositions::PutDate(DATE newVal)
{
	begin_function
	CComObjectCollectionFunctions<IPosition>(&m_coll).PutDate(newVal);
	end_function
}
