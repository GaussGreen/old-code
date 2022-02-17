//	handle_deals.cpp : Implementation of CDeals
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqObjects.h"
#include "comobjectcollectionfunctions.h"
#include "siriusapplication.h"
#include "xmlstreamer.h"
#include "handle_assets.h"
#include "handle_correlationmatrix.h"
#include "handle_deal.h"
#include "handle_deals.h"
#include "handle_dividendschedule.h"
#include "handle_spotschedule.h"
#include "handle_zerocurve.h"
#include "handle_zerocurves.h"

STDMETHODIMP CDeals::Evaluate(BSTR Calculate, IResult** pVal)
{			
	return CComObjectCollectionFunctions<IDeal>(&m_coll, &m_ObjectToNotional).Evaluate(Calculate, pVal);
}

STDMETHODIMP CDeals::get_Date(DATE* pVal)
{
	begin_function
	*pVal = CComObjectCollectionFunctions<IDeal>(&m_coll).GetDate();	
	end_function
}

STDMETHODIMP CDeals::GetFxUnderlyings(IAssets** pVal)
{
	return CComObjectCollectionFunctions<IDeal>(&m_coll).GetFxUnderlyings(pVal);
}

STDMETHODIMP CDeals::GetUnderlyings(IAssets** pVal)
{
	return CComObjectCollectionFunctions<IDeal>(&m_coll).GetUnderlyings(pVal);
}

STDMETHODIMP CDeals::GetVolatilityStructures(IVolatilityStructures** pVal)
{
	return CComObjectCollectionFunctions<IDeal>(&m_coll).GetVolatilityStructures(pVal);
}

STDMETHODIMP CDeals::GetZeroCurves(IZeroCurves** pVal)
{	
	return CComObjectCollectionFunctions<IDeal>(CComQIPtr<IEvaluatable>(this)).GetZeroCurves(pVal);	
}

/*static*/ HRESULT CDeals::Load(const std::string& szBook, DataSourceEnum ds, long nDate, CComPtr<IDeals>& spDeals)
{
	HRESULT								hr;
	FileSystemEnum						fs;
		
	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();
	if (hr = spDeals.CoCreateInstance(CLSID_Deals)) return hr;							
	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		_ConnectionPtr						spConnection;
		std::map<CAdapt<CComBSTR>, CComPtr<IDeal> >	map;							// map of DealID to deal pointer
				
		// non-position part
		{
			std::stringstream					ssQuery;
			_RecordsetPtr						prs;							// ADO recordset
			CComVariant							vDealID;						// deal field ID
			CComVariant							vValue;							// deal field value						
			ssQuery << "sp_user_get_deals_main '" << szBook << "', " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "'";
			try {
				_Module.GetConnection(spConnection);
				prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
				if (prs->adoEOF) return CParameterMap::ReturnErrorRS(IDS_BOOK_ID_INVALID, szBook, IID_IDeals);
				while (!prs->adoEOF){			
					vDealID = prs->GetFields()->GetItem(0L)->GetValue();
					vValue = prs->GetFields()->GetItem(1L)->GetValue();
					if (vDealID.vt != VT_BSTR || vValue.vt != VT_BSTR) return E_FAIL;						
					// parse the XML into a variant, create a deal object with that variant value and add it to the deals collection
					CComPtr<IDeal>	spDeal;
					CComVariant		v;					
					if (hr = spDeal.CoCreateInstance(CLSID_Deal)) return hr;
					ATLASSERT(false);	// ToDo - why is the next line commented?!
					// if (hr = CXmlStreamer::GetVariant((char*)_bstr_t(vValue.bstrVal), v)) return hr;					
					if (hr = g_pApplication->GetObjectManager().ImplementPutValue(v, spDeal.p)) return hr;					
					if (hr = spDeals->Add(CComVariant(), spDeal)) return hr;
					// We need an explicit insertion because spDeals is not, in general, the maintained colletion
					g_pApplication->GetObjectManager().InsertObject(spDeal, false);	
					map[CComBSTR(vDealID.bstrVal)] = spDeal;	// we will need this later
					prs->MoveNext();
				}
			} catch (_com_error& e){
				throw estring(e);
			}
		}
				
		// position part
		{
			std::stringstream					ssQuery;
			_RecordsetPtr						prs;							// ADO recordset
			CComVariant							vDealID;
			CComVariant							vPositionID;
			CComVariant							vNotional;
			CComVariant							vValue;
			ssQuery << "sp_user_get_deals_ref '" << szBook << "', " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "'";
			try {
				prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
				while (!prs->adoEOF){
					// Note that the stored procuedure is set up to return the DealID's in order,
					// followed by the position and its details.
					vDealID = prs->GetFields()->GetItem(0L)->GetValue();
					vPositionID = prs->GetFields()->GetItem(1L)->GetValue();
					vNotional = prs->GetFields()->GetItem(2L)->GetValue();
					vValue = prs->GetFields()->GetItem(3L)->GetValue();
					if (vDealID.vt != VT_BSTR || vPositionID.vt != VT_BSTR || vNotional.vt != VT_R8 || vValue.vt != VT_BSTR) return E_FAIL;						
					CComPtr<IDeal> spDeal = map[CComBSTR(vDealID.bstrVal)];
					if (spDeal){ //an error is probably due to synchronisation - this does not matter
						CComPtr<IPositions>	spPositions;
						CComVariant			vPosition;						
						if (spDeal->get_Positions(&spPositions)) return E_FAIL;							
						ATLASSERT(false);	// ToDo - why is the next line commented?						
						//if (hr = CXmlStreamer::GetVariant((char*)_bstr_t(vValue.bstrVal), vPosition)) return hr;
						if (vPosition.vt != VT_DISPATCH) return E_FAIL;
						ATLASSERT(false);	// ToDo - do we need to insert vPosition.pdispVal into the maintained collection?						
						if (hr = spPositions->Add(CComVariant(), (IPosition*)vPosition.pdispVal)) return hr;
						ATLASSERT(false);	// ToDo - does put_Notional work with the parameter set to vPosition?
						if (hr = spPositions->put_Notional(vPosition, vNotional.dblVal)) return hr;
					}
					prs->MoveNext();
				}
			} catch (_com_error& e){
				throw estring(e);
			}
		}
		return S_OK;
	} else if (fs == fsNTFS){										
		return CComObjectCollectionFunctions<IDeal>().ImplementLoadCollection_NTFS(spDeals, IDS_HEADER_BOOK, szBook, _Module.GetLocation(), "", ds, nDate, CDeal::Load);
	}
	return E_FAIL;
}

STDMETHODIMP CDeals::PutDataSource(DataSourceEnum newVal)
{
	begin_function
	CComObjectCollectionFunctions<IDeal>(&m_coll).PutDataSource(newVal);
	end_function
}

STDMETHODIMP CDeals::PutDate(DATE newVal)
{
	begin_function
	CComObjectCollectionFunctions<IDeal>(&m_coll).PutDate(newVal);
	end_function
}
