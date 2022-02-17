//	handle_deal.cpp : Implementation of CDeal
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionfunctions.h"
#include "handle_deal.h"
#include "siriusapplication.h"
#include "handle_position.h"
#include "MlEqDate.h"
#include "ComObjectCollectionSerialisableKey.h"

/*static*/ const long					CDeal::s_nMaxPositionsInDeal = 128;				// this is constrained by sp_user_add_deal

implement_member_variable_rekey(Deal, Deals, DATE, Date, long, Date);
implement_member_string_rekey(Deal, Deals, Name, Name);
implement_member_variable_rekey(Deal, Deals, DataSourceEnum, DataSource, DataSourceEnum, DataSource);
implement_serialisable(Deal);


STDMETHODIMP CDeal::Evaluate(BSTR Calculate, IResult** pVal)
{
	begin_function
	return m_h->Evaluate(Calculate).CopyTo(pVal);
	end_function
}

STDMETHODIMP CDeal::GetFxUnderlyings(IAssets** pVal)
{
	begin_function
	return m_h->GetPositions()->GetFxUnderlyings(pVal);	
	end_function
}

STDMETHODIMP CDeal::get_PositionNotional(IPosition* Position, double *pVal)
{
	begin_function
	if (!m_h->GetPositions()) return E_POINTER;
	return m_h->GetPositions()->get_Notional(CComVariant(Position), pVal);
	end_function;
}

STDMETHODIMP CDeal::get_Positions(IPositions** pVal)
{
	begin_function
	return m_h->GetPositions().CopyTo(pVal);
	end_function
}

//	returns a collection of assets on which the deal depends
STDMETHODIMP CDeal::GetUnderlyings(IAssets** pVal)
{
	begin_function
	return m_h->GetPositions()->GetUnderlyings(pVal);
	end_function
}

STDMETHODIMP CDeal::get_Value(VARIANT* pVal)
{			
	return g_pApplication->GetObjectManager().ImplementGetValue(dynamic_cast<IDeal*>(this), pVal);
}

STDMETHODIMP CDeal::GetVolatilityStructures(IVolatilityStructures** pVal)
{
	begin_function
	return m_h->GetPositions()->GetVolatilityStructures(pVal);
	end_function
}

STDMETHODIMP CDeal::GetZeroCurves(IZeroCurves** pVal)
{	
	begin_function
	return m_h->GetPositions()->GetZeroCurves(pVal);
	end_function
}

STDMETHODIMP CDeal::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = {&IID_IDeal};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

/*static*/ HRESULT CDeal::Load(const std::string& szName, DataSourceEnum ds, long nDate, CComPtr<IDeal>& spDeal)
{		
	FileSystemEnum						fs;	
		
	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();
		
	// load the XML for the deal
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		
		HRESULT							hr;	
		_ConnectionPtr					spConnection;
		
		if (hr = spDeal.CoCreateInstance(CLSID_Deal)) return hr;
		// Set all the parts of the deal apart from the positions collection.
		{
			std::stringstream			ssQuery;						// SQL query to execute
			_RecordsetPtr				prs;							// ADO recordset
			CComVariant					v;
			CComVariant					vField;			
			ssQuery << "sp_user_get_deal_main '" << estring(szName) << "', " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "'";
			try {
				_Module.GetConnection(spConnection);
				prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
			} catch (_com_error& e){
				throw estring(e);
			}
			// the first recordset contains the deal non-position XML data				
			if (prs->adoEOF) return CParameterMap::ReturnErrorRS(IDS_DEAL_ID_INVALID, szName, IID_IDeal);
			vField = prs->GetFields()->GetItem(CComVariant(0, VT_I4))->GetValue();
			if (vField.vt != VT_BSTR) return CParameterMap::ReturnErrorR(IDS_UNHANDLED_EXCEPTION);
			// parse the XML into a variant
			if (hr = CXmlStreamer::GetVariant((char*)_bstr_t(vField.bstrVal), v)) return hr;			
			if (hr = g_pApplication->GetObjectManager().ImplementPutValue(v, spDeal.p)) return hr;
			spDeal->put_DataSource(ds);
		}
		
		// Set the positions collection.
		CComPtr<IPositions>				spPositions;
		if (hr = spPositions.CoCreateInstance(CLSID_Positions)) return hr;	
		{
			CComPtr<IPositions>			spMaintainedPositions;
			std::stringstream			ssQuery;						// SQL query to execute
			_RecordsetPtr				prs;							// ADO recordset
						
			if (hr = _Module.GetSiriusApplication()->get_Positions(&spMaintainedPositions)) return hr;
			ssQuery << "sp_user_get_deal_ref '" << estring(szName) << "', " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "'";
			try {
				prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
			} catch (_com_error& e){
				throw estring(e);
			}		
			while (!prs->adoEOF){				
				CComPtr<IPosition>			spPosition;
				CComVariant					vPosition;
				CComVariant					vNotional;
				
				vPosition = prs->GetFields()->GetItem(0L)->GetValue();
				vNotional = prs->GetFields()->GetItem(1L)->GetValue();

				// Note that the next line emulates the NTFS case where positions are only loaded
				// if not already in memory.
				if (spMaintainedPositions->get_Item(CComObjectCollectionSerialisableKey(estring(vPosition), nDate, ds), &spPosition)) propagate_error;
				if (spPositions->Add(CComVariant(), spPosition)) propagate_error;
				if (spPositions->put_Notional(CComVariant(spPosition.p), vNotional.dblVal)) propagate_error;
				prs->MoveNext();
			}
		}
		return spDeal->put_Positions(spPositions);
	} else if (fs == fsNTFS){
		CComObjectCollectionFunctions<IDeal>::ImplementLoad_NTFS(IID_IDeal, _Module.GetLocation(), "", szName, ds, nDate, spDeal);
		return S_OK;
	}
	return E_FAIL;	
}

STDMETHODIMP CDeal::PutDataSource(/*[in]*/ DataSourceEnum newVal)
{
	begin_function
	return m_h->GetPositions()->PutDataSource(newVal);
	end_function
}

STDMETHODIMP CDeal::PutDate(DATE newVal)
{
	begin_function
	return m_h->GetPositions()->PutDate(newVal);
	end_function
}

STDMETHODIMP CDeal::put_PositionNotional(IPosition* Position, double newVal)
{
	begin_function
	if (!m_h->GetPositions()) return E_POINTER;
	return m_h->GetPositions()->put_Notional(CComVariant(Position), newVal);
	end_function
}

STDMETHODIMP CDeal::put_Positions(IPositions* newVal)
{
	begin_function
	m_h->PutPositions(newVal);
	end_function
}

STDMETHODIMP CDeal::put_Value(VARIANT newVal)
{	
	begin_function
	HRESULT hr = g_pApplication->GetObjectManager().ImplementPutValue(newVal, dynamic_cast<IDeal*>(this));
	return hr;
	end_function	
}

STDMETHODIMP CDeal::Save(BSTR* pVal)
{
	HRESULT								hr;
	FileSystemEnum						fs;		
	long								nCurrentDate = MlEqDate::GetCurrentDate();
				
	begin_function
	check_publishing_enabled
	if (!m_h->GetDate()) throw "You cannot save the deal '" + m_h->GetName() + "' since it does not have a date.";
	
	// Prevent the saving of deals where the positions dates and data sources do not match the corresponding deal properties.
	{
		long nPositions = 0L;
		m_h->GetPositions()->get_Count(&nPositions);
		for (long nPosition = 1; nPosition < nPositions; nPosition++){
			CComPtr<IPosition>  spPosition;
			DATE				date;
			DataSourceEnum		ds;
			CComBSTR			sName;
			
			if (m_h->GetPositions()->get_Item(CComVariant(nPosition), &spPosition)) continue;
			spPosition->get_DataSource(&ds);
			spPosition->get_Date(&date);
			spPosition->get_Name(&sName);

			long				nDate = (long)date;
			if (!nDate) throw "You cannot save the deal '" + m_h->GetName() + "' since one of its positions ('" + estring(sName) + "') does not have a date.";
			if (nDate != m_h->GetDate()) throw "The position '" + estring(sName) + "' has the date set to " + MlEqDate(nDate).GetString() + " which differs from the deal date " + MlEqDate(m_h->GetDate()).GetString();
			if (ds != m_h->GetDataSource()) throw "The position '" + estring(sName) + "' has the data source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "' which differs from the deal data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "'";
		}
	}

	// All the positions must have the same date as the deal or we throw an error.
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		CParameterMap						pmThis;							// representation of the non-position part of this object
		estring								szPositionList;					// comma separated list of position names
		estring								szNotionalList;					// comma separated list of notional values	
				
		if (m_h->GetDataSource() != Last) throw "Deals with data sources other than 'Last' cannot be saved. The deal '" + m_h->GetName() + "' has the source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "'";
		
		// To start with, save all the positions in the deal.
		if (hr = m_h->GetPositions()->Save(NULL)) return hr;

		// Get the prototype pmThis.
		{
			CComVariant						v;
			DataSourceEnum					ds = m_h->GetDataSource();
			// All saved records need to be data source invariant.
			m_h->PutDataSource(NoDataSource);
			if (hr = g_pApplication->GetObjectManager().ImplementGetValue(dynamic_cast<IDeal*>(this), &v) || pmThis.SetValue(v)) return hr;
			m_h->PutDataSource(ds);
		}

		// Build szPositionList and szNotionalList; removing all the position data out of the deal.
		{
			std::set<estring>				setNames;						// use this to check for duplicates and the number of positions in the deal
			long							nRow = pmThis.GetRows() - 1;
			while (nRow >= 0){
				CComPtr<IPosition>	spPosition;
				CComVariant			vHandle;
				if (hr = pmThis.GetValue(nRow, 0, &vHandle)) return hr;
				if (!g_pApplication->GetObjectManager().GetObject(vHandle, spPosition)){
					// This row is a handle.
					CComBSTR		sName;
					double			fNotional;				
					long			nDate;
					
					if (spPosition->get_Name(&sName)) propagate_error;
					{
						DATE			date;
						if (spPosition->get_Date(&date)) propagate_error;
						nDate = (long)date;
					}
					if (hr = pmThis.GetValue(nRow, 1, &fNotional)) return hr;
					
					estring			szName(sName);
					if (setNames.find(szName) != setNames.end()){
						return CParameterMap::ReturnErrorRS(IDS_DUPLICATE_POSITION_ID, szName);
					} else {
						setNames.insert(szName);
						if (setNames.size() > s_nMaxPositionsInDeal){
							return CParameterMap::ReturnErrorRS(IDS_TOO_MANY_POSITIONS_IN_DEAL, estring(s_nMaxPositionsInDeal));
						}
					}
					
					if (szPositionList.size()) szPositionList += ",";
					szPositionList += sName;				
					if (szNotionalList.size()) szNotionalList += ",";
					szNotionalList += fNotional;
					if (hr = pmThis.RemoveRow(nRow)) return hr;
					if (nRow == pmThis.GetRows()) nRow--;
				} else {
					// this row is not a handle - leave it alone!
					nRow--;
				}
			}
		}
		
		// Save the deal shell.	
		xmlstreamer						ssXML;
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;						// SQL query to execute
		_RecordsetPtr					prs;							// ADO recordset				

		if (hr = CXmlStreamer::GetXML(pmThis.GetValue(), ssXML)) return hr;
		ssQuery << "sp_user_add_deal '" << m_h->GetName() << "', " << m_h->GetDate() << ", '" << m_h->GetBook() << "', '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) << "', '" << _Module.GetLocation() << "', '" << szPositionList << "', '" << szNotionalList << "', '" << (char*)ssXML << "'";
		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
		} catch (_com_error& e){
			throw estring(e);
		}
		
		// Set the return type the returned value of the stored procedure.
		CComVariant vValue(prs->GetFields()->GetItem(CComVariant(0, VT_I4))->GetValue());
		if (vValue.vt != VT_BSTR){
			ATLASSERT(false);
			return CParameterMap::ReturnErrorR(IDS_UNHANDLED_EXCEPTION);
		}
		if (nCurrentDate != MlEqDate::GetCurrentDate()) throw "The date changed during the save deal process. This might have created unexpected position / deal relationships. You should consider resaving this deal.";
		return CComBSTR(vValue.bstrVal).CopyTo(pVal);
	} else if (fs == fsNTFS){		
		xmlstreamer							ssXML;		
		// We relax the condition that only 'Last' data sources can be saved (or the UpdatePL process - if implemented for deals - gets tricky).
		if (m_h->GetDataSource() == NoDataSource) throw "The deal '" + m_h->GetName() + "' has the data source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "' and so cannot be saved";		
		if (hr = m_h->GetPositions()->Save(NULL)) return hr;
		if (hr = CXmlStreamer::GetXML(CComPtr<IDeal>(this), ssXML)) return hr;
		if (hr = CComObjectCollectionFunctions<IDeal>::ImplementSave_NTFS(ssXML, IID_IDeal, _Module.GetLocation(), m_h->GetName(), m_h->GetDataSource(), m_h->GetDate(), pVal)) return hr;
		if (nCurrentDate != MlEqDate::GetCurrentDate()) throw "The date changed during the save deal process. This might have created unexpected position / deal relationships. You should consider resaving this deal.";
		return S_OK;		
	} else {
		return E_FAIL;
	}
	end_function
}