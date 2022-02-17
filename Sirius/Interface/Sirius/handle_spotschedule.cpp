//	handle_spotschedule.cpp : Implementation of CSpotSchedule
//
//	Author:					  David Cuin
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionfunctions.h"
#include "handle_spotschedule.h"
#include "siriusapplication.h"
#include "xmlstreamer.h"
#include "ComObjectCollectionSerialisableKey.h"

implement_member_variable_rekey(SpotSchedule, SpotSchedules, DATE, Date, long, Date);
implement_member_variable_rekey(SpotSchedule, SpotSchedules, DataSourceEnum, DataSource, DataSourceEnum, DataSource);
implement_serialisable(SpotSchedule);

STDMETHODIMP CSpotSchedule::Clear()
{
	m_h->clear();
	return S_OK;
}

HRESULT	CSpotSchedule::FinalConstruct(void)
{
	return S_OK;
}

STDMETHODIMP CSpotSchedule::get_Dates(VARIANT *pVal)
{	
	HRESULT								hr;
	std::vector<long>					vector;		
	CParameterMap						pmDates;
	
	begin_function
	m_h->GetDates(vector);
	if (hr = pmDates.SetValue(vector)) return hr;
	return pmDates.GetValue(pVal);
	end_function
}

STDMETHODIMP CSpotSchedule::get_FirstDate(DATE *pVal)
{		
	begin_function		
	*pVal = m_h->GetFirstDate();	
	end_function
}

STDMETHODIMP CSpotSchedule::get_FirstValue(double *pVal)
{
	begin_function
	*pVal = m_h->GetFirstValue();
	end_function	
}

STDMETHODIMP CSpotSchedule::get_LastDate(DATE *pVal)
{	
	begin_function
	*pVal = m_h->GetLastDate();
	end_function
}

STDMETHODIMP CSpotSchedule::get_LastValue(double *pVal)
{
	begin_function
	*pVal = m_h->GetLastValue();
	end_function	
}

STDMETHODIMP CSpotSchedule::get_Name(BSTR* pVal)
{
	begin_function
	estring sz = m_h->GetName();
	return sz.GetBSTR(pVal);
	end_function
}

STDMETHODIMP CSpotSchedule::get_Value(VARIANT *pVal)
{	
	begin_function
	HRESULT								hr;
	std::string							szDataSource = CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource());
	std::vector<long>					vectorDates;
	std::vector<double>					vectorValues;
	CParameterMap						pmDates, pmValues;
		
	unmap_parameter(m_h->GetName(), pmName);
	unmap_parameter(szDataSource, pmDataSource);
	unmap_parameter(m_h->GetDate(), pmDate);
	m_h->GetDatesAndValues(vectorDates, vectorValues);			
	if (hr = pmDates.SetValue(vectorDates)) return hr;	
	if (hr = pmValues.SetValue(vectorValues)) return hr;	
	return CParameterMap::VariableArgumentListToArray(pVal, 5, pmName, pmDate, pmDataSource, pmDates, pmValues);
	end_function
}

STDMETHODIMP CSpotSchedule::get_ValueAt(DATE Date, double *pVal)
{	
	begin_function
	*pVal = m_h->GetValueAt(Date);
	end_function	
}

STDMETHODIMP CSpotSchedule::get_Values(VARIANT *pVal)
{
	HRESULT								hr;
	std::vector<double>					vector;
	CParameterMap						pmValues;	

	begin_function
	m_h->GetValues(vector);
	if (hr = pmValues.SetValue(vector)) return hr;	
	return pmValues.GetValue(pVal);
	end_function
}

STDMETHODIMP CSpotSchedule::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_ISpotSchedule };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

/*static*/ HRESULT CSpotSchedule::Load(const std::string& szIdentifier, DataSourceEnum ds, long nDate, CComPtr<ISpotSchedule>& spSpotSchedule)
{
	begin_function
	HRESULT								hr;	
	FileSystemEnum						fs;	
	bool								bDateDefaulted = false;			// true if the input date, nDate, was zero so we have to default to the current date		
	long								nLastBwlDate = 0;				// Most recent spot schedule value obtained from bwl.

	// The meaning of nDate is (normally) to return a spot schedule where the
	// value at nDate is acquired from Sirius and values before nDate are
	// acquired from Bwl. Also, the sirius data source is used to acquire data
	// between the last value before nDate that bwl returns and nDate. The spot schedule
	// is trimmed so any data after nDate are removed.

	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate){
		bDateDefaulted = true;
		nDate = CComObjectCollectionSerialisableDefaulter::GetDate();
	}

	// USD.USD case
	if (!estring::CompareNoCase(szIdentifier, "USD.USD")){
		// create a schedule with a single value of 1.0 on the current date
		CComPtr<ISpotSchedule>		sp;
		if (hr = sp.CoCreateInstance(CLSID_SpotSchedule)) return hr;
		if (hr = sp->put_ValueAt(nDate, 1.0)) return hr;
		map_com_to_analytic(sp, SpotSchedule, h);
		h->PutName(szIdentifier);
		h->PutDataSource(ds);
		h->PutDate(nDate);
		spSpotSchedule = sp;
		return S_OK;
	}
	
	// Create the primitive spSpotSchedule by loading data from Bwl.
	if (_Module.GetUseBwlForSpotSchedule()){
		CComPtr<IDispatch> spObject;	
		try {
			_Module.LoadFromBwl(CLSID_SpotSchedule, szIdentifier, spObject);
		} catch (...){
			// do nothing
		}
		spSpotSchedule = dynamic_cast<ISpotSchedule*>(spObject.p);		
	}
	if (!spSpotSchedule) spSpotSchedule.CoCreateInstance(CLSID_SpotSchedule);
	map_com_to_analytic(spSpotSchedule, SpotSchedule, hSpotSchedule);	
	hSpotSchedule->RemoveOnOrAfter(nDate);
	try {
		// To compute nLastBwlDate we can't use MlEqSchedule::GetLastDate since it could attempt to extrapolate.
		// We are interested in the real, not any inferred, data.
		if (!hSpotSchedule->size()){
			nLastBwlDate = 0;
		} else {
			std::map<long, double>::const_iterator it = hSpotSchedule->end();
			nLastBwlDate = (--it)->first;			
		}
	} catch (...){
		nLastBwlDate = 0;
	}	
	
	// Get the rest of the spot schedule data from the native Sirius database.
	if ((fs = _Module.GetFileSystem()) == fsNTFS){
		// Spot schedules are not stored historically even though they are, in a sense, historic.
		CComPtr<ISpotSchedule> spSiriusSpotSchedule;
		CComObjectCollectionFunctions<ISpotSchedule>::ImplementLoad_NTFS(IID_ISpotSchedule, _Module.GetLocation(), "", szIdentifier, ds, 0, spSiriusSpotSchedule);
		map_com_to_analytic(spSiriusSpotSchedule, SpotSchedule, hSiriusSpotSchedule);
		hSpotSchedule->Merge(*hSiriusSpotSchedule);
	} else if (fs == fsSQLServer){
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;						// SQL query to execute
		_RecordsetPtr					prs;							// ADO recordset
										
		if (nLastBwlDate){
			// get the range nLastBwlDate + 1 to nDate from Sirius
			ssQuery << "sp_user_get_spots '" << szIdentifier << "', '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "', " << nLastBwlDate + 1 << ", " << nDate;
		} else {
			// get all the dates from 0 to nDate from Sirius
			ssQuery << "sp_user_get_spots '" << szIdentifier << "', '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "', " << 0 << ", " << nDate;
		}
		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
			while (!prs->adoEOF){
				long nDate = prs->GetFields()->GetItem(0L)->GetValue().lVal;
				double fValue = prs->GetFields()->GetItem(1L)->GetValue().dblVal;
				hSpotSchedule->insert(std::pair<long, double>(nDate, fValue));	// note that this overwrites any values loaded from bwl
				prs->MoveNext();
			}
			if (!hSpotSchedule->size()){
				throw "Cannot find any spot values for '" + szIdentifier + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
			}
			if (!bDateDefaulted && !hSpotSchedule->HasValueAt(nDate)){
				// If the caller has requested a specific date and that date
				// is not available then we should error.				
				throw "Cannot find spot value for '" + szIdentifier + "' on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
			}			
		} catch (_com_error& e){
			throw estring(e);
		}
	} else {
		return E_FAIL;
	}

	hSpotSchedule->PutName(szIdentifier);
	hSpotSchedule->PutDataSource(ds);
	hSpotSchedule->PutDate(nDate);
	end_function
}

STDMETHODIMP CSpotSchedule::Multiply(/*[in]*/ double Amount)
{
	begin_function
	m_h->Multiply(Amount);
	end_function
}

STDMETHODIMP CSpotSchedule::put_FirstValue(double newVal)
{
	begin_function	
	m_h->PutFirstValue(newVal);	
	end_function
}

STDMETHODIMP CSpotSchedule::put_LastValue(double newVal)
{
	begin_function	
	m_h->PutLastValue(newVal);	
	end_function
}

STDMETHODIMP CSpotSchedule::put_Value(VARIANT newVal)
{	
	/*parameter list is 0  - Name
						1  - DateOpt
					    2  - DataSource
						3+ - Schedule Data*/
							
	begin_function
	HRESULT								hr;
	std::vector<CParameterMap>			vpm;
	
	if (hr = CParameterMap::ArrayToVector(newVal, NULL, &vpm)) return hr;
	if (vpm.size() < 4) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_ISpotSchedule);	
	map_parameter(vpm[0], estring, Name);
	map_optional_parameter(vpm[1], long, DateOpt, 0L);	// we do the defaulting later
	map_enum_parameter(vpm[2], DataSourceEnum, DataSource)
	
	// the remaining elements in vpm make up the schedule
	for (long n = 4; n < vpm.size(); n++){
		if (vpm[3].GetRows() != vpm[n].GetRows()) return CParameterMap::ReturnErrorR(IDS_VECTOR_SIZES_DIFFER, IID_ISpotSchedule);
		vpm[3].AddToRHS(vpm[n]);
	}
	if (vpm[3].GetCols() > 2) return CParameterMap::ReturnErrorR(IDS_TOO_MANY_COLUMNS, IID_ISpotSchedule);
	if (vpm[3].GetCols() < 2) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_ISpotSchedule);
	
	// set the member variables
	if (hr = vpm[3].GetValue(0, 1, *m_h)) return hr;
	m_h->PutDataSource(DataSource);
	m_h->PutDate(DateOpt ? DateOpt : MlEqDate::GetCurrentDate());
	m_h->PutName(Name);
	end_function
}

STDMETHODIMP CSpotSchedule::put_ValueAt(DATE Date, double newVal)
{		
	(*m_h)[Date] = newVal;	
	return S_OK;
}

STDMETHODIMP CSpotSchedule::Save(BSTR *pVal)
{
	begin_function	
	HRESULT								hr;	
	xmlstreamer							ssXML;							// XML representation of vData			
	FileSystemEnum						fs = _Module.GetFileSystem();
	
	check_publishing_enabled
	if (!m_h->GetName().CompareNoCaseAndSpace("USD.USD")) return CParameterMap::ReturnErrorR(IDS_CANNOT_SAVE_USD_SPOT, IID_ISpotSchedule);	
	if (!m_h->GetDate()) throw "You cannot save the spot schedule '" + m_h->GetName() + "' since it does not have a date.";
	
	if (fs == fsSQLServer){
		if (m_h->GetDataSource() != Last) throw "Zero curves with data sources other than 'Last' cannot be saved. The zero curve '" + m_h->GetName() + "' has the source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "'";
		
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;
		_RecordsetPtr					prs;		
		ssQuery << "sp_user_add_spot '" << m_h->GetName() << "', " << m_h->GetDate() << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) << "', '" << _Module.GetLocation() << "', " << m_h->GetValueAt(m_h->GetDate());
		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
		} catch (_com_error& e){
			throw estring(e);
		}
		
		// set the return type the returned value of the stored procedure
		CComVariant vValue(prs->GetFields()->GetItem(0L)->GetValue());
		if (vValue.vt != VT_BSTR){
			ATLASSERT(false);
			return CParameterMap::ReturnErrorR(IDS_UNHANDLED_EXCEPTION);
		}
		if (pVal) return CComBSTR(vValue.bstrVal).CopyTo(pVal);
		return S_OK;		
	} else if (fs == fsNTFS){					
		// Spot schedules are not stored historically even though they are, in a sense, historic.
		CComPtr<ISpotSchedule>	spSpotSchedule;
		MlEqSpotScheduleHandle	hSpotSchedule;
		try {
			CSpotSchedule::Load(m_h->GetName(), m_h->GetDataSource(), m_h->GetDate(), spSpotSchedule);
		} catch (...){
			// do nothing
		}
		if (!spSpotSchedule){		
			spSpotSchedule.CoCreateInstance(CLSID_SpotSchedule);
			map_com_to_analytic(spSpotSchedule, SpotSchedule, h);
			hSpotSchedule = h;
			hSpotSchedule->PutName(m_h->GetName());
			hSpotSchedule->PutDataSource(m_h->GetDataSource());
			hSpotSchedule->PutDate(m_h->GetDate());
		} else {
			map_com_to_analytic(spSpotSchedule, SpotSchedule, h);
			hSpotSchedule = h;
		}		
		if (_Module.GetUseBwlForSpotSchedule()){
			// only merge in today's spot value (this emulates the SQL server case)
			(*hSpotSchedule)[m_h->GetDate()] = m_h->GetValueAt(m_h->GetDate());
		} else {
			// merge in all the values in m_h (note that m_h will take precedence)
			hSpotSchedule->Merge(*m_h);
		}
		if (hr = CXmlStreamer::GetXML(spSpotSchedule, ssXML)) return hr;
		return CComObjectCollectionFunctions<ISpotSchedule>::ImplementSave_NTFS(ssXML, IID_ISpotSchedule, _Module.GetLocation(), m_h->GetName(), m_h->GetDataSource(), 0/*we don't save spot schedules historically in the NTFS case*/, pVal);
	} else {
		return E_FAIL;
	}
	return E_FAIL;	
	end_function
}