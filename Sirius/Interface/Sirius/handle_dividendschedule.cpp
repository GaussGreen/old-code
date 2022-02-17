//	handle_dividendschedule.cpp : Implementation of CDividendSchedule
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionfunctions.h"
#include "handle_dividendschedule.h"
#include "mleqzerocurve.h"
#include "siriusapplication.h"
#include "xmlstreamer.h"
#include "comobjectcollectionserialisablekey.h"

implement_member_variable_rekey(DividendSchedule, DividendSchedules, DATE, Date, long, Date);
implement_member_variable_rekey(DividendSchedule, DividendSchedules, DataSourceEnum, DataSource, DataSourceEnum, DataSource);
implement_serialisable(DividendSchedule);


STDMETHODIMP CDividendSchedule::CalibrateYield(double Spot, IZeroCurve* ZeroCurve, BSTR Tenor, BSTR AddAt)
{
	begin_function	
	map_bare_com_to_analytic(ZeroCurve, ZeroCurve, hzc);
	m_h->CalibrateYield(Spot, hzc, estring(Tenor), estring(AddAt));
	end_function
}

STDMETHODIMP CDividendSchedule::DiscreteShift(double Growth, BSTR Tenor)
{
	begin_function
	m_h->DiscreteShift(Growth, estring(Tenor));
	end_function
}

STDMETHODIMP CDividendSchedule::get_DayCountConvention(DayCountConventionEnum *pVal)
{
	begin_function
	*pVal = m_h->GetDayCountConvention();
	end_function
}

STDMETHODIMP CDividendSchedule::get_FirstDate(DividendTypeEnum DividendType, DATE* pVal)
{
	begin_function
	switch (DividendType){
	case Discrete:
		*pVal = m_h->Discrete.GetFirstDate();
		break;
	case Continuous:
		*pVal = m_h->Continuous.GetFirstDate();
		break;
	default:
		throw "Unknown or unsupported dividend type";
	}
	end_function
}

STDMETHODIMP CDividendSchedule::get_LastDate(DividendTypeEnum DividendType, DATE* pVal)
{
	begin_function
	switch (DividendType){
	case Discrete:
		*pVal = m_h->Discrete.GetLastDate();
		break;
	case Continuous:
		*pVal = m_h->Continuous.GetLastDate();
		break;
	default:
		throw "Unknown or unsupported dividend type";
	}
	end_function
}

STDMETHODIMP CDividendSchedule::get_Name(BSTR* pVal)
{
	begin_function
	estring sz = m_h->GetName();
	return sz.GetBSTR(pVal);
	end_function
}
	
STDMETHODIMP CDividendSchedule::GetPresentValue(IZeroCurve* ZeroCurve, DATE From, DATE To, double *pVal)
{		
	begin_function	
	map_bare_com_to_analytic(ZeroCurve, ZeroCurve, hZeroCurve)		
	*pVal = m_h->GetPresentValue(hZeroCurve, From, To);
	end_function	
}

STDMETHODIMP CDividendSchedule::get_Schedule(VARIANT* pVal)
{	
	HRESULT								hr;
	CParameterMap						pmValue;

	begin_function
	if (hr = pmValue.SetValue(m_h->Discrete, m_h->Continuous)) return hr;
	return pmValue.GetValue(pVal);
	end_function
}

STDMETHODIMP CDividendSchedule::get_Value(VARIANT* pVal)
{
	HRESULT								hr;
	CParameterMap						pmValue;
	CParameterMap						pmYieldCurve;
	std::string							szDataSource = CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource());
	
	begin_function
	// name
	unmap_parameter(m_h->GetName(), pmName);
	// date
	unmap_parameter(m_h->GetDate(), pmDate);
	// data source
	unmap_parameter(szDataSource, pmDataSource);	
	// yield curve
	if (m_spZeroCurve){
		pmYieldCurve.SetValue(m_spZeroCurve.p);
	}
	// day count convention	
	unmap_parameter(CEnumMap::GetString("DayCountConventionEnum", LIBID_Sirius, m_h->GetDayCountConvention()), pmDcc);
	// schedule data
	if (hr = pmValue.SetValue(m_h->Discrete, m_h->Continuous)) return hr;
	if (!pmValue.IsBlank()){
		return CParameterMap::VariableArgumentListToArray(pVal, 6, pmName, pmDate, pmDataSource, pmDcc, pmYieldCurve, pmValue);
	} else {
		return CParameterMap::VariableArgumentListToArray(pVal, 5, pmName, pmDate, pmDataSource, pmDcc, pmYieldCurve);
	}
	end_function
}

STDMETHODIMP CDividendSchedule::get_ValueAt(DividendTypeEnum DividendType, DATE Date, double* pVal)
{
	begin_function
	switch (DividendType){
	case Discrete:
		*pVal = m_h->Discrete.GetValueAt(Date);
		break;
	case Continuous:
		*pVal = m_h->Continuous.GetValueAt(Date);
		break;
	default:
		throw "Unknown or unsupported dividend type";
	}
	end_function
}

STDMETHODIMP CDividendSchedule::GetYield(DATE From, DATE To, double *Yield)
{
	begin_function
	*Yield = m_h->GetYield(From, To);
	end_function
}

STDMETHODIMP CDividendSchedule::get_YieldCurve(IZeroCurve** pVal)
{
	begin_function
	return m_spZeroCurve.CopyTo(pVal);
	end_function
}

STDMETHODIMP CDividendSchedule::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IDividendSchedule };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

/*static*/ HRESULT CDividendSchedule::Load(const std::string& szIdentifier, DataSourceEnum ds, long nDate, CComPtr<IDividendSchedule>& spDividendSchedule)
{	
	FileSystemEnum						fs;
			
	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();

	// load the XML for the dividend schedule
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;						// SQL query to execute
		_RecordsetPtr					prs;							// ADO recordset
		CComVariant						vField;							// dividend schedule field value
						
		ssQuery << "sp_user_get_dividendschedule '" << szIdentifier << "', " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "'";
		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
			if (prs->adoEOF) throw "Cannot find dividend schedule '" + szIdentifier + "' on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
			vField = prs->GetFields()->GetItem(CComVariant(0, VT_I4))->GetValue();
		} catch (_com_error& e){
			throw estring(e);
		}
		if (vField.vt != VT_BSTR) return E_FAIL;				
		return CXmlStreamer::GetObject((char*)_bstr_t(vField.bstrVal), ds, spDividendSchedule);
	} else if (fs == fsNTFS){		
		CComObjectCollectionFunctions<IDividendSchedule>::ImplementLoad_NTFS(IID_IDividendSchedule, _Module.GetLocation(), "", szIdentifier, ds, nDate, spDividendSchedule);
		return S_OK;
	} else {
		return E_FAIL;
	}	
}

STDMETHODIMP CDividendSchedule::put_DayCountConvention(DayCountConventionEnum newVal)
{
	begin_function
	m_h->PutDayCountConvention(newVal);	
	end_function
}

STDMETHODIMP CDividendSchedule::put_Schedule(VARIANT newVal)
{
	HRESULT								hr;
	CParameterMap						pm;
	std::map<long, double>				discrete, continuous;
	
	begin_function
	if (hr = pm.SetValue(newVal)) return hr;
	if (hr = pm.GetValue(0, 1, 2, discrete, continuous)) return hr;
	m_h->PutDiscrete(discrete);
	m_h->PutContinuous(continuous);
	end_function
}

STDMETHODIMP CDividendSchedule::put_Value(VARIANT newVal)
{
	std::vector<CParameterMap>			vpm;
	HRESULT								hr;
	std::map<long, double>				discrete, continuous;

	/*parameter list is 0 - name
						1 - DateOpt
						2 - data source
						3 - day count convention
						4 - yield curve handle (understood as providing additional continuous yield data)
						5 - array of dates (understood as Ex-Dates)
						6 - value of discrete dividend paid at a given date
						7 - value of continuous yield between the given date and the previous date*/
	
	begin_function
	// resize the vpm array, dealing with optional parameters and cases where
	// the discrete / continuous values are concatenated to the dates
	if (hr = CParameterMap::ArrayToVector(newVal, NULL, &vpm)) return hr;
	if (!vpm.size()){
		return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IDividendSchedule);
	} else if (vpm.size() < 6){
		// contains no schedule data - this is OK
		vpm.resize(5);
	} else if (vpm.size() == 6){
		// vpm[5] must be either 2 or 3 columns (or zero columns!)
		if (!vpm[5].IsBlank() && vpm[5].GetCols() != 2 && vpm[5].GetCols() != 3) return CParameterMap::ReturnErrorR(IDS_COLUMNS_INVALID, IID_IDividendSchedule);
		if (!vpm[5].IsBlankOrDouble(-1, -1, IDS_INVALID_PARAMETER_VALUE)) return E_FAIL;
	} else if (vpm.size() == 7 || vpm.size() == 8){
		if (!vpm[5].IsBlankOrDouble(-1, -1, IDS_INVALID_PARAMETER_VALUE)) return E_FAIL;
		if (vpm.size() == 7 && !vpm[6].IsBlankOrDouble(-1, -1, IDS_INVALID_PARAMETER_VALUE)) return E_FAIL;
		if (vpm.size() == 8 && vpm[6].IsBlank()){
			// Nasty horsey case where the discrete dividend array is empty but the contiuous schedule is not.			
			vpm[6].SetSize(vpm[7].GetRows(), 1);
			vpm[6].SetValue(-1, -1, 0.0);
		}
		for (long n = 6; n < vpm.size(); n++){
			vpm[5].AddToRHS(vpm[n]);
		}
		vpm.resize(6);
	} else {
		return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IDividendSchedule);
	}
	
	// get temporary assignments
	map_parameter(vpm[0], estring, Name);
	map_optional_parameter(vpm[1], long, DateOpt, MlEqDate::GetCurrentDate())
	map_enum_parameter(vpm[2], DataSourceEnum, DataSource);
	map_optional_enum_parameter(vpm[3], DayCountConventionEnum, DayCountConventionOpt, ActualActual);
	map_optional_com_object_parameter(vpm[4], ZeroCurve, YieldCurveOpt);
	map_com_to_analytic(YieldCurveOpt, ZeroCurve, hZeroCurve);				
	if (vpm.size() > 5 && (hr = vpm[5].GetValue(0, 1, 2, discrete, continuous))) return hr;
	MlEqDividendScheduleHandle h = new MlEqDividendSchedule(Name, DataSource, DateOpt, DayCountConventionOpt, hZeroCurve, discrete, continuous);
		
	// All is well if this point is reached - set the member variables.
	m_h = h;
	m_spZeroCurve = YieldCurveOpt;
	end_function
}

STDMETHODIMP CDividendSchedule::put_ValueAt(DividendTypeEnum DividendType, DATE Date, double newVal)
{
	begin_function
	switch (DividendType){
	case Discrete:
		m_h->PutDiscreteValueAt(Date, newVal);
		break;
	case Continuous:
		m_h->PutContinuousValueAt(Date, newVal);
		break;
	default:
		throw "Unknown or unsupported dividend type";
	}
	end_function
}

STDMETHODIMP CDividendSchedule::put_YieldCurve(IZeroCurve* newVal)
{
	begin_function	
	m_spZeroCurve = newVal;
	map_com_to_analytic(m_spZeroCurve, ZeroCurve, hYieldCurve);
	m_h->PutForeignCurve(hYieldCurve);
	end_function
}

STDMETHODIMP CDividendSchedule::Refresh(VARIANT_BOOL Recursive)
{
	begin_function
	g_pApplication->GetObjectManager().Refresh(this, Recursive ? true : false);
	end_function
}

STDMETHODIMP CDividendSchedule::Reset()
{
	begin_function
	m_h->Reset();
	end_function		
}

STDMETHODIMP CDividendSchedule::Save(BSTR *pVal)
{
	HRESULT								hr;	
	xmlstreamer							ssXML;							// XML representation of vData	
	FileSystemEnum						fs;
	
	begin_function
	check_publishing_enabled
	if (!m_h->GetDate()) throw "You cannot save the dividend schedule '" + m_h->GetName() + "' since it does not have a date.";

	hr = CXmlStreamer::GetXML(CComPtr<IDividendSchedule>(this), ssXML);
	if (hr) return hr;

	// save it using the appropriate stored procedure	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		std::stringstream				ssQuery;						// SQL query to execute	
		_RecordsetPtr					prs;							// ADO recordset
		_ConnectionPtr					spConnection;
		
		if (m_h->GetDataSource() != Last) throw "Dividend Schedules with data sources other than 'Last' cannot be saved. The dividend schedule '" + m_h->GetName() + "' has the source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "'";
		ssQuery << "sp_user_add_dividendschedule '" << m_h->GetName() << "', '" << m_h->GetDate() << "', '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) << "', '" << _Module.GetLocation() << "', '" << (char*)ssXML << "'";
		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
		} catch (_com_error& e){
			throw estring(e);
		}
		// set the return type the returned value of the stored procedure
		CComVariant vValue(prs->GetFields()->GetItem(CComVariant(0, VT_I4))->GetValue());
		if (vValue.vt != VT_BSTR){
			ATLASSERT(false);
			return CParameterMap::ReturnErrorR(IDS_UNHANDLED_EXCEPTION);
		}		
		return CComBSTR(vValue.bstrVal).CopyTo(pVal);		
	} else if (fs == fsNTFS){		
		// We relax the condition that only 'Last' data sources can be saved (or the UpdatePL process gets tricky)
		if (m_h->GetDataSource() == NoDataSource) throw "The dividend schedule '" + m_h->GetName() + "' has the data source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "' and so cannot be saved";
		return CComObjectCollectionFunctions<IDividendSchedule>::ImplementSave_NTFS(ssXML, IID_IDividendSchedule, _Module.GetLocation(), m_h->GetName(), m_h->GetDataSource(), m_h->GetDate(), pVal);
	}
	end_function
}

STDMETHODIMP CDividendSchedule::Shift(double Yield)
{
	begin_function
	m_h->Shift(Yield);
	end_function	
}

STDMETHODIMP CDividendSchedule::Stick(void)
{
	begin_function
	m_h->Stick();
	end_function
}