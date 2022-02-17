//	handle_zerocurve.cpp : Implementation of CZeroCurve
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionfunctions.h"
#include "handle_zerocurve.h"
#include "siriusapplication.h"
#include "xmlstreamer.h"
#include "handle_interpolator.h"
#include "MlEqObjects.h"
#include "ComObjectCollectionSerialisableKey.h"

implement_member_variable_rekey(ZeroCurve, ZeroCurves, DATE, Date, long, ReferenceDate);
implement_member_variable_rekey(ZeroCurve, ZeroCurves, DataSourceEnum, DataSource, DataSourceEnum, DataSource);
implement_serialisable(ZeroCurve);

STDMETHODIMP CZeroCurve::Default(BSTR Currency, BSTR Request, BSTR* pVal)
{
	begin_function
	map_parameter(Currency, estring, szCurrency);
	map_parameter(Request, estring, szRequest);
	estring sz = MlEqZeroCurve::Default(szCurrency, szRequest);
	return sz.GetBSTR(pVal);
	end_function	
}

HRESULT CZeroCurve::FinalConstruct(void)
{
	return S_OK;
}

STDMETHODIMP CZeroCurve::GetCashRate(BSTR TermOrDate, double* pVal)
{
	begin_function
	*pVal = m_h->GetCashRate(estring(TermOrDate));
	end_function
}

STDMETHODIMP CZeroCurve::GetDiscountFactor(DATE Maturity, double* pVal)
{			
	begin_function
	*pVal = m_h->GetDiscountFactor(Maturity);
	end_function	
}

STDMETHODIMP CZeroCurve::GetForwardDiscountFactor(DATE From, DATE To, double* pVal)
{			
	begin_function
	*pVal = m_h->GetDiscountFactor(From, To);
	end_function	
}

std::string CZeroCurve::GetGDAHandle(void) const
{
	throw "CZeroCurve::GetGDAHandle is not implemented in this release";
}

STDMETHODIMP CZeroCurve::get_FxSpot(double *pVal)
{
	begin_function
	*pVal = m_h->GetFxSpot();
	end_function	
}

STDMETHODIMP CZeroCurve::get_InterpolatorType(BSTR* pVal)
{
	begin_function	
	return estring(m_h->GetInterpolator()).GetBSTR(pVal);
	end_function
}

STDMETHODIMP CZeroCurve::get_Name(BSTR* pVal)
{
	begin_function
	return estring(m_h->GetName()).GetBSTR(pVal);
	end_function
}

STDMETHODIMP CZeroCurve::get_Value(VARIANT* pVal)
{				
	begin_function	
	std::stringstream					ss;
	GDA::XMLOutStream					xml(ss);
	estring								sz;
	std::vector<std::string>			vector;
	std::string							szDataSource = CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource());
	
	xml << m_h->GetYieldCurve();
	sz.Set(ss);
	sz.Split("\n", &vector);
	unmap_parameter(vector, pmGdaData);
	unmap_parameter(szDataSource, pmDataSource);
	return CParameterMap::VariableArgumentListToArray(pVal, 2, pmGdaData, pmDataSource);				
	end_function
}

STDMETHODIMP CZeroCurve::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IZeroCurve };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

/*static*/ HRESULT CZeroCurve::Load(const std::string& szIdentifier, DataSourceEnum ds, long nDate, CComPtr<IZeroCurve>& spZeroCurve)
{	
	FileSystemEnum						fs;
			
	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();
	
	// load the XML for the zero curve	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;						// SQL query to execute
		_RecordsetPtr					prs;							// ADO recordset
		CComVariant						vField;							// zero curve field value

		ssQuery << "sp_user_get_zerocurve '" << szIdentifier << "', " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "'";
		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
			if (prs->adoEOF) throw "Cannot find zero curve '" + szIdentifier + "' on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
			vField = prs->GetFields()->GetItem(CComVariant(0, VT_I4))->GetValue();
		} catch (_com_error& e){
			throw estring(e);
		}
		if (vField.vt != VT_BSTR) return E_FAIL;
		
		// parse the XML into a variant and set this zero curve object to that variant
		return CXmlStreamer::GetObject((char*)_bstr_t(vField.bstrVal), ds, spZeroCurve);
	} else if (fs == fsNTFS){
		CComObjectCollectionFunctions<IZeroCurve>::ImplementLoad_NTFS(IID_IZeroCurve, _Module.GetLocation(), "", szIdentifier, ds, nDate, spZeroCurve);
		return S_OK;
	} else {
		return E_FAIL;
	}	
}

STDMETHODIMP CZeroCurve::put_FxSpot(double newVal)
{
	begin_function
	m_h->PutFxSpot(newVal);
	end_function
}

STDMETHODIMP CZeroCurve::put_InterpolatorType(BSTR newVal)
{
	begin_function
	estring								szInterpolatorType(newVal);
	m_h->PutInterpolator(szInterpolatorType);
	end_function
}

STDMETHODIMP CZeroCurve::put_Value(VARIANT newVal)
{	
	begin_function
			
	HRESULT											hr;
	std::vector<CParameterMap>						vpm;
	GDA::Functor_const_ref							hYieldCurve;		// either created or abstracted from a Gda handle
				
	if (hr = CParameterMap::ArrayToVector(newVal, NULL, &vpm)) return hr;		
	if (vpm.size() != 2) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IZeroCurve);

	/*parameter list is 0 - GdaHandle / GdaData
						1 - DataSource*/
			
	map_parameter(vpm[0], std::vector<std::string>, GdaHandle)
	map_enum_parameter(vpm[1], DataSourceEnum, DataSource)

	if (!GdaHandle.size()) return CParameterMap::ReturnErrorR(IDS_NO_DATA_FOUND, IID_IZeroCurve);
	estring::trim(&GdaHandle[0]);
			
	if (!GdaHandle[0].size()){
		// Probably because there is no encoding string.
		GdaHandle.erase(GdaHandle.begin());
		if (!GdaHandle.size()) return CParameterMap::ReturnErrorR(IDS_NO_DATA_FOUND, IID_IZeroCurve);
		estring::trim(&GdaHandle[0]);
	}
		
	if (!estring::CompareNoCase("<?xml", estring::left(&GdaHandle[0], 5)) ||
		!estring::CompareNoCase("<HDElement", estring::left(&GdaHandle[0], 10))){
		// this is GdaData
		std::stringstream					ss;
		GDA::XMLInStream					xml(ss);
		GDA::HDElement						hdeXml;						// this is the HDElement streamed into hdeXml
		for (long n = 0; n < GdaHandle.size(); n++){
			ss << GdaHandle[n] << '\n';		// we don't care about the extra '\n' that this will append!
		}
		xml >> hdeXml;		
		if (hdeXml.isSerialisation()){
			hdeXml = hdeXml.lookup("State").lookup("Parameters");
		}
		hYieldCurve = GDA::GetLibraryInstance()->getFactory()->create("Curve.YieldCurve", hdeXml, GDA::GetDefaultContext());
	} else if (GdaHandle.size() == 1){
		// assume a Gda handle has been passed in		
		try {
			GDA::Context_ref					ctx = GDA::GetLibraryInstance()->getDefaultNamespace()->getContext();
			GDA::HDElement						hdeYieldCurve = ctx->resolve(GdaHandle[0].data());
			hYieldCurve = GDA::Functor_const_ref(hdeYieldCurve.asConstObject());
		} catch (...){
			return CParameterMap::ReturnErrorRS(IDS_NO_OBJECT_WITH_HANDLE, GdaHandle[0], IID_IZeroCurve);
		}
	} else {
		return CParameterMap::ReturnErrorR(IDS_COULD_NOT_PARSE_XML, IID_IZeroCurve);		
	}
	
	m_h->PutYieldCurve(hYieldCurve);
	m_h->PutDataSource(DataSource);
	end_function
}

STDMETHODIMP CZeroCurve::Reset()
{	
	begin_function
	m_h->Reset();
	end_function	
}

STDMETHODIMP CZeroCurve::Save(BSTR *pVal)
//	pVal - returned, nullable
{
	HRESULT								hr;	
	xmlstreamer							ssXML;							// XML representation of vData		
	FileSystemEnum						fs;
	
	begin_function
	check_publishing_enabled
	if (!m_h->GetReferenceDate()) throw "You cannot save the zero curve '" + m_h->GetName() + "' since it does not have a date.";
			
	hr = CXmlStreamer::GetXML(CComPtr<IZeroCurve>(this), ssXML);	
	if (hr) return hr;

	// save it using the appropriate stored procedure	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;						// SQL query to execute	
		_RecordsetPtr					prs;							// ADO recordset
		
		if (m_h->GetDataSource() != Last) throw "Zero curves with data sources other than 'Last' cannot be saved. The zero curve '" + m_h->GetName() + "' has the source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "'";
		ssQuery << "sp_user_add_zerocurve '" << m_h->GetName() << "', '" << m_h->GetReferenceDate() << "', '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) << "', '" << _Module.GetLocation() << "', '" << (char*)ssXML << "'";
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
		// we relax the condition that only 'Last' data sources can be saved (or the UpdatePL process gets tricky)
		return CComObjectCollectionFunctions<IZeroCurve>::ImplementSave_NTFS(ssXML, IID_IZeroCurve, _Module.GetLocation(), m_h->GetName(), m_h->GetDataSource(), m_h->GetReferenceDate(), pVal);
	}
	end_function
}

STDMETHODIMP CZeroCurve::Shift(double Amount)	
{
	begin_function		
	m_h->PutShift(Amount);	
	end_function
}

STDMETHODIMP CZeroCurve::ShiftByInterpolator(IInterpolator* Interpolator)
{
	begin_function	
	map_bare_com_to_analytic(Interpolator, Interpolator, h);		
	m_h->PutShift(h);
	end_function
}

STDMETHODIMP CZeroCurve::Stick(void)
{
	begin_function
	m_h->Stick();
	end_function
}