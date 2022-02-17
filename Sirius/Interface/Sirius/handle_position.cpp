//	handle_position.cpp : Implementation of CPosition
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionfunctions.h"
#include "handle_position.h"
#include "siriusapplication.h"
#include "MlEqDate.h"
#include "comobjectcollectionserialisablekey.h"

implement_member_variable_rekey(Position, Positions, DATE, Date, long, Date);
implement_member_string_rekey(Position, Positions, Name, Name);
implement_member_variable_rekey(Position, Positions, DataSourceEnum, DataSource, DataSourceEnum, DataSource);
implement_serialisable(Position);


STDMETHODIMP CPosition::Evaluate(BSTR Calculate, IResult** pVal)
{
	begin_function
	return m_h->Evaluate(Calculate).CopyTo(pVal);
	end_function
}

STDMETHODIMP CPosition::GetFxUnderlyings(IAssets** pVal)
{	
	begin_function
	return m_h->GetProducts()->GetFxUnderlyings(pVal);	
	end_function
}

STDMETHODIMP CPosition::get_ProductNotional(IProduct* Product, double *pVal)
{
	begin_function
	if (!m_h->GetProducts()) return E_POINTER;
	return m_h->GetProducts()->get_Notional(CComVariant(Product), pVal);
	end_function
}

STDMETHODIMP CPosition::get_Products(IProducts** pVal)
{
	begin_function
	return m_h->GetProducts().CopyTo(pVal);
	end_function
}

STDMETHODIMP CPosition::get_Result(IResult** pVal)
{
	begin_function
	return m_h->GetResult().CopyTo(pVal);
	end_function
}

//	returns a collection of assets on which the position depends
STDMETHODIMP CPosition::GetUnderlyings(IAssets** pVal)
{	
	begin_function
	return m_h->GetProducts()->GetUnderlyings(pVal);
	end_function
}

STDMETHODIMP CPosition::get_Value(VARIANT *pVal)
{		
	return g_pApplication->GetObjectManager().ImplementGetValue(dynamic_cast<IPosition*>(this), pVal);
}

STDMETHODIMP CPosition::GetVolatilityStructures(IVolatilityStructures** pVal)
{
	begin_function
	return m_h->GetProducts()->GetVolatilityStructures(pVal);
	end_function
}

STDMETHODIMP CPosition::GetZeroCurves(IZeroCurves** pVal)
{	
	begin_function
	return m_h->GetProducts()->GetZeroCurves(pVal);
	end_function
}

STDMETHODIMP CPosition::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IPosition };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return E_FAIL;
}

/*static*/ HRESULT CPosition::Load(const std::string& szName, DataSourceEnum ds, long nDate, CComPtr<IPosition>& spPosition)
{		
	FileSystemEnum						fs;	
		
	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();

	// load the XML for the position	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;						// SQL query to execute
		_RecordsetPtr					prs;							// ADO recordset
		CComVariant						vField;							// position field value
						
		ssQuery << "sp_user_get_position '" << szName << "', " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "'";
		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
			if (prs->adoEOF) throw "Cannot find position '" + szName + "' on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
			vField = prs->GetFields()->GetItem(CComVariant(0, VT_I4))->GetValue();						
		} catch (_com_error& e){
			throw estring(e);
		}
		if (vField.vt != VT_BSTR) return E_FAIL;
		if (CXmlStreamer::GetObject((char*)_bstr_t(vField.bstrVal), ds, spPosition)) propagate_error;
	} else if (fs == fsNTFS){		
		CComObjectCollectionFunctions<IPosition>::ImplementLoad_NTFS(IID_IPosition, _Module.GetLocation(), "", szName, ds, nDate, spPosition);
	} else {
		return E_FAIL;
	}

	// Iterate through all the products in the position and set their data sources to ds.
	CComPtr<IProducts> spProducts;
	spPosition->get_Products(&spProducts);
	CComObjectCollectionFunctions<IProduct>().PutDataSource(spProducts, ds);
	return S_OK;
}

STDMETHODIMP CPosition::PutDataSource(/*[in]*/ DataSourceEnum newVal)
{
	begin_function
	return m_h->GetProducts()->PutDataSource(newVal);
	end_function
}

STDMETHODIMP CPosition::PutDate(DATE newVal)
{
	begin_function
	return m_h->GetProducts()->PutDate(newVal);
	end_function
}

STDMETHODIMP CPosition::put_Products(IProducts* newVal)
{
	begin_function
	m_h->PutProducts(newVal);
	end_function
}

STDMETHODIMP CPosition::put_ProductNotional(IProduct* Product, double newVal)
{
	begin_function
	if (!m_h->GetProducts()) return E_POINTER;
	return m_h->GetProducts()->put_Notional(CComVariant(Product), newVal);
	end_function
}

STDMETHODIMP CPosition::put_Result(IResult* newVal)
{
	begin_function
	m_h->PutResult(newVal);
	end_function
}

STDMETHODIMP CPosition::put_Value(VARIANT newVal)
{			
	return g_pApplication->GetObjectManager().ImplementPutValue(newVal, dynamic_cast<IPosition*>(this));
}

STDMETHODIMP CPosition::Save(BSTR* pVal)
{
	HRESULT								hr;	
	xmlstreamer							ssXML;							// XML representation of vData		
	FileSystemEnum						fs;
	
	begin_function
	check_publishing_enabled
	if (!m_h->GetDate()) throw "You cannot save the position '" + m_h->GetName() + "' since it does not have a date.";
	
	if (hr = CXmlStreamer::GetXML(CComPtr<IPosition>(this), ssXML)) return hr;
	
	// save it using the appropriate stored procedure	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){		
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;						// SQL query to execute	
		_RecordsetPtr					prs;							// ADO recordset
				
		if (m_h->GetDataSource() != Last) throw "Positions with data sources other than 'Last' cannot be saved. The position '" + m_h->GetName() + "' has the source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "'";
		ssQuery << "sp_user_add_position '" << m_h->GetName() << "', " << m_h->GetDate() << ", '" << m_h->GetBook() << "', '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, Last) << "', '" << _Module.GetLocation() << "', '" << (char*)ssXML << "'";
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
		// we relax the condition that only 'Last' data sources can be saved (or the UpdatePL process gets tricky)
		if (m_h->GetDataSource() == NoDataSource) throw "The position '" + m_h->GetName() + "' has the data source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "' and so cannot be saved";
		return CComObjectCollectionFunctions<IPosition>::ImplementSave_NTFS(ssXML, IID_IPosition, _Module.GetLocation(), m_h->GetName(), m_h->GetDataSource(), m_h->GetDate(), pVal);
	} else {
		return E_FAIL;
	}
	end_function
}
