//	handle_correlationmatrix.cpp : Implementation of CCorrelationMatrix
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionfunctions.h"
#include "handle_correlationmatrix.h"
#include "xmlstreamer.h"
#include "siriusapplication.h"
#include "MlEqDate.h"
#include "comobjectcollectionserialisablekey.h"

/*static*/ estring CCorrelationMatrix::s_szCorrelationMatrixName = "Global";

implement_member_variable_rekey(CorrelationMatrix, CorrelationMatrices, DATE, Date, long, Date);
implement_member_variable_rekey(CorrelationMatrix, CorrelationMatrices, DataSourceEnum, DataSource, DataSourceEnum, DataSource);
implement_serialisable(CorrelationMatrix);

STDMETHODIMP CCorrelationMatrix::Add(ICorrelationMatrix* CorrelationMatrix, VARIANT_BOOL Replace)
{
	CCorrelationMatrix*	pcm = dynamic_cast<CCorrelationMatrix*>(CorrelationMatrix);
	if (!pcm) return E_POINTER;
	begin_function
	m_h->Add(pcm->m_h, Replace ? true : false);
	end_function
}

STDMETHODIMP CCorrelationMatrix::Clear()
{
	begin_function
	m_h->Clear();
	end_function
}

STDMETHODIMP CCorrelationMatrix::CopyToClipboard()
{
	CParameterMap						pmValue;
	GetParameterMap(&pmValue);
	return pmValue.CopyToClipboard();	
}

HRESULT CCorrelationMatrix::FinalConstruct(void)
{
	m_h->PutName(s_szCorrelationMatrixName);
	return S_OK;
}

// This function makes a database call if necessary.
double CCorrelationMatrix::GetCorrelation(const std::string& sz1, const std::string& sz2)
{
	try {
		double f = m_h->GetCorrelation(sz1, sz2);
		return f;
	} catch (const std::string& szError){
		// Attempt to load a value from the database. We add this value to
		// this correlation matrix if successful.
		CComPtr<IDispatch>					spObject;
		CComQIPtr<ICorrelationMatrix>		spCorrelationMatrix;
		try {
			_Module.Load(CLSID_CorrelationMatrix, sz1 + "," + sz2, "", m_h->GetDate(), true, spObject);
		} catch (...){
			throw szError;
		}
		if (!(spCorrelationMatrix = spObject)) throw szError;
		map_com_to_analytic(spCorrelationMatrix, CorrelationMatrix, h);
		m_h->Add(h, false);
		return m_h->GetCorrelation(sz1, sz2);
	}
}

STDMETHODIMP CCorrelationMatrix::GetCorrelation(BSTR Name1, BSTR Name2, double *pVal)
{	
	begin_function
	*pVal = 0.0;
	*pVal = GetCorrelation(estring(Name1), estring(Name2));
	end_function	
}

STDMETHODIMP CCorrelationMatrix::get_Matrix(VARIANT* pVal)
{
	begin_function
	CParameterMap						pmValue;
	GetParameterMap(&pmValue);
	return pmValue.GetValue(pVal);
	end_function
}

STDMETHODIMP CCorrelationMatrix::get_Name(BSTR* pVal)
{
	begin_function
	return estring::GetBSTR(m_h->GetName(), pVal);
	end_function
}

void CCorrelationMatrix::GetParameterMap(CParameterMap* ppmValue)
{
	std::map<MlEqCorrelationMatrix::helper, double>::const_iterator		it;	
	std::vector<std::string>											vectorColumnHeadings;
	std::map<std::string, long>											mapColumnHeadings;
	std::vector<std::string>											vectorRowHeadings;
	std::map<std::string, long>											mapRowHeadings;
	CParameterMap														pmRowHeadings, pmColumnHeadings;
		
	ppmValue->SetSize(m_h->GetMap().size(), m_h->GetMap().size());		// maximum upper limit
	for (it = m_h->GetMap().begin(); it != m_h->GetMap().end(); it++)
	{		
		double f = it->second;
		
		if (mapRowHeadings[it->first.L()]){
			// Found lower part of the correlation in the row heading.
			if (!mapColumnHeadings[it->first.H()]){				
				// Add the higher part to the column heading
				vectorColumnHeadings.push_back(it->first.H());
				mapColumnHeadings[it->first.H()] = vectorColumnHeadings.size();
			}
			ppmValue->SetValue(mapRowHeadings[it->first.L()] - 1, mapColumnHeadings[it->first.H()] - 1, it->second);
		} else if (mapColumnHeadings[it->first.L()]){
			// Found lower part of the correlation in the column heading			
			if (!mapRowHeadings[it->first.H()]){
				// Add the higher part to the row heading
				vectorRowHeadings.push_back(it->first.H());
				mapRowHeadings[it->first.H()] = vectorRowHeadings.size();
			}
			ppmValue->SetValue(mapRowHeadings[it->first.H()] - 1, mapColumnHeadings[it->first.L()] - 1, it->second);			
		} else if (mapRowHeadings[it->first.H()]){
			// Found higher part of the correlation in the row heading
			if (!mapColumnHeadings[it->first.L()]){
				// Add the lower part to the column heading
				vectorColumnHeadings.push_back(it->first.L());
				mapColumnHeadings[it->first.L()] = vectorColumnHeadings.size();
			}
			ppmValue->SetValue(mapRowHeadings[it->first.H()], mapColumnHeadings[it->first.L()], it->second);
		} else if (mapColumnHeadings[it->first.H()]){
			// Found higher part of the correlation in the column heading
			if (!mapRowHeadings[it->first.L()]){				
				// Add the lower part to the row heading
				vectorRowHeadings.push_back(it->first.L());
				mapRowHeadings[it->first.L()] = vectorRowHeadings.size();
			}
			ppmValue->SetValue(mapRowHeadings[it->first.L()] - 1, mapColumnHeadings[it->first.H()] - 1, it->second);
		} else {
			// Nothing found. We add the lower part to the row headings and the higher part to the column headings.
			vectorRowHeadings.push_back(it->first.L());
			vectorColumnHeadings.push_back(it->first.H());
			mapRowHeadings[it->first.L()] = vectorRowHeadings.size();
			mapColumnHeadings[it->first.H()] = vectorColumnHeadings.size();
			ppmValue->SetValue(mapRowHeadings[it->first.L()] - 1, mapColumnHeadings[it->first.H()] - 1, it->second);
		}
	}
	ppmValue->SetRows(vectorRowHeadings.size());
	ppmValue->SetColumns(vectorColumnHeadings.size());

	// now go through pmValue and insert as many values in the blank elements as possible (using symmetry property)
	for (long nRow = 0; nRow < ppmValue->GetRows(); nRow++){
		for (long nCol = 0; nCol < ppmValue->GetCols(); nCol++){
			if (ppmValue->IsBlank(nRow, nCol)){				
				MlEqCorrelationMatrix::helper h(vectorRowHeadings[nRow], vectorColumnHeadings[nCol]);
				std::map<MlEqCorrelationMatrix::helper, double>::const_iterator		it = m_h->GetMap().find(h);
				if (it != m_h->GetMap().end()){
					// found a correlation value
					ppmValue->SetValue(nRow, nCol, it->second);
				}
			}
		}
	}

	// set pmColumnHeadings and pmRowHeadings		
	pmRowHeadings.SetValue(vectorRowHeadings);
	pmColumnHeadings.SetValue(vectorColumnHeadings);
	pmColumnHeadings.Transpose();

	// add these to pmValue
	pmRowHeadings.AddToRHS(*ppmValue);
	pmColumnHeadings.InsertColumn(0);
	pmColumnHeadings.AddToEnd(pmRowHeadings);
	ppmValue->Attach(&pmColumnHeadings);

	// Set the dimensions in the top left cell. It looks nicer in Excel if we have something in this cell, however meaningless!	
	std::stringstream ss;
	ss << "[" << ppmValue->GetRows() - 1<< " x " << ppmValue->GetCols() - 1 << "]";
	ppmValue->SetValue(0, 0, ss);
}

STDMETHODIMP CCorrelationMatrix::get_Value(VARIANT *pVal)
{
	begin_function
	std::string							szDataSource = CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource());
	CParameterMap						pmValue;
	
	unmap_parameter(szDataSource, pmDataSource);
	unmap_parameter(m_h->GetDate(), pmDate);
	GetParameterMap(&pmValue);
	return CParameterMap::VariableArgumentListToArray(pVal, 3, pmDate, pmDataSource, pmValue);
	end_function
}

STDMETHODIMP CCorrelationMatrix::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_ICorrelationMatrix };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

/*static*/ HRESULT CCorrelationMatrix::Load(const std::string& szIdentifier, DataSourceEnum ds, long nDate, CComPtr<ICorrelationMatrix>& spCorrelationMatrix)
{	
	FileSystemEnum						fs;	
			
	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();
	ATLASSERT(!spCorrelationMatrix);
		
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;						// SQL query to execute
		_RecordsetPtr					prs;							// ADO recordset		
				
		ssQuery << "sp_user_get_correlation " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "', '" << szIdentifier << "'";
		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
		} catch (_com_error& e){
			throw estring(e);
		}		
		while (!prs->adoEOF){						
			CComVariant					vItem1;
			CComVariant					vItem2;
			CComVariant					vCorrelation;
			if (!spCorrelationMatrix){
				// Create the new matrix
				spCorrelationMatrix.CoCreateInstance(CLSID_CorrelationMatrix);
				spCorrelationMatrix->put_Date(nDate);
				spCorrelationMatrix->put_DataSource(ds);
			}
			vItem1 = prs->GetFields()->GetItem(0L)->GetValue();
			vItem2 = prs->GetFields()->GetItem(1L)->GetValue();
			vCorrelation = prs->GetFields()->GetItem(2L)->GetValue();
			ATLASSERT(vItem1.vt == VT_BSTR && vItem2.vt == VT_BSTR && vCorrelation.vt == VT_R8);
			spCorrelationMatrix->SetCorrelation(vItem1.bstrVal, vItem2.bstrVal, vCorrelation.dblVal);
			prs->MoveNext();
		}
		
		if (!spCorrelationMatrix){
			if (szIdentifier.size()){
				throw "Cannot find any correlation values for '" + szIdentifier + "' on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
			} else {
				throw "Cannot find any correlation values on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
			}
		}
	} else if (fs == fsNTFS){
		CComObjectCollectionFunctions<ICorrelationMatrix>::ImplementLoad_NTFS(IID_ICorrelationMatrix, _Module.GetLocation(), "", s_szCorrelationMatrixName, ds, nDate, spCorrelationMatrix);
		return S_OK;
	} else {
		return E_FAIL;
	}
	return S_OK;
}

MlEqCorrelationMatrixHandle CCorrelationMatrix::PutMap(CParameterMap* ppm, const std::string& szName, DataSourceEnum ds, long nDate) const
//	This function corrupts ppm
{
	CParameterMap						pmRowHeadings, pmColumnHeadings;
	MlEqCorrelationMatrixHandle			h(new MlEqCorrelationMatrix);

	// set map
	if (ppm->IsBlank()) throw CStringResource(IDS_NO_DATA_FOUND);
	ppm->SetBlank(0, 0);
	if (ppm->GetColumn(0, &pmRowHeadings)) propagate_error;
	if (ppm->GetRow(0, &pmColumnHeadings)) propagate_error;
	if (ppm->RemoveRow(0)) propagate_error;
	if (ppm->RemoveColumn(0)) propagate_error;
	if (!ppm->IsBlankOrDouble(-1, -1, IDS_INVALID_PARAMETER)) throw "Unhandled exception in CCorrelationMatrix::PutMap";
	std::string szRow, szColumn;
	
	// Process the correlations from the bottom left to top right.
	// In this way, we support upper triangular form for symmetric matrices.
	for (long nRow = pmRowHeadings.GetRows() - 1; nRow >= 0 ; nRow--){
		if (pmRowHeadings.GetValue(nRow, 0, &szRow)) propagate_error;
		for (long nCol = 0; nCol < pmColumnHeadings.GetCols(); nCol++){
			if (pmColumnHeadings.GetValue(0, nCol, &szColumn)) propagate_error;
			double fCorrelation;
			// We allow duplicate correlations (we simply overwrite any existing data).
			if (ppm->GetValue(nRow, nCol, &fCorrelation)){
				std::string sz;
				ppm->GetValue(nRow, nCol, &sz);
				throw "Invalid parameter '" + sz + "'";
			}
			h->SetCorrelation(szRow, szColumn, fCorrelation);
		}
	}
	
	h->PutDataSource(ds);
	h->PutDate(nDate);
	h->PutName(szName);
	return h;
}

STDMETHODIMP CCorrelationMatrix::put_Matrix(VARIANT newVal)
{
	begin_function
	HRESULT								hr;
	CParameterMap						pm;
	
	if (hr = pm.SetValue(newVal)) return hr;
	MlEqCorrelationMatrixHandle h = PutMap(&pm, m_h->GetName(), m_h->GetDataSource(), m_h->GetDate());
	m_h = h;
	end_function
}

STDMETHODIMP CCorrelationMatrix::put_Value(VARIANT newVal)
{
	/*parameter list is 0 - DateOpt
						1 - DataSource
						2 - Matrix */
	
	begin_function
	HRESULT													hr;
	std::vector<CParameterMap>								vpm;
	
	if (hr = CParameterMap::ArrayToVector(newVal, NULL, &vpm)) return hr;
	if (vpm.size() != 3) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_ICorrelationMatrix);
	map_optional_parameter(vpm[0], long, DateOpt, MlEqDate::GetCurrentDate())
	map_enum_parameter(vpm[1], DataSourceEnum, DataSourceOpt)
	MlEqCorrelationMatrixHandle	h = PutMap(&vpm[2], s_szCorrelationMatrixName, DataSourceOpt, DateOpt);
	
	// All is well if this point is reached - set the member variables
	m_h = h;
	
	// Merge this correlation matrix into the corresponding matrix in
	// the object manager's maintained collection.
	g_pApplication->GetCorrelationMatrixHandle(m_h->GetDataSource(), m_h->GetDate(), true)->Add(m_h, false);
	end_function
}

STDMETHODIMP CCorrelationMatrix::Save(BSTR* pVal)
{
	HRESULT								hr;	
	xmlstreamer							ssXML;							// XML representation of vData			
	FileSystemEnum						fs;
	
	begin_function
	check_publishing_enabled
	if (!m_h->GetDate()) throw "You cannot save the correlation matrix '" + m_h->GetName() + "' since it does not have a date.";
			
	// save it using the appropriate stored procedure
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		_ConnectionPtr					spConnection;
		if (m_h->GetDataSource() != Last) throw "Correlation Matrices with data sources other than 'Last' cannot be saved. This correlation matrix has the source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "'";
		for (std::map<MlEqCorrelationMatrix::helper, double>::const_iterator it = m_h->GetMap().begin(); it != m_h->GetMap().end(); it++){
			std::stringstream				ssQuery;						// SQL query to execute
			ssQuery << "sp_user_add_correlation " << m_h->GetDate() << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) << "' ,'" << _Module.GetLocation() << "', '" << it->first.L() << "', '" << it->first.H() << "', " << it->second;
			try {
				_Module.GetConnection(spConnection);
				spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
			} catch (_com_error& e){
				throw estring(e);
			}
		}
		return estring(s_szCorrelationMatrixName).GetBSTR(pVal);
	} else if (fs == fsNTFS){
		// We relax the condition that only 'Last' data sources can be saved (or the UpdatePL process gets tricky)
		// To emulate the SQL server case, we need to load the current saved matrix, merge in this one then resave.
		if (m_h->GetDataSource() == NoDataSource) throw "The correlation matrix '" + m_h->GetName() + "' has the data source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "' and so cannot be saved";
		CComPtr<ICorrelationMatrix>	spCorrelationMatrix;
		try {
			CCorrelationMatrix::Load("", m_h->GetDataSource(), m_h->GetDate(), spCorrelationMatrix);
		} catch (...){
			// Do nothing
		}
		if (spCorrelationMatrix){
			spCorrelationMatrix->Add(this, VARIANT_FALSE);
			if (hr = CXmlStreamer::GetXML(spCorrelationMatrix, ssXML)) return hr;
		} else {
			if (hr = CXmlStreamer::GetXML(CComPtr<ICorrelationMatrix>(this), ssXML)) return hr;
		}
		return CComObjectCollectionFunctions<ICorrelationMatrix>::ImplementSave_NTFS(ssXML, IID_ICorrelationMatrix, _Module.GetLocation(), s_szCorrelationMatrixName, m_h->GetDataSource(), m_h->GetDate(), pVal);
	} else {
		return E_FAIL;
	}
	end_function
}

STDMETHODIMP CCorrelationMatrix::SetCorrelation(BSTR Name1, BSTR Name2, double fCorrelation)
{
	m_h->SetCorrelation(estring(Name1), estring(Name2), fCorrelation);
	return S_OK;
}