//	SiriusApplication.cpp : Implementation of CSiriusApplication
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionfunctions.h"
#include "siriusapplication.h"
#include "handle_product.h"
#include "handle_position.h"
#include "handle_deal.h"
#include "handle_deals.h"
#include "handle_positions.h"
#include "handle_zerocurve.h"
#include "handle_zerocurves.h"
#include "handle_spotschedule.h"
#include "handle_spotschedules.h"
#include "handle_dividendschedule.h"
#include "handle_dividendschedules.h"
#include "handle_volatilitystructure.h"
#include "handle_volatilitystructures.h"
#include "handle_correlationmatrix.h"
#include "handle_correlationmatrices.h"
#include "handle_asset.h"
#include "handle_assets.h"
#include "comobjectcollectionserialisablekey.h"


////////////////////////////////////////////////////////////////////////////
//	standard installation macros
//
BEGIN_IMPLEMENT_COLLECTION()		// this macro block implements HRESULT CSiriusApplication::InstallStandardObjects(void)	
	IMPLEMENT_COLLECTION(LIBID_Sirius, ARPropVolDataCollection, ARPropVolData, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, ARPropFitVolDataCollection, ARPropFitVolData, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, Arrays, Array, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, Dates, Date, true)
	IMPLEMENT_COLLECTION(LIBID_Sirius, DateSchedules, DateSchedule, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, Deals, Deal, true)	
	IMPLEMENT_COLLECTION(LIBID_Sirius, HermiteVolDataCollection, HermiteVolData, false)	
	IMPLEMENT_COLLECTION(LIBID_Sirius, HullAndWhiteCollection, HullAndWhite, false)	
	IMPLEMENT_COLLECTION(LIBID_Sirius, HullVolDataCollection, HullVolData, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, Interpolators, Interpolator, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, JumpWingVolDataCollection, JumpWingVolData, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, JumpWingFitVolDataCollection, JumpWingFitVolData, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, Matrices, Matrix, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, MonteCarloCollection, MonteCarlo, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, MultiplierVolDataCollection, MultiplierVolData, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, Parameters, Parameter, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, ParameterLists, ParameterList, true)
	IMPLEMENT_COLLECTION(LIBID_Sirius, Positions, Position, true)
	IMPLEMENT_COLLECTION(LIBID_Sirius, Products, Product, true)	
	IMPLEMENT_COLLECTION(LIBID_Sirius, RamFitVolDataCollection, RamFitVolData, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, RamVolDataCollection, RamVolData, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, Results, Result, false)		
	IMPLEMENT_COLLECTION(LIBID_Sirius, Scenarios, Scenario, true)
	IMPLEMENT_COLLECTION(LIBID_Sirius, StrikesCollection, Strikes, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, SVLFitVolDataCollection, SVLFitVolData, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, SVLVolDataCollection, SVLVolData, false)
	IMPLEMENT_COLLECTION(LIBID_Sirius, VolDataCollection, VolData, false)	
	IMPLEMENT_COLLECTION(LIBID_Sirius, PdeCollection, Pde, false)
END_IMPLEMENT_COLLECTION()


////////////////////////////////////////////////////////////////////////////
//	Clear
//
//	Delete all the handles in the object manager.
//
STDMETHODIMP CSiriusApplication::Clear(void)
{
	begin_function
	m_ObjectManager.ClearCollections();
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	CreateObject
//
//	Creates an object of a given name with input arguments. This function
//  provides a VB analogue of MLCreate.
//
STDMETHODIMP CSiriusApplication::CreateObject(BSTR ObjectName, VARIANT Arg1, VARIANT Arg2, VARIANT Arg3, VARIANT Arg4, VARIANT Arg5, VARIANT Arg6, VARIANT Arg7, VARIANT Arg8, VARIANT Arg9, VARIANT Arg10, VARIANT Arg11, VARIANT Arg12, VARIANT Arg13, VARIANT Arg14, VARIANT Arg15, VARIANT Arg16, VARIANT Arg17, VARIANT Arg18, VARIANT Arg19, VARIANT Arg20, VARIANT Arg21, VARIANT Arg22, VARIANT Arg23, VARIANT Arg24, VARIANT Arg25, VARIANT Arg26, VARIANT Arg27, VARIANT Arg28, VARIANT Arg29, VARIANT Arg30, VARIANT Arg31, VARIANT Arg32, IDispatch** pVal)
{
	begin_function	
	HRESULT								hr;	
	CComPtr<IDispatch>					spObject;
	
	if (hr = CreateObject(estring(ObjectName), &spObject, 32, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7, Arg8, Arg9, Arg10, Arg11, Arg12, Arg13, Arg14, Arg15, Arg16, Arg17, Arg18, Arg19, Arg20, Arg21, Arg22, Arg23, Arg24, Arg25, Arg26, Arg27, Arg28, Arg29, Arg30, Arg31, Arg32)) return hr;
	return spObject.CopyTo(pVal);
	end_function
}
HRESULT CSiriusApplication::CreateObject(const std::string& szObjectName, IDispatch** pVal, int nArgs, /*VARIANT*/...)
{		
	HRESULT								hr;
	CLSID								clsid;	
	CComVariant							vData;
	va_list								vl;
	CComPtr<IDispatch>					spObject;

	va_start(vl, nArgs);
	if (hr = CParameterMap::VariableArgumentListToArray(&vData, nArgs, vl)) return hr;
	va_end(vl);		
	if (hr = g_pApplication->GetObjectManager().NameToCLSID(szObjectName, &clsid)) return hr;
	if (hr = spObject.CoCreateInstance(clsid)) return hr;
	if (hr = CComDispatchDriverEx(spObject).PutPropertyByName(L"Value", &vData)) return hr;
	return spObject.CopyTo(pVal);
}


///////////////////////////////////////////////////////////////////////////////
//	Enums
//
//	This function returns all the enumerator values associated with an
//	enumator name.
//
STDMETHODIMP CSiriusApplication::Enums(BSTR EnumName, VARIANT* pVal)
{			
	std::vector<std::string>			aszNames;
	
	begin_function
	map_parameter(EnumName, estring, szEnumName);
	CEnumMap::GetEnumList(szEnumName, LIBID_Sirius, &aszNames);
	unmap_parameter(aszNames, pm);
	pm.GetColumn(0, VT_BSTR, pVal);
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	EnumToString
//
//	This function maps a enumerator value to its corresponding string
//
STDMETHODIMP CSiriusApplication::EnumToString(BSTR EnumName, long Value, BSTR* psName)
{
	estring								szName;
	HRESULT								hr;
		
	begin_function
	map_parameter(EnumName, estring, szEnumName);
	if (!psName) return E_POINTER;
	if (hr = CEnumMap::GetString(szEnumName, LIBID_Sirius, Value, &szName)) return hr;
	return szName.GetBSTR(psName);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	FinalConstruct
//
//	Called when an instance of this object is created.
//	There should only be one object since it's a singleton.
//
HRESULT	CSiriusApplication::FinalConstruct(void)
{
	HRESULT								hr;
	CComPtr<IAssets>					spAssets;	
	CComPtr<ICorrelationMatrices>		spCorrelationMatrices;
	CComPtr<IDividendSchedules>			spDividendSchedules;
	CComPtr<IZeroCurves>				spZeroCurves;
	CComPtr<IVolatilityStructures>		spVolatilityStructures;
	CComPtr<ISpotSchedules>				spSpotSchedules;
	
	m_bErrorLogging = false;
	m_bFunctionsDisabled = false;
	
	if (hr = m_spMarketData.CoCreateInstance(CLSID_MarketData)) return hr;
		
	// load the market data objects			
	if (hr = m_spMarketData->get_Assets((IAssets**)&spAssets)) return hr;
	if (hr = m_spMarketData->get_CorrelationMatrices((ICorrelationMatrices**)&spCorrelationMatrices)) return hr;
	if (hr = m_spMarketData->get_DividendSchedules((IDividendSchedules**)&spDividendSchedules)) return hr;
	if (hr = m_spMarketData->get_VolatilityStructures((IVolatilityStructures**)&spVolatilityStructures)) return hr;
	if (hr = m_spMarketData->get_ZeroCurves((IZeroCurves**)&spZeroCurves)) return hr;
	if (hr = m_spMarketData->get_SpotSchedules((ISpotSchedules**)&spSpotSchedules)) return hr;
	if (hr = m_ObjectManager.InstallMarketDataComponent(LIBID_Sirius, CLSID_Asset, IID_IAsset, spAssets, CLSID_Assets, true)) return hr;
	if (hr = m_ObjectManager.InstallMarketDataComponent(LIBID_Sirius, CLSID_CorrelationMatrix, IID_ICorrelationMatrix, spCorrelationMatrices, CLSID_CorrelationMatrices, false)) return hr;
	if (hr = m_ObjectManager.InstallMarketDataComponent(LIBID_Sirius, CLSID_VolatilityStructure, IID_IVolatilityStructure, spVolatilityStructures, CLSID_VolatilityStructures, false)) return hr;
	if (hr = m_ObjectManager.InstallMarketDataComponent(LIBID_Sirius, CLSID_ZeroCurve, IID_IZeroCurve, spZeroCurves, CLSID_ZeroCurves, false)) return hr;
	if (hr = m_ObjectManager.InstallMarketDataComponent(LIBID_Sirius, CLSID_DividendSchedule, IID_IDividendSchedule, spDividendSchedules, CLSID_DividendSchedules, true)) return hr;
	if (hr = m_ObjectManager.InstallMarketDataComponent(LIBID_Sirius, CLSID_SpotSchedule, IID_ISpotSchedule, spSpotSchedules, CLSID_SpotSchedules, true)) return hr;
	
	// load the product type objects (any error in LoadProductTypes is reported in a dialog)
	LoadProductTypes();
			
	// load the standard objects
	hr = InstallStandardObjects();	
	return hr;
}


/////////////////////////////////////////////////////////////////////////////
//	FinalRelease
//
//	called when the final reference on the object has been released
//
void CSiriusApplication::FinalRelease(void)
{
	ATLTRACE("CSiriusApplication::FinalRelease() called\n");	
}


/////////////////////////////////////////////////////////////////////////////
//	GetAssets
//
//	Returns a pointer to the currently selected assets collection.
//
/*static*/ HRESULT CSiriusApplication::GetAssets(CComPtr<IAssets>& spAssets)
{
	HRESULT								hr;
	CComPtr<ISiriusApplication>			spApplication;
	CComPtr<IMarketData>				spMarketData;	
						
	if (hr = _Module.GetSiriusApplication()->get_MarketData(&spMarketData)) return hr;
	if (hr = spMarketData->get_Assets(&spAssets)) return hr;	
	return S_OK;
}
/*static*/ CAssets* CSiriusApplication::GetAssets(void)
{
	CComPtr<IAssets>					spMaintainedAssets;				// maintained assets collection	

	if (g_pApplication->GetAssets(spMaintainedAssets)) propagate_error;
	return dynamic_cast<CAssets*>(spMaintainedAssets.p);	
}


/////////////////////////////////////////////////////////////////////////////
//	GetCorrelation
//
//	Returns the correlation number defined between two assets. We make
//  database calls if necessary.
//
double CSiriusApplication::GetCorrelation(const std::string& szAsset1, const std::string& szAsset2, const std::string& szDataSourceOpt, long nDate) const
{
	// If the data source is blank or is set to Last or PL then first try obtaining the correlation from
	// the appropriate correlation matrix.
	HRESULT								hr;
	DataSourceEnum						ds;
	std::string							szError;
		
	if (!szDataSourceOpt.size()){
		ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
		hr = S_OK;
	} else {
		hr = CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, szDataSourceOpt, NoDataSource, &ds);
	}
	if (hr == S_OK){
		// Get the appropriate correlation matrix from the market data set.
		MlEqCorrelationMatrixHandle h = GetCorrelationMatrixHandle(ds, nDate, true);
		try {
			double f = h->GetCorrelation(szAsset1, szAsset2);
			return f;
		} catch (const std::string& sz){
			szError.assign(sz);
		}
	}

	// Try obtaining the correlation from the database.
	CComPtr<IDispatch>					spObject;
	CComQIPtr<ICorrelationMatrix>		spCorrelationMatrix;
	try {
		_Module.Load(CLSID_CorrelationMatrix, szAsset1 + "," + szAsset2, szDataSourceOpt, nDate, true, spObject);
	} catch (...){
		throw szError;
	}
	if (!(spCorrelationMatrix = spObject)) throw szError;
	map_com_to_analytic(spCorrelationMatrix, CorrelationMatrix, h);
	double f = h->GetCorrelation(szAsset1, szAsset2);
	return f;
}


/////////////////////////////////////////////////////////////////////////////
//	GetCorrelationMatrix
//
//	Returns the correlation matrix defined on the given date and data source
//  in the market data collection. Returns zero if there is no matrix defined.
//
/*static*/ MlEqCorrelationMatrixHandle CSiriusApplication::GetCorrelationMatrixHandle(DataSourceEnum ds, long nDate, bool bCreateIfNecessary)
{
	map_com_to_analytic(GetCorrelationMatrix(ds, nDate, bCreateIfNecessary), CorrelationMatrix, h);
	return h;
}
/*static*/ CComPtr<ICorrelationMatrix> CSiriusApplication::GetCorrelationMatrix(DataSourceEnum ds, long nDate, bool bCreateIfNecessary)
// bCreateIfNecessary - true if we create the matrix if one is not available in the market data set.
{
	CComPtr<IMarketData>				spMarketData;
	CComPtr<ICorrelationMatrices>		spCorrelationMatrices;
	CComPtr<ICorrelationMatrix>			spCorrelationMatrix;
	CCorrelationMatrices*				pCorrelationMatrices;
	
	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();
	
	CComObjectCollectionSerialisableKey key(CCorrelationMatrix::s_szCorrelationMatrixName, nDate, ds);
	
	if (_Module.GetSiriusApplication()->get_MarketData(&spMarketData)) propagate_error;
	if (spMarketData->get_CorrelationMatrices(&spCorrelationMatrices)) propagate_error;
	pCorrelationMatrices = dynamic_cast<CCorrelationMatrices*>(spCorrelationMatrices.p);
	if (pCorrelationMatrices->IsInCollection(key)){
		if (spCorrelationMatrices->get_Item(key, &spCorrelationMatrix)){
			propagate_error;
		}
	} else {
		if (bCreateIfNecessary){
			spCorrelationMatrix.CoCreateInstance(CLSID_CorrelationMatrix);
			spCorrelationMatrix->put_DataSource(ds);
			spCorrelationMatrix->put_Date(nDate);
			spCorrelationMatrices->Add(CComVariant(), spCorrelationMatrix);
		} else {
			throw "Cannot find the correlation matrix on date '" + MlEqDate(nDate).GetString() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
		}
	}
	return spCorrelationMatrix;
}


/////////////////////////////////////////////////////////////////////////////
//	GetCreateParameterNames
//
//	Returns the parameter names for the MLCreate function associated with an
//	input object.
//
HRESULT CSiriusApplication::GetCreateParameterNames(CComDispatchDriverEx& ddObject, const std::vector<CParameterMap>& vpm, std::vector<std::string>* pvector, std::map<std::string, long>* pmap, std::map<std::string, long>* pmapCase) const
//	pvector - vector of parameter names (returned, nullable)
//	pmap - map of lower case and trimmed parameter names and the corresponding parameter number; one-based on the leftmost element (returned, nullable)
//	pmapCase - as for pmap but we return case sensitive values (returned, nullable)
{
	HRESULT								hr;
	std::string							szObjectName;
	CComPtr<ITypeInfo>					spti;
	CComPtr<ITypeLib>					sptl;
	MEMBERID							memid;	
	unsigned short						nFind = 1;	
	OLECHAR*							sFunctionName;
	long								nParameter = 0;
	
	if (pvector) pvector->clear();
	if (pmap) pmap->clear();
	if (pmapCase) pmapCase->clear();
	if (hr = CParameterMap::GetObjectName(ddObject, &szObjectName)) return hr;
	if (hr = ddObject.p->GetTypeInfo(0, LOCALE_USER_DEFAULT, &spti)) return hr;
	if (hr = spti->GetContainingTypeLib(&sptl, NULL)) return hr;
	spti = NULL;
	sFunctionName = estring::AnsiToUnicode("MLCreate" + szObjectName);
	HRESULT hrHash = LHashValOfNameSys(SYS_WIN32, LOCALE_USER_DEFAULT, sFunctionName);
	sptl->FindName(sFunctionName, hrHash, &spti, &memid, &nFind);
	delete sFunctionName;	
	if (!nFind || !spti){
		// Could not find a simple function. Try to find szObjectName + "TypeEnum" by examining a property of ddObject called szObjectName + "Type".
		std::string						szEnum;
		long							nEnum;
		if (vpm.size() < 1 || vpm[0].GetValue(&nEnum)) return E_FAIL;				
		if (hr = CEnumMap::GetString(szObjectName + "TypeEnum", LIBID_Sirius, nEnum, &szEnum)) return hr;		
		sFunctionName = estring::AnsiToUnicode("MLCreate" + szEnum + szObjectName);
		nParameter++;
		if (pvector) pvector->push_back(szObjectName + "Type");		
		if (pmap || pmapCase){
			estring sz(szObjectName + "Type");
			sz.trim();
			if (pmapCase)(*pmapCase)[sz] = nParameter;
			sz.lc();			
			if (pmap)(*pmap)[sz] = nParameter;			
		}
		HRESULT hrHash = LHashValOfNameSys(SYS_WIN32, LOCALE_USER_DEFAULT, sFunctionName);
		nFind = 1;
		sptl->FindName(sFunctionName, hrHash, &spti, &memid, &nFind);
		delete sFunctionName;		
	}
	if (!nFind || !spti) return E_FAIL;

	BSTR*								pb = new BSTR[256];
	unsigned int						nNames;
	
	spti->GetNames(memid, pb, 256, &nNames);
	ATLASSERT(nNames < 256);			// Need to make the (hard-coded) array larger. Getting the size correct is possible, but not computationally worthwhile.
	for (unsigned int nName = 1; nName < nNames - 1 /*i.e. we don't include the final parameter which is the function return parameter name*/; nName++){
		nParameter++;
		estring sz(pb[nName]);
		::SysFreeString(pb[nName]);
		if (pvector) pvector->push_back(sz);
		if (pmap || pmapCase){
			sz.trim();
			if (pmapCase)(*pmapCase)[sz] = nParameter;
			sz.lc();			
			if (pmap)(*pmap)[sz] = nParameter;			
		}
	}
	delete pb;
	return S_OK;
}


/////////////////////////////////////////////////////////////////////////////
//	GetFunctionsDisabled
//
//	Returns the disabled / enabled state of the ML[...] functions 
//
STDMETHODIMP CSiriusApplication::get_FunctionsDisabled(VARIANT_BOOL *pVal)
{
	*pVal = m_bFunctionsDisabled;
	return S_OK;
}
bool CSiriusApplication::GetFunctionsDisabled(void) const
{
	return m_bFunctionsDisabled;
}


////////////////////////////////////////////////////////////////////////////
//	GetLocations
//
//	Retrieve the list of trading locations from a database.
//	The return value is the connection pointer associated with the database used.
//
_ConnectionPtr CSiriusApplication::GetLocations(FileSystemEnum fs, const CSiriusComModule::connection_parameters* pcp, std::vector<estring>* pasz, std::map<estring, estring>* pmap, bool bReturnLowerCase, bool bThrow) const
//	fs - input file system
//	pcp - (nullable), specifies an alternative database to use
//	pasz, pmap - receiver variables (returned, nullable), the RHS of the map is always the proper case
{
	if (fs == fsSQLServer){
		_ConnectionPtr						spConnection;
		std::stringstream					ssQuery;						// SQL query to execute
		_RecordsetPtr						prs;							// ADO recordset
		CComVariant							vField_ID;						// position field ID
		CComVariant							vField_Value;					// position field value		
				
		if (pasz) pasz->clear();
		if (pmap) pmap->clear();
		ssQuery << "sp_user_get_locations";
		try {
			_Module.GetConnection(spConnection, pcp);
			prs = spConnection->Execute(L"sp_user_get_locations", NULL, -1);
			if (prs->adoEOF){
				if (bThrow)	throw "No locations defined in the database";	
				return spConnection;
			}
			while (!prs->adoEOF){				
				estring sz((CComVariant)prs->GetFields()->GetItem(0L)->GetValue());
				estring sz_lc(sz);
				if (bReturnLowerCase) sz_lc.lc();
				if (pmap) (*pmap)[sz_lc] = sz;
				if (pasz) pasz->push_back(sz_lc);
				prs->MoveNext();
			}
			return spConnection;
		} catch (_com_error& e){
			if (bThrow) throw estring(e);			
		}	
	} else if (fs == fsNTFS){
		if (bThrow) throw "Locations are not explicitly defined in the local file system case";		
	}
	return NULL;
}
STDMETHODIMP CSiriusApplication::get_Locations(VARIANT* pVal)
{
	begin_function
	std::vector<estring>				asz;	
	GetLocations(_Module.GetFileSystem(), NULL, &asz, NULL, false, true);
	unmap_parameter(asz, pm);
	return pm.GetValue(pVal);
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	GetMarketData
//
//	returns an interface to the currently selected market data object
//
STDMETHODIMP CSiriusApplication::get_MarketData(IMarketData** pVal)
{	
	return m_spMarketData.CopyTo(pVal);
}


/////////////////////////////////////////////////////////////////////////////
//	GetOneParameter
//
//	Returns the value of a named parameter of an MLCreate[...] function of
//	a given object.
//
HRESULT CSiriusApplication::GetOneParameter(CComDispatchDriverEx& ddObject, std::string szParameter, CComVariant* pv)
{
	HRESULT											hr;
	std::vector<CParameterMap>						vpm;	
	std::map<std::string, long>						map;	
	std::map<std::string, long>::const_iterator		it;	
		
	if (hr = CParameterMap::ArrayToVector(*pv, NULL, &vpm)) return hr;
	if (hr = GetCreateParameterNames(ddObject, vpm, NULL, &map, NULL)) return hr;
	estring::lc(&szParameter);
	estring::trim(&szParameter);
	if ((it = map.find(szParameter)) == map.end() || it->second - 1 < 0 || it->second - 1 >= vpm.size()) return E_FAIL;
	return vpm[it->second - 1].GetValue(pv);	
}


/////////////////////////////////////////////////////////////////////////////
//	GetObject
//
//	Returns an object associated with an input handle. We allow cases where
//  the 'handle' is infact a disaptch pointer to an object of the desired 
//  type, or comprises the variant data associated with an object.
//
STDMETHODIMP CSiriusApplication::GetObject(VARIANT Handle, BSTR ObjectName, IDispatch** pVal)
//	Handle - this is normally an Excel-style handle but it could be a pointer
//  to an object or even an array containing a scalar. We have to allow for all these
//  possibilities.
{
	begin_function
	CComPtr<IDispatch>				spDispatch;
	HRESULT							hr;
	CParameterMap					pmHandle;
	estring							szHandle;
	
	map_optional_parameter(ObjectName, estring, szObjectName, "");
	if (pmHandle.SetValue(Handle)) return CParameterMap::ReturnErrorR(IDS_NO_OBJECT_WITH_HANDLE);
	
	if (!pmHandle.IsScalar()){
		// Handle could a represent the object's data.				
		if (!m_ObjectManager.GetObject(pmHandle.GetValue(), szObjectName, spDispatch)){
			return spDispatch.CopyTo(pVal);
		}
		return CParameterMap::ReturnErrorR(IDS_NO_OBJECT_WITH_HANDLE);
	}
	if (pmHandle.GetElementPtr(0, 0)->vt == VT_DISPATCH){				
		return CComPtr<IDispatch>(pmHandle.GetElementPtr(0, 0)->pdispVal).CopyTo(pVal);	// A bit longwinded but guarantees reference counting!
	}
	if (pmHandle.GetValue(0, 0, &szHandle)){
		return CParameterMap::ReturnErrorR(IDS_NO_OBJECT_WITH_HANDLE);
	}
	if (hr = m_ObjectManager.GetObject(szHandle, szObjectName, spDispatch)) return hr;
	return spDispatch.CopyTo(pVal);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	GetObjects
//
//	Returns an array corresponding to the objects associated with an input
//	array of handles and / or object pointers.
//
//	The output array is the same size as the (non-blank region of) the input.
//
STDMETHODIMP CSiriusApplication::GetObjects(VARIANT Handles, BSTR ObjectName, VARIANT* pObjectsArray)
{		
	begin_function
	map_optional_parameter(ObjectName, estring, szObjectName, "");	
	if (Handles.vt == VT_DISPATCH){
		// Pass straight through subject to szObjectName constraint.
		estring szName;
		CParameterMap::GetObjectName(CComPtr<IDispatch>(Handles.pdispVal), &szName);
		if (m_ObjectManager.IsSiriusSingular(szName)){
			if (!szObjectName.size() || !szName.CompareNoCase(szObjectName)){
				return ::VariantCopy(pObjectsArray, &Handles);				
			} else {
				throw "Invalid object. Expected '" + szName + "', found '" + szObjectName + "'";
			}
		} else if (m_ObjectManager.IsSiriusCollection(szName)){
			CComDispatchDriverEx dd(Handles.pdispVal);
			CComVariant vCount;
			dd.GetPropertyByName(L"Count", &vCount);
			CParameterMap pm;
			pm.SetSize(vCount.lVal, 1);						
			for (CComVariant vItem = 1L; vItem.lVal <= vCount.lVal; vItem.lVal++){
				CComVariant vObject;
				dd.GetPropertyByName(L"Item", &vItem, 1, &vObject);
				ATLASSERT(vObject.vt == VT_DISPATCH);
				pm.SetValue(vItem.lVal - 1, 0, vObject);
			}
			pm.GetColumn(0, VT_DISPATCH, pObjectsArray);
			return S_OK;
		}
	}

	HRESULT								hr;	
	CParameterMap						pmHandles;
	CParameterMap						pmObjects;						// objects to return

	if (hr = pmHandles.SetValue(Handles)) return hr;
	pmObjects.SetSize(pmHandles.GetRows(), pmHandles.GetCols());
	for (long nRow = 0; nRow < pmHandles.GetRows(); nRow++){
		for (long nCol = 0; nCol < pmHandles.GetCols(); nCol++){
			CComPtr<IDispatch> spObject;			
			if (hr = m_ObjectManager.GetObject(*pmHandles.GetElementPtr(nRow, nCol), szObjectName, spObject)) return hr;
			pmObjects.SetValue(nRow, nCol, spObject.p);
		}
	}
	return pmObjects.GetValue(pObjectsArray);
	end_function	
}


/////////////////////////////////////////////////////////////////////////////
//	GetObjectManager
//
//	Returns a constant reference to the object manager object.
//
CObjectManager& CSiriusApplication::GetObjectManager(void)
{
	return m_ObjectManager;
}


/////////////////////////////////////////////////////////////////////////////
//	GetObjectTypes
//
//	return an array of installed object names and their GUIDs
//
STDMETHODIMP CSiriusApplication::GetObjectTypes(VARIANT *pVal)
{	
	return m_ObjectManager.GetObjectListOfType(CSiriusSingularObject::StandardObject | CSiriusSingularObject::MarketDataComponent, pVal);
}


/////////////////////////////////////////////////////////////////////////////
//	GetProductTypes
//
//	return an array of installed product type names
//
STDMETHODIMP CSiriusApplication::GetProductTypes(VARIANT *pVal)
{		
	return m_ObjectManager.GetObjectListOfType(CSiriusSingularObject::ProductType, pVal);
}


/////////////////////////////////////////////////////////////////////////////
//	GetProperties
//
//	Returns a list of possible parameter values for MLGetValue.
//
HRESULT CSiriusApplication::GetProperties(CComPtr<IDispatch> spObject, VARIANT* pResult)
{
	HRESULT								hr;	
	CParameterMap						pm;
	std::map<std::string, long>			map1, map2;	
	std::vector<CParameterMap>			vpm;	
	CComDispatchDriverEx				ddObject(spObject);
	CComVariant							v;

	if (hr = m_ObjectManager.GetObjectProperties(INVOKE_PROPERTYGET, spObject, NULL, &pm)) return hr;	
	if (hr = pm.GetValue(0, -1, map1)) return hr;

	// Attempt to get the properties that are parameter names of an appropriate MLCreate[...] function.
	if (!ddObject.GetPropertyByName(L"Value", &v) && !CParameterMap::ArrayToVector(v, NULL, &vpm) && !GetCreateParameterNames(ddObject, vpm, NULL, NULL, &map2)){
		CParameterMap::MergeMap(&map1, map2);
	}
	if (hr = pm.SetValue(map1)) return hr;
	if (hr = pm.RemoveColumn(1)) return hr;
	return pm.GetValue(pResult);
}
								  

/////////////////////////////////////////////////////////////////////////////
//	GetSetup
//
//	returns a pointer to the setup singleton object
//
STDMETHODIMP CSiriusApplication::get_Setup(ISetup** pVal)
{
	HRESULT								hr;
	CComPtr<ISetup>						spSetup;
	if (hr = spSetup.CoCreateInstance(CLSID_Setup)) return hr;
	return spSetup.CopyTo(pVal);
}


/////////////////////////////////////////////////////////////////////////////
//	GetVersion
//
//	Returns the version number of the module.
//
HRESULT CSiriusApplication::GetVersion(const std::string& szModuleName, std::string* psz) const
{
	//extracts product version from file
	std::stringstream ssVersion;
	if (HMODULE hModule = ::GetModuleHandle(szModuleName.c_str())){
		TCHAR szDllName[MAX_PATH + 1];
		GetModuleFileName(hModule, szDllName, sizeof(szDllName)); 
		DWORD dwHnd;
		if (DWORD dwVerInfoSize = GetFileVersionInfoSize(szDllName, &dwHnd)){
			void* pBuffer;
			if((pBuffer = malloc(dwVerInfoSize))==0L){
				return E_FAIL;
			}
			GetFileVersionInfo(szDllName, dwHnd, dwVerInfoSize, pBuffer); 
			VS_FIXEDFILEINFO *pFixedInfo;
			UINT uVersionLen;
			VerQueryValue(pBuffer,_T("\\"),(void**)&pFixedInfo,(UINT *)&uVersionLen);
			ssVersion << HIWORD(pFixedInfo->dwProductVersionMS)
			   << "." << LOWORD(pFixedInfo->dwProductVersionMS)
			   << "." << HIWORD(pFixedInfo->dwProductVersionLS)
			   << "." << LOWORD(pFixedInfo->dwProductVersionLS);
			free(pBuffer); 
		}
	}
	psz->assign(ssVersion.str());
	return S_OK;
}

STDMETHODIMP CSiriusApplication::GetVersion(BSTR ModuleName, BSTR* pVersion)
{	
	HRESULT								hr;
	estring								szVersion;
	
	if (hr = GetVersion(estring(ModuleName), &szVersion)) return hr;
	return szVersion.GetBSTR(pVersion);
}


/////////////////////////////////////////////////////////////////////////////
//	GetZeroCurves
//
//	Returns a pointer to the currently selected zero curve collection.
//
/*static*/ HRESULT CSiriusApplication::GetZeroCurves(CComPtr<IZeroCurves>& sp)
{
	HRESULT								hr;	
	CComPtr<IMarketData>				spMarketData;	
						
	if (hr = _Module.GetSiriusApplication()->get_MarketData(&spMarketData)) return hr;
	if (hr = spMarketData->get_ZeroCurves(&sp)) return hr;	
	return S_OK;
}


/////////////////////////////////////////////////////////////////////////////
//	InsertCorrelation
//
//	Ensures that a correlation pair is present in a correlation matrix in
//  the collection maintained by the object manager.
//
/*static*/ void CSiriusApplication::InsertCorrelation(const std::string& sz1, const std::string& sz2, DataSourceEnum ds, long nDate, bool bThrow)
{			
	CComPtr<ICorrelationMatrix>			spCorrelationMatrix = GetCorrelationMatrix(ds, nDate, true);
	CCorrelationMatrix*					pCorrelationMatrix;

	try {
		if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
		if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();
		pCorrelationMatrix = dynamic_cast<CCorrelationMatrix*>(spCorrelationMatrix.p);
		pCorrelationMatrix->GetCorrelation(sz1, sz2);	// this will add the sz1, sz2 pair to the matrix if necessary
	} catch (const std::string& sz){
		if (bThrow) throw sz;
	}
}


/////////////////////////////////////////////////////////////////////////////
//	InsertObject
//
//	Inserts an object into the appropriate collection. A suitable handle
//	is returned.
//
STDMETHODIMP CSiriusApplication::InsertObject(IDispatch* Object, VARIANT_BOOL ExternalHandle, BSTR* pVal)
{		
	begin_function
	map_parameter(ExternalHandle, bool, bExternal)	
	return m_ObjectManager.InsertObject(CComPtr<IDispatch>(Object), bExternal).GetBSTR(pVal);	
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	InterfaceSupportsErrorInfo
//
//	used to return errors by VB clients etc.
//
STDMETHODIMP CSiriusApplication::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_ISiriusApplication };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}


///////////////////////////////////////////////////////////////////////////////
//	Load
//
//	Loads an object from the database
//
STDMETHODIMP CSiriusApplication::Load(BSTR ObjectName, BSTR Identifier, VARIANT DataSourceOpt, VARIANT DateOpt, IDispatch** pVal)
//	ObjectName - either the name of the type of object to load, or its ProgID
{	
	CComPtr<IDispatch>					spObject;
	HRESULT								hr;
	
	begin_function
	map_parameter(ObjectName, estring, szObjectName);
	map_optional_parameter(Identifier, estring, szIdentifier, "");
	map_optional_parameter(DataSourceOpt, estring, szDataSource, "");
	map_optional_parameter(DateOpt, DATE, date, 0.0);
	if (hr = _Module.Load(szObjectName, szIdentifier, szDataSource, date, true, spObject)) return hr;
	return spObject.CopyTo(pVal);
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	LoadProductTypes
//
//	attempt to load the product types from the appropriate file
//
HRESULT CSiriusApplication::LoadProductTypes(void)
{
	HRESULT								hr;
	std::string							szFileName = _Module.GetProductTypesFileName();
	std::string							szLIBID;
	std::string							szCLSID;
	std::string							szIID;
	std::string							szName;
	
	if (!szFileName.size()){
		CParameterMap::DisplayError(IDS_INVALID_PRODUCT_TYPE_FILE, MB_ICONSTOP);
		return E_FAIL;
	}
	m_ObjectManager.Clear(CSiriusSingularObject::ProductType);
	if (m_ObjectManager.InstallProductType(LIBID_Sirius, CLSID_Swap, IID_ISwap)) ATLASSERT(false);
	std::ifstream file;
	file.open(szFileName.c_str(), std::ios::in);
	if (file.fail()){
		CParameterMap::DisplayError(IDS_INVALID_PRODUCT_TYPE_FILE, MB_ICONSTOP);
		return E_FAIL;
	}
	while (!file.fail() && !file.eof()){
		estring szLine;
		std::getline(file, szLine);
		szLine.trim();
		if (szLine == "<ProductTypes>" || szLine == "</ProductTypes>" || !szLine.length()) continue;
		szLIBID = szLine.esubstr("LIBID=\"", "\"");
		szCLSID = szLine.esubstr("CLSID=\"", "\"");
		szIID = szLine.esubstr("IID=\"", "\"");
		szName = szLine.esubstr(">", "<");
		estring::StripWhiteSpace(&szCLSID);
		estring::StripWhiteSpace(&szIID);
		estring::StripWhiteSpace(&szName);
		if (!szCLSID.length() || !szIID.length() || !szName.length()) continue;
		if (hr = m_ObjectManager.InstallProductType(szLIBID, szCLSID, szIID, szName)){
			CParameterMap::DisplayError(MB_ICONSTOP);
			return hr;
		}
	}
	return S_OK;
}
															

/////////////////////////////////////////////////////////////////////////////
//	PutFunctionsDisabled
//
//	Toggles the disabling flag of all ML[...] functions
//
STDMETHODIMP CSiriusApplication::put_FunctionsDisabled(VARIANT_BOOL newVal)
{
	m_bFunctionsDisabled = newVal ? true : false;
	return S_OK;
}
void CSiriusApplication::PutFunctionsDisabled(bool b)
{
	m_bFunctionsDisabled = b;
}


///////////////////////////////////////////////////////////////////////////////
//	PutMarketData
//
//	switches the market data object attached to this SiriusApplication
//
STDMETHODIMP CSiriusApplication::put_MarketData(IMarketData* newVal)
{
	begin_function
	if (!newVal){
		throw CStringResource(IDS_CANNOT_NULLIFY_MARKET_DATA);
	}	
	// ToDo - if we enable this then we will have to synchronise with
	// the maintained collection objects in the object manager AND adjust
	// the code in CComObjectCollectionSerialisable.h to support multiple
	// instances of the same object key across different market data objects.
	throw "Changing the market data object is not supported in this release. Speak to David Cuin if you nequire it.";
	m_spMarketData = newVal;
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	CreateObjectFromXml
//
//	Makes an object out of the input Xml and returns a handle to it.
//
STDMETHODIMP CSiriusApplication::CreateObjectFromXml(/*[in]*/ BSTR Xml, /*[out, retval]*/ IDispatch** pVal)
{
	begin_function	
	CComVariant							vObject;
	CComBSTR							sHandle;
		
	if (CXmlStreamer::GetVariant(estring(Xml).data(), vObject)) propagate_error;
	if (vObject.vt != VT_DISPATCH) throw "The input Xml does not represent an object";	
	return CComPtr<IDispatch>(vObject.pdispVal).CopyTo(pVal);	
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	Save
//
//	saves a handle to the database
//
STDMETHODIMP CSiriusApplication::Save(BSTR Handle)
{
	begin_function	
	CComPtr<IDispatch>					spObject;	
	estring								szHandle(Handle);		
	DISPID								dispid;
	HRESULT								hr;
				
	szHandle.trim();
	if (!szHandle.size()) return CParameterMap::ReturnErrorR(IDS_INVALID_HANDLE_PARAMETER, IID_ISiriusApplication);
	if (hr = GetObject(szHandle, L"", &spObject)) return hr;
	
	// find the Save function
	if (CComDispatchDriverEx(spObject).GetIDOfName(L"Save", &dispid)){
		std::string szName;
		CParameterMap::GetObjectName(spObject, &szName);
		return CParameterMap::ReturnErrorRS(IDS_SAVE_NOT_FOUND, szName, IID_ISiriusApplication);
	}
	// then invoke it		
	return CComDispatchDriverEx(spObject).Invoke0(dispid);
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	StringToEnum
//
//	This function maps a string value to its corresponding enumerator
//
STDMETHODIMP CSiriusApplication::StringToEnum(BSTR EnumName, BSTR Name, long* pnValue)
{
	HRESULT								hr;
	long								nValue;
	
	begin_function
	map_parameter(EnumName, estring, szEnumName);
	map_parameter(Name, estring, szName);
	if (!pnValue) return E_POINTER;
	if (hr = CEnumMap::GetEnum(szEnumName, LIBID_Sirius, szName, &nValue)) return hr;
	*pnValue = nValue;
	end_function
}


//////////////////////////////////////////////////////////////////////////////
//	UpdatePL
//
//	Updates the PL data source tables with the data currently defined in
//	the tables associated with the data source Last.
//
STDMETHODIMP CSiriusApplication::UpdatePL(DATE CopyDateOpt)
{
	FileSystemEnum						fs;
	long								nDate = CopyDateOpt;
	
	begin_function
	check_publishing_enabled
	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;						// SQL query to execute	
		_RecordsetPtr					prs;							// ADO recordset
		
		if (!nDate){
			ssQuery << "sp_user_update_pl '" << _Module.GetLocation() << "'";
		} else {
			ssQuery << "sp_user_update_pl '" << _Module.GetLocation() << "', " << nDate;
		}
		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
		} catch (_com_error& e){
			throw estring(e);
		}
	} else if (fs == fsNTFS){						
		HRESULT							hr;		
	
		// Emulate the stored procedure 'sp_user_update_pl' by loading 
		// objects on today's date from the data source 'Last' and saving
		// them to the 'PL' data source.
		
		if (!nDate) nDate = MlEqDate::GetCurrentDate();		
		if (hr = CComObjectCollectionFunctions<ISpotSchedule>::UpdatePL_NTFS(0, CSpotSchedules::Load)) return hr;
		if (hr = CComObjectCollectionFunctions<IDividendSchedule>::UpdatePL_NTFS(nDate, CDividendSchedules::Load)) return hr;
		if (hr = CComObjectCollectionFunctions<ICorrelationMatrix>::UpdatePL_NTFS(nDate, CCorrelationMatrices::Load)) return hr;
		if (hr = CComObjectCollectionFunctions<IZeroCurve>::UpdatePL_NTFS(nDate, CZeroCurves::Load)) return hr;
		if (hr = CComObjectCollectionFunctions<IVolatilityStructure>::UpdatePL_NTFS(nDate, CVolatilityStructures::Load)) return hr;
		if (hr = CComObjectCollectionFunctions<IAsset>::UpdatePL_NTFS(nDate, CAssets::Load)) return hr;
		if (hr = CComObjectCollectionFunctions<IPosition>::UpdatePL_NTFS(nDate, CPositions::Load)) return hr;
	}
	end_function
}


//////////////////////////////////////////////////////////////////////////////
//	Error handling functions
//
//  These functions should be the only COM interface implementation functions
//	that do not use begin_function and end_function. They also can never
//	generate a string error as doing so would cause infinite recursion.
//
STDMETHODIMP CSiriusApplication::AddError(BSTR Error)
{
	estring								szError(Error);			
	szError.trim();
	Fire_OnError(Error); 
	if (!m_bErrorLogging || !szError.size()) return E_FAIL;		
	m_vectorErrors.push_back(szError);	
	return S_OK;
}
STDMETHODIMP CSiriusApplication::get_FirstError(BSTR *pVal)
{	
	if (!m_vectorErrors.size()) return E_FAIL;
	return m_vectorErrors[0].GetBSTR(pVal);	
}
STDMETHODIMP CSiriusApplication::get_Errors(VARIANT *pVal)
{
	CParameterMap						pm;
	HRESULT								hr;
		
	if (hr = pm.SetValue(m_vectorErrors)) return hr;
	return pm.GetValue(pVal);	
}
STDMETHODIMP CSiriusApplication::get_LastError(BSTR* pVal)
{
	if (!m_vectorErrors.size()) return E_FAIL;
	return m_vectorErrors[m_vectorErrors.size() - 1].GetBSTR(pVal);
}
void CSiriusApplication::GetLastError(VARIANT* pResult)
{
	if (!m_vectorErrors.size()) return;
	m_vectorErrors[m_vectorErrors.size() - 1].GetValue(pResult);
}
STDMETHODIMP CSiriusApplication::StartErrors(void)
{	
	m_vectorErrors.clear();
	m_bErrorLogging = true;		
	return S_OK;
}
STDMETHODIMP CSiriusApplication::StopErrors(void)
{
	m_bErrorLogging = false;
	return S_OK;
}