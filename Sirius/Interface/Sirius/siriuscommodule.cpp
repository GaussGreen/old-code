//	siriuscommodule.cpp: implementation of the CSiriusComModule class.
//
//	Author:				 David Cuin
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "siriuscommodule.h"
#include "resource.h"
#include "xmlstreamer.h"
#include "siriusapplication.h"
#include "excelinterface.h"
#include "StdString.h"
#include "mleqobjects.h"
#include "handle_asset.h"
#include "handle_assets.h"
#include "handle_correlationmatrix.h"
#include "handle_deal.h"
#include "handle_deals.h"
#include "handle_dividendschedule.h"
#include "handle_dividendschedules.h"
#include "handle_position.h"
#include "handle_positions.h"
#include "handle_spotschedule.h"
#include "handle_spotschedules.h"
#include "handle_volatilitystructure.h"
#include "handle_volatilitystructures.h"
#include "handle_zerocurve.h"
#include "handle_zerocurves.h"
#include "comobjectcollectionserialisabledefaulter.h"
#include "comobjectcollectionserialisablekey.h"

/*static*/ const char						CSiriusComModule::s_szODBCDriverName[] = "SQL Server";
/*static*/ const char						CSiriusComModule::s_sz_registry_path[] = "Software\\Sirius\\Setup";
/*static*/ const char						CSiriusComModule::s_sz_registry_product_types[] = "ProductTypesFile";
/*static*/ const char						CSiriusComModule::s_sz_registry_sirius_sql_server[] = "SQLServerDSN";
/*static*/ const char						CSiriusComModule::s_sz_registry_file_root[] = "FileRoot";
/*static*/ const char						CSiriusComModule::s_sz_registry_file_system[] = "FileSystem";
/*static*/ const char						CSiriusComModule::s_sz_registry_display_sirius_status[]= "DisplaySiriusStatus";
/*static*/ const char						CSiriusComModule::s_sz_registry_product_only[]= "DisplayProductOnly";
/*static*/ const char						CSiriusComModule::s_sz_registry_enable_publishing[] = "EnablePublishing";
/*static*/ const char						CSiriusComModule::s_sz_registry_sirius_user_name[] = "UserName";
/*static*/ const char						CSiriusComModule::s_sz_registry_sirius_password[] = "Password";
/*static*/ const char						CSiriusComModule::s_sz_registry_database_mode[] = "DatabaseMode";
/*static*/ const char						CSiriusComModule::s_sz_registry_use_bwl_for_spot_schedule[] = "UseBwlForSpotSchedule";
/*static*/ const char						CSiriusComModule::s_sz_registry_config_file_path[] = "ConfigFilesPath";
/*static*/ const char						CSiriusComModule::s_sz_registry_use_yesterday_close[] = "UseYesterdaysClose";
/*static*/ const char						CSiriusComModule::s_sz_registry_location[] = "Location";
/*static*/ const char						CSiriusComModule::s_sz_registry_location_first_try[] = "LocationFirstTry";
/*static*/ const char						CSiriusComModule::s_sz_registry_location_first_try_enabled[] = "LocationFirstTryEnabled";
/*static*/ const char						CSiriusComModule::s_sz_default_data_source[] = "DefaultDataSource";
/*static*/ const char						CSiriusComModule::s_sz_registry_bwl_sql_server[] = "Bwl SQL Server";
/*static*/ const char						CSiriusComModule::s_sz_registry_bwl_user_name[] = "Bwl UserName";
/*static*/ const char						CSiriusComModule::s_sz_registry_bwl_password[] = "Bwl Password";
/*static*/ const char						CSiriusComModule::s_sz_registry_bwl_interpolation[] = "Bwl Interpolation";
/*static*/ const char						CSiriusComModule::s_sz_registry_bwl_spot_schedule_source[] = "BwlSpotScheduleSource";


/////////////////////////////////////////////////////////////////////////////
//	Construction
//
CSiriusComModule::CSiriusComModule(void) : m_spConnection(NULL)
{	
	m_hImageList = NULL;
	m_hFont = NULL;	
	RegistryToMembers(&m_state, false);	// 'false' is passed to prevent the attaching of Taurus at this point
	m_state_init = m_state;			
}


/////////////////////////////////////////////////////////////////////////////
//	Load[...]
//
//	These are the set of functions that deal with loading of data from ANY
//	database.
//
//	All requests from Sirius must be routed through either Load or LoadFromBwl.
//
//	LoadFromBwl is used by CSpotSchedule::Load since bwl can be regarded as the
//  'in house' source for Spot Schedules.
//
//	Load is used in all other instances, and, in particular, in
//	CComObjectCollection::get_Item to support 'just-in-time' acquisition of
//	assets.
//
//	As such, LoadFromBwl and Load should be the only public load functions of
//	this class. 
//
//  YOU SHOULD NOT ADJUST THIS ARCHITECTURE WITHOUT SPEAKING FIRST TO DAVID CUIN.
//
/*public*/ HRESULT CSiriusComModule::Load(const std::string& szObjectName, const std::string& szIdentifier, const std::string& szDataSource, DATE date, bool bAllowCloning, CComPtr<IDispatch>& spObject)
//  szObjectName - can also be of the form Library.Object.Version
{
	CLSID							clsid;
	
	if (g_pApplication->GetObjectManager().NameToCLSID(szObjectName, &clsid) && g_pApplication->GetObjectManager().CollectionNameToCollectionCLSID(szObjectName, &clsid)) return CParameterMap::ReturnErrorRS(IDS_LOAD_OBJECT, szObjectName);
	return Load(clsid, szIdentifier, szDataSource, date, bAllowCloning, spObject);	
}
/*public*/ HRESULT CSiriusComModule::Load(CLSID clsid, const std::string& szIdentifier, const std::string& szDataSource, DATE date, bool bAllowCloning, CComPtr<IDispatch>& spObject)
{
	HRESULT								hr;
	
	if (!estring::CompareNoCase("bwldirect", szDataSource)){
		if (hr = LoadFromBwl(clsid, szIdentifier, spObject)) return hr;
	} else if (clsid == CLSID_Asset){
		// The asset case is separated in order to implement the on-the-fly creation
		// of quanto assets.
		// Note that the _LoadAsset function itself calls into _Load
		if (hr = _LoadAsset(szIdentifier, szDataSource, date, bAllowCloning, spObject)) return hr;
	} else if (m_state.m_dbm != UseTaurus && clsid == CLSID_SpotSchedule && !estring::CompareNoCase("bwl", szDataSource)){
		// ToDo - remove this case once Sirius Server is deprecated. The next line enforces this.
		bool bDummy = m_state.m_dbm == UseSiriusServer;
		if (hr = LoadFromBwl(clsid, szIdentifier, spObject)) return hr;
	} else {
		if (hr = _Load(clsid, szIdentifier, szDataSource, date, spObject)) return hr;
	}

	// Add the object to the appropriate maintained collection.
	if (g_pApplication->GetObjectManager().IsCLSIDSiriusCollection(clsid)){
		// Don't do anything with collections as their constituents have already been added.
		return S_OK;	
	} else {
		try {
			g_pApplication->GetObjectManager().InsertObject(spObject, false /*though an MLLoad function might adjust this to an explicit spreadsheet handle*/);			
			// For the correlation matrix case, set the output spObject to the appropriate global correlation.
			if (clsid == CLSID_CorrelationMatrix){
				// Merge the correlation matrix acquired (denoted as spObject) into the corresponding matrix
				// in the object manager's maintained collection.
				CComPtr<ICorrelationMatrix> spCorrelationMatrix = dynamic_cast<ICorrelationMatrix*>(spObject.p);
				map_com_to_analytic(spCorrelationMatrix, CorrelationMatrix, hCorrelationMatrix);				
				spObject = g_pApplication->GetCorrelationMatrix(hCorrelationMatrix->GetDataSource(), hCorrelationMatrix->GetDate(), true);
				return S_OK;	
			}
		} catch (...){
			// Do nothing - InsertObject might fail if we have nominated spObject to not
			// be insertable into the maintained collection. I also don't want to test
			// that in advance of calling InsertObject since I want all that
			// functionality hidden.
		}
		return S_OK;
	}
}
/*public*/ HRESULT CSiriusComModule::LoadFromBwl(CLSID clsid, const std::string& szIdentifier, CComPtr<IDispatch>& spObject)
{	
	CheckBwlConnection();
	HRESULT								hr;
	
	if (clsid == CLSID_SpotSchedule){
		estring							szInterpolation = GetBwlInterpolationRule();
		std::vector<std::string>		aszDataSource;
		CComVariant						vSeries;
		
		GetBwlSpotScheduleSource(&aszDataSource);
		if (!aszDataSource.size()) aszDataSource.resize(1);	// default data source is blank
		for (std::vector<std::string>::const_iterator it = aszDataSource.begin(); it != aszDataSource.end(); it++){
			CComVariant av[] = {szInterpolation.GetValue(), 0L, estring::GetValue(*it), L"", 1L, 0L, 0L, L"Close", estring(szIdentifier), L""};
			if (!(hr = CComDispatchDriverEx(m_spBwlEngine).InvokeN(L"LoadSeries", av, 10, &vSeries))) break;
		}
		if (hr){
			// failed for all the data sources
			throw "Cannot find spot schedule '" + szIdentifier + "' with data source(s) '" + _Module.GetBwlSpotScheduleSource() + "'";
		}
		CComPtr<ISpotSchedule>		spSpotSchedule;
		if (hr = spSpotSchedule.CoCreateInstance(CLSID_SpotSchedule)) return hr;
		CComVariant					v;				
		if (hr = CParameterMap::VariableArgumentListToArray(&v, 4, CComVariant(estring(szIdentifier)), CComVariant(0L), CComVariant("Last"), vSeries)) return hr;				
		if (hr = spSpotSchedule->put_Value(v)) return hr;
		spObject = spSpotSchedule;
		return S_OK;
	} else {
		CComBSTR sProgID;
		CParameterMap::ProgIDFromCLSID(clsid, sProgID);
		throw "Objects of type '" + estring(sProgID) + "' cannot not be loaded from Bwl";
	}
}
/*protected*/ HRESULT CSiriusComModule::_Load(CLSID clsid, const std::string& szIdentifier, const std::string& szDataSource, DATE date, CComPtr<IDispatch>& spObject)
{
	// We suspect the data source is of the form A (B) or B (A).
	// If A is "SiriusDirect" (and then B, if given, is a Sirius-style data source) OR
	// A is a Sirius-style data source and B is blank, then call _LoadDirect.

	// This function farms out to only _LoadFromServer, _LoadFromTaurus or _LoadDirect.
			
	bool								bDirect = false;				// true if we always call _LoadDirect
	DataSourceEnum						dsDirect = NoDataSource;		// data source to use if we call _LoadDirect
	bool								bSiriusDataSourceValid = true;	// true if we szDataSource is valid in the context of a call to _LoadDirectcould make a call to _LoadDirect withszDataSource is invalid for a _LoadDirect call
	
	if (szDataSource.size()){
		// Implement a no-questions-asked _LoadFromTaurus if the data source is set to "Taurus"
		// (Very useful for single stock market data templates).
		if (!estring::CompareNoCaseAndSpace(szDataSource, "Taurus")){
			return _LoadFromTaurus(clsid, szIdentifier, szDataSource, date, spObject);
		}
	
		// Strings that will pass bDirect are: SiriusDirect (DataSourceEnum)
		//									   DataSourceEnum (SiriusDirect).
		//
		// Strings that will pass bSiriusDataSourceValid are: Sirius (DataSourceEnum)
		//													  SiriusDirect (DataSourceEnum)
		//													  DataSourceEnum
		//													  DataSourceEnum (Sirius)
		//												      DataSourceEnum (SiriusDirect)
		//													  Sirius
		//													  SiriusDirect
		estring							szA, szB;
		estring::SplitBracket(szDataSource, &szA, &szB);
		if (!szA.size()){
			szA.assign(szB);
			szB.clear();
		}
		if (szA.size()){
			if (!CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, szA, &dsDirect)){
				// Of the form DataSourceEnum (B) or DataSourceEnum
				if (!szB.CompareNoCase("siriusdirect")) bDirect = true;
				if (szB.size() && szB.CompareNoCase("siriusdirect") && szB.CompareNoCase("sirius")) bSiriusDataSourceValid = false;
			} else if (!CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, szB, &dsDirect)){
				// Of the form A (DataSourceEnum)
				if (!szA.CompareNoCase("siriusdirect")) bDirect = true;
				if (szA.size() && szA.CompareNoCase("siriusdirect") && szA.CompareNoCase("sirius")) bSiriusDataSourceValid = false;
			} else if (szB.size()){
				// Of the form A (B) where nether A nor B are DataSourceEnum
				bSiriusDataSourceValid = false;
			} else {
				// Of the form A where A is not a DataSourceEnum
				if (!szA.CompareNoCase("siriusdirect")) bDirect = true;
				if (szA.CompareNoCase("siriusdirect") && szA.CompareNoCase("sirius")) bSiriusDataSourceValid = false;
			}
		}
	}
	
	// Consider branching to _LoadFromServer or _LoadFromTaurus
	if (!bDirect){
		switch (m_state.m_dbm){
		case TaurusExceptPositions:
			// This is a temporary option - ToDo - remove and delete TaurusExceptPositions from enums.h
			if (clsid == CLSID_Positions || clsid == CLSID_Position || clsid == CLSID_Deal || clsid == CLSID_Deals || clsid == CLSID_Product || clsid == CLSID_Products){
				return _LoadFromServer(clsid, szIdentifier, szDataSource, date, spObject);
			} else {
				return _LoadFromTaurus(clsid, szIdentifier, szDataSource, date, spObject);
			}
		case DirectExceptPositions:		
			if (clsid == CLSID_Positions || clsid == CLSID_Position){
				// ToDo - replace with _LoadFromTaurus
				return _LoadFromServer(clsid, szIdentifier, szDataSource, date, spObject);
			} else {
				break;
			}
		case UseTaurus:
			return _LoadFromTaurus(clsid, szIdentifier, szDataSource, date, spObject);
		case UseSiriusServer:									
			return _LoadFromServer(clsid, szIdentifier, szDataSource, date, spObject);		
		case RetrieveSpotSchedulesOnly:
			// Only get spot schedules from an external source (in this case we use Bwl in direct mode)
			if (clsid != CLSID_SpotSchedule) return E_FAIL;
			break;
		case DirectDatabaseConnection:	
			// Direct mode.
			break;
		case NoDatabaseConnection:
			// No database connection will always fail
			throw "Cannot load '" + szIdentifier + "' becuase the connection mode is set to 'NoDatabaseConnection'";
		}
	}

	// _LoadDirect call from here
	if (!bSiriusDataSourceValid){
		throw "The data source '" + szDataSource + "' is invalid for a native Sirius database call";
	} 
	return _LoadDirect(clsid, szIdentifier, dsDirect, date, spObject);
}
/*protected*/ HRESULT CSiriusComModule::_LoadAsset(const std::string& szIdentifier, const std::string& szDataSource, DATE date, bool bAllowCloning, CComPtr<IDispatch>& spAsset)
{
	HRESULT								hr;
	CComPtr<IAsset>						spNaturalAsset;
	CurrencyEnum						ceNatural = NoCurrency, ceComposite = NoCurrency, cePay = NoCurrency;
	CComPtr<IAsset>						spNewAsset;	
	CComPtr<IAssets>					spMaintainedAssets;				// maintained assets collection
	CAssets*							pMaintainedAssets;
	std::string							szException;
	DataSourceEnum						ds = NoDataSource;	

	if (CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, szDataSource, &ds)){		
		ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();		
	}
						
	try {
		if (!_Load(CLSID_Asset, szIdentifier, szDataSource, date, spAsset)) return S_OK;
	} catch (const std::string& sz){
		szException.assign(sz);		
	} catch (const CComBSTR& s){
		szException = estring(s);
	}

	
	// Check if szIdentifier of the form A (B) > C, A (B) or A > C.
	// attempt to construct the asset on the fly.			
	estring szNaturalAsset = estring::left(&szIdentifier, "(").left(">");			szNaturalAsset.trim();
	estring szCompositeCcy = estring::mid(&szIdentifier, "(", ")");					szCompositeCcy.trim();
	estring szPayCcy = estring::mid(&szIdentifier, ">");							szPayCcy.trim();
	if (!szNaturalAsset.size() || !szCompositeCcy.size() &&  !szPayCcy.size()) throw szException;
	
	// First check that we don't already have the asset (we may have misspelt it!)
	estring szProperName(szNaturalAsset);
	if (szCompositeCcy.size()) szProperName += " (" + szCompositeCcy + ")";
	if (szPayCcy.size()) szProperName += " > " + szPayCcy;	
	try {
		if (!_Load(CLSID_Asset, szProperName, szDataSource, date, spAsset)) return S_OK;
	} catch (...){
		// do nothing
	}
	
	if (g_pApplication->GetAssets(spMaintainedAssets)) propagate_error;
	pMaintainedAssets = dynamic_cast<CAssets*>(spMaintainedAssets.p);

	// Since the asset is of the form A (B) > C, A (B) or A > C we should be able to construct it on the fly.		
	// We first check to see if the asset with name szNaturalAsset is already in the assets collection, and,
	// if so we use that. This initially might seem inconsistent, but it mirrors the behaviour of loading a
	// basket when the constituents have already been loaded.	
	if (ds != NoDataSource){
		if (pMaintainedAssets->IsInCollection(CComObjectCollectionSerialisableKey(szNaturalAsset, date, ds))){
			if (hr = spMaintainedAssets->get_Item(CComObjectCollectionSerialisableKey(szNaturalAsset, date, ds), &spNaturalAsset)){
				// Should never happen since we have just established that the item is in the collection!
				ATLASSERT(false);
			}
		}
	}
			
	// Now attempt to load the natural asset.
	if (!spNaturalAsset){
		CComPtr<IDispatch> spObject;
		try {
			if (!_Load(CLSID_Asset, szNaturalAsset, szDataSource, date, spObject)){
				spNaturalAsset = dynamic_cast<IAsset*>(spObject.p);
				//
				// Unless we are specificly allowing cloning, we don't add baskets to the maintained collection
				// at this point. To do so causes problems. For example, imagine we have a basket A (B) > C and
				// we attempt a get_Item("A (B) > D"). This function will be called and the natural basket A
				// will be loaded. If this is added to the maintained collection then the A (B) > D may be
				// cloned from A rather than from A (B) > C. This would be unexpected.
				//
				// The problem relates specificly to baskets since a basket load A will generate an asset with
				// name A (B) > C or A > B; i.e. the returned name differs from the load ID.
				//
				VARIANT_BOOL bIsBasket;
				if (bAllowCloning || !spNaturalAsset->IsBasket(&bIsBasket) && !bIsBasket){
					if (hr = spMaintainedAssets->Add(CComVariant(), spNaturalAsset)) return hr;
				}
			}
		} catch (...){
			// do nothing
		}
	}		
	
	// If all else fails then attempt to load an asset that looks like the natural asset.
	// We only attempt this if cloning is enabled.
	// Furthermore, we only allow an asset to be created in this way if it is a basket
	// (which, of course, we don't know in advance).
	if (!spNaturalAsset && bAllowCloning){
		CComPtr<IDispatch> spObject;
		try {
			if (!_Load(CLSID_Asset, szNaturalAsset + "*", szDataSource, date, spObject)){				
				spNaturalAsset = dynamic_cast<IAsset*>(spObject.p);
				if (hr = spMaintainedAssets->Add(CComVariant(), spNaturalAsset)) return hr;
				VARIANT_BOOL bIsBasket = VARIANT_FALSE;
				spNaturalAsset->IsBasket(&bIsBasket);
				if (bIsBasket){
					// We have to throw away spNaturalAsset (but we can add it to the
					// maintained set for later use).
					spNaturalAsset = NULL;
				}
			}
		} catch (...){
			// do nothing
		}
	}				
	if (!spNaturalAsset) throw szException;

	// Check the overspecified input asset case (e.g. .FTSE (GBP) > GBP is overspelt GBP)	
	spNaturalAsset->get_Currency(&ceNatural);	
	if (szCompositeCcy.size() && (hr = CEnumMap::GetEnum("CurrencyEnum", LIBID_Sirius, szCompositeCcy, &ceComposite))) return hr;
	if (szPayCcy.size() && (hr = CEnumMap::GetEnum("CurrencyEnum", LIBID_Sirius, szPayCcy, &cePay))) return hr;
	if ((!ceComposite || ceNatural == ceComposite) && (!cePay || ceNatural == cePay)){
		// over specified asset - the natural asset is sufficient
		spAsset = spNaturalAsset;
		return S_OK;
	}
		
	// Create the new asset from spNaturalAsset	
	if (!bAllowCloning) throw szException;
	CAsset::CloneAsset(spNaturalAsset, ceComposite, cePay, spNewAsset);	
	spAsset = spNewAsset;	
	return S_OK;
}
/*protected*/ HRESULT CSiriusComModule::_LoadDirect(CLSID clsid, const std::string& szIdentifier, DataSourceEnum ds, DATE date, CComPtr<IDispatch>& spObject) const
{
	pin_date(ds, date);

	#define begin_implement_load			if (false){																\
												ATLASSERT(false);
										
	#define implement_load(Object)			} else if (clsid == CLSID_##Object){									\
												CComPtr<I##Object>	sp;												\
												HRESULT			    hr;												\
												if (hr = C##Object::Load(szIdentifier, ds, date, sp)) return hr;	\
												spObject = sp;

	#define end_implement_load				} else {																							\
												std::string			szObjectName;																\
												if (g_pApplication->GetObjectManager().CLSIDToName(clsid, &szObjectName)){						\
													g_pApplication->GetObjectManager().CollectionCLSIDToCollectionName(clsid, &szObjectName);	\
												}																								\
												return CParameterMap::ReturnErrorRS(IDS_LOAD_OBJECT, szObjectName);								\
											}																									\
											return S_OK;

	check_load_enabled	
	begin_implement_load
	implement_load(Asset)
	implement_load(Assets)
	implement_load(CorrelationMatrix)
	implement_load(Deal)
	implement_load(Deals)
	implement_load(DividendSchedule)
	implement_load(DividendSchedules)
	implement_load(Position)
	implement_load(Positions)
	implement_load(SpotSchedule)
	implement_load(SpotSchedules)
	implement_load(VolatilityStructure)
	implement_load(VolatilityStructures)
	implement_load(ZeroCurve)
	implement_load(ZeroCurves)
	end_implement_load

	#undef begin_implement_load
	#undef implement_load
	#undef end_implement_load
}

/*protected*/ HRESULT CSiriusComModule::_LoadFromServer(CLSID clsid, const std::string& szIdentifier, const std::string& szDataSource, DATE date, CComPtr<IDispatch>& spObject)
{	
	pin_date(szDataSource, date);
	
	HRESULT								hr;
	CComBSTR							sProgID;	
	CComVariant							av[4];
	static bool							bDontAttemptReattach = false;			// set to true if we don't attempt to reattach the server if the pointer is invalidated
	bool								bUseVariant;
	CComBSTR							sMethod;								// either 'LoadVariant' or 'LoadXML'
	CComVariant							vRet;									// return of method call
				
	if (!m_spSiriusServer){
		if (hr = TaurusAttach(false)) throw "Sirius encountered an error when attempting to attach the Sirius Server";
		if (!m_spSiriusServer) return E_POINTER;
	}
	CParameterMap::ProgIDFromCLSID(clsid, sProgID);	
	av[3] = sProgID;
	av[2] = (CComVariant)estring(szIdentifier);
	av[1] = (CComVariant)estring(szDataSource);
	av[0] = date;

	if (clsid == CLSID_CorrelationMatrix || clsid == CLSID_SpotSchedule || clsid == CLSID_ZeroCurve){
		bUseVariant = true;		// faster
		sMethod = L"LoadVariant";
	} else {
		bUseVariant = false;	// slower but unavoidable if the object contains pointers to other objects
		sMethod = L"LoadXML";
	}
	
	if (hr = CComDispatchDriverEx(m_spSiriusServer).InvokeN(sMethod, av, 4, &vRet)){		
		long nCode = HRESULT_CODE(hr);
		if (nCode == RPC_S_SERVER_UNAVAILABLE){
			if (bDontAttemptReattach){
				throw "The external sirius pointer has been invalidated";
			} else {
				// Attempt to reattach the server and retry
				if (TaurusAttach(false)) bDontAttemptReattach = true;
				hr = CComDispatchDriverEx(m_spSiriusServer).InvokeN(sMethod, av, 4, &vRet);
			}
		}
	}
	if (hr) propagate_error_ex(hr);
	
	if (bUseVariant){
		CComPtr<IDispatch>			sp;
		sp.CoCreateInstance(clsid);		
		if (hr = CComDispatchDriverEx(sp).PutPropertyByName(L"Value", &vRet)) return hr;
		spObject = sp;
		return S_OK; 
	} else {		
		CComVariant					v;
		if (hr = CXmlStreamer::GetVariant(estring(vRet).data(), v)) return hr;
		if (v.vt != VT_DISPATCH) return E_FAIL;
		spObject = v.pdispVal;
		return spObject ? S_OK : E_FAIL;		
	}			
}

/*protected*/ HRESULT CSiriusComModule::_LoadFromTaurus(CLSID clsid, const std::string& szIdentifier, const std::string& szDataSource, DATE date, CComPtr<IDispatch>& spObject)
{		
	pin_date(szDataSource, date);

	HRESULT								hr;
	CComBSTR							sProgID;
	CComVariant							av[5];
	static bool							bDontAttemptReattach = false;			// Set to true if we don't attempt to reattach the server if the pointer is invalidated.
	CComVariant							vSiriusItem;
	CComVariant							vObject;
		
	if (!m_spTaurus){
		if (TaurusAttach(true)) throw "Sirius encountered an error when attempting to attach Taurus";
		if (!m_spTaurus) return E_POINTER;
	}

	CParameterMap::ProgIDFromCLSID(clsid, sProgID);

	av[4] = estring(_Module.GetLocation());
	av[3] = sProgID;
	av[2] = estring(szIdentifier);
	av[1] = estring(szDataSource);
	av[0] = date;
	
	if (hr = CComDispatchDriverEx(m_spTaurus).InvokeN(L"Load", av, 5, &vSiriusItem)){
		long nCode = HRESULT_CODE(hr);
		if (nCode == RPC_S_SERVER_UNAVAILABLE){
			if (bDontAttemptReattach){
				throw "The external Taurus pointer has been invalidated";
			} else {
				// Attempt to reattach the server and retry
				if (TaurusAttach(true)) bDontAttemptReattach = true;
				hr = CComDispatchDriverEx(m_spTaurus).InvokeN(L"Load", av, 5, &vSiriusItem);
			}
		}
	}
	if (hr) propagate_error;			
	if (CComDispatchDriverEx(vSiriusItem.pdispVal).Invoke0(L"GetObject", &vObject)) propagate_error;	
	if (vObject.vt != VT_DISPATCH || !(spObject = vObject.pdispVal)) throw "No data found";
	return S_OK;
}


/////////////////////////////////////////////////////////////////////////////
//	State variable get / put
//
//	BwlInterpolationRule
//
std::string	CSiriusComModule::GetBwlInterpolationRule(void) const
{
	return m_state.m_szBwlInterpolationRule;
}
void CSiriusComModule::PutBwlInterpolationRule(const std::string& sz)
{
	m_state.m_szBwlInterpolationRule = sz;
}
//
//	BwlPassword
//
std::string CSiriusComModule::GetBwlPassword(void) const
{
	return m_state.m_szBwlPassword;
}
void CSiriusComModule::PutBwlPassword(const std::string& sz)
{
	m_state.m_szBwlPassword = sz;	
}
//
//	BwlSpotScheduleSource
//
void CSiriusComModule::GetBwlSpotScheduleSource(std::vector<std::string>* pasz) const
{
	estring::Split(GetBwlSpotScheduleSource(), ",", pasz);
}
std::string CSiriusComModule::GetBwlSpotScheduleSource(void) const
{
	return m_state.m_szBwlSpotScheduleSource;
}
void CSiriusComModule::PutBwlSpotScheduleSource(const std::string& sz)
{
	m_state.m_szBwlSpotScheduleSource = sz;
}
//
//	BwlSQLServerDSN
//
std::string CSiriusComModule::GetBwlSQLServerDSN(void) const
{
	return m_state.m_szBwlSQLServerDSN;
}
void CSiriusComModule::PutBwlSQLServerDSN(const std::string& sz)
{
	m_state.m_szBwlSQLServerDSN = sz;	
	ResetBwlDatabaseConnection();
}
//
//	BwlUserName
//
std::string CSiriusComModule::GetBwlUserName(void) const
{
	return m_state.m_szBwlUserName;	
}
void CSiriusComModule::PutBwlUserName(const std::string& sz)
{
	m_state.m_szBwlUserName = sz;
}
//
//	DatabaseMode
//
DatabaseModeEnum CSiriusComModule::GetDatabaseMode(void) const
{	
	return m_state.m_dbm;	
}
void CSiriusComModule::PutDatabaseMode(DatabaseModeEnum dbm)
{
	m_state.m_dbm = dbm;
	ResetAllConnections();
	CExcelInterface::UpdateSiriusStatusCaption();
}
//
//	DefaultDataSource
//
DataSourceEnum CSiriusComModule::GetDefaultDataSource(void) const
{
	return m_state.m_dsDefault;
}
void CSiriusComModule::PutDefaultDataSource(DataSourceEnum dsDefault)
{
	m_state.m_dsDefault = dsDefault;
	m_state.m_szDataSourceDefault = CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, dsDefault);
}
const std::string& CSiriusComModule::GetDefaultDataSourceStr(void) const
{
	return m_state.m_szDataSourceDefault;
}
//
//	DisplayProductOnly
//
bool CSiriusComModule::GetDisplayProductOnly(void) const
{
	return m_state.m_bDisplayProductOnly;
}
void CSiriusComModule::PutDisplayProductOnly(bool b)
{
	m_state.m_bDisplayProductOnly = b;
}
//
//	DisplaySiriusStatus
//
bool CSiriusComModule::GetDisplaySiriusStatus(void) const
{
	return m_state.m_bDisplaySiriusStatus;	
}
void CSiriusComModule::PutDisplaySiriusStatus(bool b)
{
	m_state.m_bDisplaySiriusStatus = b;	
	CExcelInterface::UpdateSiriusStatusCaption();
}
//
//	EnablePublishing
//
bool CSiriusComModule::GetEnablePublishing(void) const
{
	return m_state.m_bEnablePublishing;
}
void CSiriusComModule::PutEnablePublishing(bool b)
{
	m_state.m_bEnablePublishing = b;
}
//
//	FileSystem
//
FileSystemEnum CSiriusComModule::GetFileSystem(void) const
{
	return m_state.m_fs;
}
void CSiriusComModule::PutFileSystem(FileSystemEnum fs)
{	
	m_state.m_fs = fs;	
	CExcelInterface::UpdateSiriusStatusCaption();	
}
//
//	FileSystemRoot
//
std::string CSiriusComModule::GetFileSystemRoot(void) const
{
	return m_state.m_szFileSystemRoot;
}
void CSiriusComModule::PutFileSystemRoot(const std::string& sz)
{
	m_state.m_szFileSystemRoot = sz;	
	if (m_state.m_szFileSystemRoot.size() && m_state.m_szFileSystemRoot[m_state.m_szFileSystemRoot.size() - 1] == '\\'){
		m_state.m_szFileSystemRoot.resize(m_state.m_szFileSystemRoot.size() - 1);
	}
}
//
//	Location
//
std::string CSiriusComModule::GetLocation(void) const
{
	return m_state.m_szLocation;
}
void CSiriusComModule::PutLocation(std::string* psz, const connection_parameters* pcp, FileSystemEnum* pfs)
{		
	psz->assign(m_state.m_szLocation = ValidateLocation(*psz, pcp, pfs));
}
//
//	LocationFirstTry
//
std::string CSiriusComModule::GetLocationFirstTry(void) const
{
	return m_state.m_szLocationFirstTry;
}
void CSiriusComModule::PutLocationFirstTry(std::string* psz, const connection_parameters* pcp, FileSystemEnum* pfs)
{
	psz->assign(m_state.m_szLocationFirstTry = ValidateLocation(*psz, pcp, pfs));
}
//	
//	LocationFirstTryEnabled
//
bool CSiriusComModule::GetLocationFirstTryEnabled(void) const
{
	return m_state.m_bLocationFirstTryEnabled;
}
void CSiriusComModule::PutLocationFirstTryEnabled(bool b)
{
	m_state.m_bLocationFirstTryEnabled = b;
}
//
//	ProductTypesFileName
//
std::string CSiriusComModule::GetProductTypesFileName(void) const
{
	return m_state.m_szProductTypesFileName;
}
void CSiriusComModule::PutProductTypesFileName(const std::string& sz)
{
	m_state.m_szProductTypesFileName = sz;		
}
//
//	SiriusPassword
//
std::string CSiriusComModule::GetSiriusPassword(void) const
{
	return m_state.m_szSiriusPassword;	
}
void CSiriusComModule::PutSiriusPassword(const std::string& sz)
{
	m_state.m_szSiriusPassword = sz;
}
//
//	SiriusSQLServerDSN
//
std::string CSiriusComModule::GetSiriusSQLServerDSN(void) const
{
	return m_state.m_szSiriusSQLServerDSN;	
}
void CSiriusComModule::PutSiriusSQLServerDSN(const std::string& sz)
{
	m_state.m_szSiriusSQLServerDSN = sz;
	ResetSiriusDatabaseConnection();
	CExcelInterface::UpdateSiriusStatusCaption();
}
//
//	SiriusUserName
//
std::string CSiriusComModule::GetSiriusUserName(void) const
{
	return m_state.m_szSiriusUserName;	
}
void CSiriusComModule::PutSiriusUserName(const std::string& sz)
{
	m_state.m_szSiriusUserName = sz;
}
//
//	UseBwlForSpotSchedule
//
bool CSiriusComModule::GetUseBwlForSpotSchedule(void) const
{
	return m_state.m_bUseBwlForSpotSchedule;
}
void CSiriusComModule::PutUseBwlForSpotSchedule(bool b)
{
	m_state.m_bUseBwlForSpotSchedule = b;
}
//
//	UseYesterdaysClose
//
/*static*/ bool CSiriusComModule::GetUseYesterdaysClose(void)
{
	return _Module.m_state.m_bUseYesterdaysClose;
}
void CSiriusComModule::PutUseYesterdaysClose(bool b)
{
	m_state.m_bUseYesterdaysClose = b;
}
//
//	Read state from registry
//
void CSiriusComModule::RegistryToMembers(CSiriusComModule::state* pstate, bool bRefresh) const
{
	estring								szValue;
	
	pstate->m_szBwlInterpolationRule = "";	
	RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_bwl_interpolation, &pstate->m_szBwlInterpolationRule, false);

	pstate->m_szBwlPassword = "";
	RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_bwl_password, &pstate->m_szBwlPassword, true);	

	pstate->m_szBwlSpotScheduleSource = "";
	RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_bwl_spot_schedule_source, &pstate->m_szBwlSpotScheduleSource, false);	

	pstate->m_szBwlSQLServerDSN = "";	
	RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_bwl_sql_server, &pstate->m_szBwlSQLServerDSN, false);	
		
	pstate->m_szBwlUserName = "";
	RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_bwl_user_name, &pstate->m_szBwlUserName, false);
	
	pstate->m_dbm = NoDatabaseConnection;		
	if (!RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_database_mode, &szValue, false)){
		CEnumMap::GetEnum("DatabaseModeEnum", LIBID_Sirius, szValue, &pstate->m_dbm);
	}
	
	pstate->m_dsDefault = Last;
	if (!RegistryToMember(HKEY_CURRENT_USER, s_sz_default_data_source, &szValue, false)){
		CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, szValue, &pstate->m_dsDefault);
	}	
	pstate->m_szDataSourceDefault = CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, pstate->m_dsDefault);
	
	pstate->m_bDisplayProductOnly = false;
	if (!RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_product_only, &szValue, false)){
		pstate->m_bDisplayProductOnly = !szValue.CompareNoCase("true") ? true : false;
	}		
	
	pstate->m_bDisplaySiriusStatus = false;	
	if (!RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_display_sirius_status, &szValue, false)){
		pstate->m_bDisplaySiriusStatus = !szValue.CompareNoCase("true") ? true : false;
	}

	pstate->m_bEnablePublishing = false;
	if (!RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_enable_publishing, &szValue, false)){
		pstate->m_bEnablePublishing = !szValue.CompareNoCase("true") ? true : false;
	}

	pstate->m_fs = fsUnknown;
	if (!RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_file_system, &szValue, false)){
		CEnumMap::GetEnum("FileSystemEnum", LIBID_Sirius, szValue, &pstate->m_fs);
	}	

	pstate->m_szFileSystemRoot = "";
	RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_file_root, &pstate->m_szFileSystemRoot, false);
	if (pstate->m_szFileSystemRoot.size() && pstate->m_szFileSystemRoot[pstate->m_szFileSystemRoot.size() - 1] == '\\'){
		pstate->m_szFileSystemRoot.resize(pstate->m_szFileSystemRoot.size() - 1);
	}
	
	pstate->m_szLocation = "";
	RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_location, &pstate->m_szLocation, false);

	pstate->m_szLocationFirstTry = "";
	RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_location_first_try, &pstate->m_szLocationFirstTry, false);

	pstate->m_bLocationFirstTryEnabled = false;
	if (!RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_location_first_try_enabled, &szValue, false)){
		pstate->m_bLocationFirstTryEnabled = !szValue.CompareNoCase("true") ? true : false;
	}		

	pstate->m_szProductTypesFileName = "";
	RegistryToMember(HKEY_LOCAL_MACHINE, s_sz_registry_product_types, &pstate->m_szProductTypesFileName, false);
	
	pstate->m_szSiriusPassword = "";
	RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_sirius_password, &pstate->m_szSiriusPassword, true);	

	pstate->m_szSiriusSQLServerDSN = "";	
	RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_sirius_sql_server, &pstate->m_szSiriusSQLServerDSN, false);	
		
	pstate->m_szSiriusUserName = "";
	RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_sirius_user_name, &pstate->m_szSiriusUserName, false);
			
	pstate->m_bUseBwlForSpotSchedule = false;
	if (!RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_use_bwl_for_spot_schedule, &szValue, false)){
		pstate->m_bUseBwlForSpotSchedule = !szValue.CompareNoCase("true") ? true : false;
	}
	
	pstate->m_bUseYesterdaysClose = false;
	if (!RegistryToMember(HKEY_CURRENT_USER, s_sz_registry_use_yesterday_close, &szValue, false)){
		pstate->m_bUseYesterdaysClose = !szValue.CompareNoCase("true") ? true : false;
	}

	if (bRefresh){
		_Module.ResetAllConnections();
		CExcelInterface::UpdateSiriusStatusCaption();
	}
}
//
//	Save state to registry
//
void CSiriusComModule::MembersToRegistry(void)
{	
	if (m_state.m_szBwlInterpolationRule != m_state_init.m_szBwlInterpolationRule){
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_bwl_interpolation, m_state.m_szBwlInterpolationRule, false, false);
	}

	if (m_state.m_szBwlPassword != m_state_init.m_szBwlPassword){		
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_bwl_password, m_state.m_szBwlPassword, false, true);
	}
			
	if (m_state.m_szBwlSQLServerDSN != m_state_init.m_szBwlSQLServerDSN){	
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_bwl_sql_server, m_state.m_szBwlSQLServerDSN, false, false);
	}
	
	if (m_state.m_szBwlUserName != m_state_init.m_szBwlUserName){
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_bwl_user_name, m_state.m_szBwlUserName, false, false);
	}

	if (m_state.m_dbm != m_state_init.m_dbm){
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_database_mode, CEnumMap::GetString("DatabaseModeEnum", LIBID_Sirius, m_state.m_dbm), false, false);
	}
	
	if (m_state.m_dsDefault != m_state_init.m_dsDefault){
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_default_data_source, CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_state.m_dsDefault), false, false);
	}

	if (m_state.m_bDisplayProductOnly != m_state_init.m_bDisplayProductOnly){
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_product_only, m_state.m_bDisplayProductOnly ? "True" : "False", false, false);
	}

	if (m_state.m_bDisplaySiriusStatus != m_state_init.m_bDisplaySiriusStatus){
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_display_sirius_status, m_state.m_bDisplaySiriusStatus ? "True" : "False", false, false);
	}

	if (m_state.m_bEnablePublishing != m_state_init.m_bEnablePublishing){
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_enable_publishing, m_state.m_bEnablePublishing ? "True" : "False", false, false);
	}

	if (m_state.m_fs != m_state_init.m_fs){
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_file_system, CEnumMap::GetString("FileSystemEnum", LIBID_Sirius, m_state.m_fs), false, false);
	}

	if (m_state.m_szFileSystemRoot != m_state_init.m_szFileSystemRoot){
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_file_root, m_state.m_szFileSystemRoot, false, false);
	}

	if (m_state.m_szLocation != m_state_init.m_szLocation){		
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_location, m_state.m_szLocation, false, false);
	}

	if (m_state.m_szLocationFirstTry != m_state_init.m_szLocationFirstTry){
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_location_first_try, m_state.m_szLocationFirstTry, false, false);
	}

	if (m_state.m_bLocationFirstTryEnabled != m_state_init.m_bLocationFirstTryEnabled){
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_location_first_try_enabled, m_state.m_bLocationFirstTryEnabled ? "True" : "False", false, false);
	}

	if (m_state.m_szProductTypesFileName != m_state_init.m_szProductTypesFileName){
		MemberToRegistry(HKEY_LOCAL_MACHINE, s_sz_registry_product_types, m_state.m_szProductTypesFileName, false, false);
	}

	if (m_state.m_szSiriusPassword != m_state_init.m_szSiriusPassword){
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_sirius_password, m_state.m_szSiriusPassword, false, false);
	}

	if (m_state.m_szBwlSpotScheduleSource != m_state_init.m_szBwlSpotScheduleSource){	
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_bwl_spot_schedule_source, m_state.m_szBwlSpotScheduleSource , false, false);
	}
	
	if (m_state.m_szSiriusSQLServerDSN != m_state_init.m_szSiriusSQLServerDSN){	
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_sirius_sql_server, m_state.m_szSiriusSQLServerDSN, false, false);
	}

	if (m_state.m_szSiriusUserName != m_state_init.m_szSiriusUserName){
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_sirius_user_name, m_state.m_szSiriusUserName, false, false);
	}

	if (m_state.m_bUseBwlForSpotSchedule != m_state_init.m_bUseBwlForSpotSchedule){
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_use_bwl_for_spot_schedule, m_state.m_bUseBwlForSpotSchedule ? "True" : "False", false, false);
	}

	if (m_state.m_bUseYesterdaysClose != m_state_init.m_bUseYesterdaysClose){		
		MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_use_yesterday_close, m_state.m_bUseYesterdaysClose ? "True" : "False", false, false);
	}
}


/////////////////////////////////////////////////////////////////////////////
//	All other implementation functions
//
HRESULT CSiriusComModule::AddError(const estring& sz) const
{
	ATLASSERT(GetSiriusApplication());	
	return GetSiriusApplication()->AddError((CComBSTR)sz);
	// ToDo - This is where we could fire an error event on the Sirius Application.
	// Perhaps we could then get rid of AddError etc.
}

// This encryption routing routine uses an algorithm created by Rodney Thayer. 
// This algorithm can use keys of various sizes. With a 40-bit key (5 characters) 
// it can be freely exported from the U.S. The cipher is considered robust with 
// 128 bits of key material but can use up to 2048 bits. It is compatible with RSA’s RC4 algorithm.
/*static*/ void CSiriusComModule::Crypt(TCHAR *inp, DWORD inplen, const TCHAR* key /*= ""*/, DWORD keylen /*= 0*/)
{
    // we will consider size of sbox 256 bytes
    // (extra bytes are only to prevent any mishaps just in case)
    TCHAR Sbox[257], Sbox2[257];
    unsigned long i, j, t, x;

    // this unsecured key is to be used only when there is no input key from user
    static const TCHAR  OurUnSecuredKey[] = "Sirius" ;
    static const int OurKeyLen = std::strlen(OurUnSecuredKey);    
    TCHAR temp , k;
    i = j = k = t =  x = 0;
    temp = 0;

    // always initialize the arrays with zero
	std::memset(Sbox, 0, sizeof(Sbox));
	std::memset(Sbox2, 0, sizeof(Sbox2));

    // initialize sbox i
    for (i = 0; i < 256U; i++){
        Sbox[i] = (TCHAR)i;
    }

    j = 0;
    // whether user has sent any inpur key
    if (keylen){
        // initialize the sbox2 with user key
        for (i = 0; i < 256U ; i++){
            if (j == keylen) j = 0;                
            Sbox2[i] = key[j++];
        }    
    } else {
        // initialize the sbox2 with our key
        for (i = 0; i < 256U ; i++){
            if (j == OurKeyLen) j = 0;            
            Sbox2[i] = OurUnSecuredKey[j++];
        }
    }

    j = 0; // Initialize j
    // scramble sbox1 with sbox2
    for (i = 0; i < 256; i++){
        j = (j + (unsigned long) Sbox[i] + (unsigned long) Sbox2[i]) % 256U ;
        temp =  Sbox[i];                    
        Sbox[i] = Sbox[j];
        Sbox[j] =  temp;
    }

    i = j = 0;
    for(x = 0; x < inplen; x++){
        // increment i
        i = (i + 1U) % 256U;
        // increment j
        j = (j + (unsigned long) Sbox[i]) % 256U;

        // Scramble SBox #1 further so encryption routine will
        // will repeat itself at great interval
        temp = Sbox[i];
        Sbox[i] = Sbox[j] ;
        Sbox[j] = temp;

        // Get ready to create pseudo random  byte for encryption key
        t = ((unsigned long) Sbox[i] + (unsigned long) Sbox[j]) %  256U ;

        // get the random byte
        k = Sbox[t];

        // xor with the data and done
        inp[x] = (inp[x] ^ k);
    }    
}

void CSiriusComModule::CheckBwlConnection(void)
{
	if (!m_spBwlEngine){
		m_spBwlEngine.CoCreateInstance(CLSID_BwlEngine);
		if (!m_spBwlEngine) m_BwlUsable = Unusable;
	}		
	if (m_BwlUsable == Unknown){
		// connect to bwl		
		estring		szSystem = GetBwlSQLServerDSN();
		estring		szUser = GetBwlUserName();
		estring		szPassword = GetBwlPassword();
		CComVariant	av[] = {szSystem.GetValue(), szUser.GetValue(), szPassword.GetValue(), VARIANT_TRUE};
		m_BwlUsable = CComDispatchDriverEx(m_spBwlEngine).InvokeN(L"Open", av, 4) ? Unusable : Usable;
	}
	if (m_BwlUsable == Unusable) throw "Bwl not usable";
}

HRESULT CSiriusComModule::DeclareExcel(const VARIANT& ExcelApplicationInstance)
{		
	IID									iid;
	HRESULT								hr;	
	
	if (ExcelApplicationInstance.vt != VT_DISPATCH) return E_FAIL;
	m_spExcel = ExcelApplicationInstance.pdispVal;
	if (!(hr = CParameterMap::GetObjectIID(m_spExcel, &iid))){
		if (iid != IID_ExcelApplication) hr = E_FAIL;
	}	
	if (hr){
		m_spExcel = NULL;
	} else {		
		m_hImageList = ImageList_LoadBitmap(_pModule->m_hInstResource, (LPCTSTR)IDB_IMAGELIST, 16, 16, RGB(255, 255, 255));
		m_hFont = ::CreateFont(-13, 0, 0, 0, 400, FALSE, FALSE, 0, ANSI_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH | FF_MODERN, _T("Courier New"));	
		INITCOMMONCONTROLSEX	iccx;
		std::memset(&iccx, NULL, sizeof(iccx));
		iccx.dwSize = sizeof(iccx);
		iccx.dwICC = /*ICC_COOL_CLASSES | ICC_TREEVIEW_CLASSES |*/ ICC_USEREX_CLASSES;
		if (!::InitCommonControlsEx(&iccx)) ATLASSERT(false);	
	}
	return hr;
}

HRESULT CSiriusComModule::GetConnection(_ConnectionPtr& spConnection, const connection_parameters* pcp/* = NULL*/) const
// pcp - If null then we use member variables of this class to establish a database
//		 connecction and m_spConnection is modified.
//		 If not null then we use pcp members to establish a connection. m_spConnection
//		 if not modified.
{
	HRESULT								hr;
	
	if (pcp){
		// don't adjust m_spConnection; create a new connection
		if (hr = spConnection.CreateInstance(__uuidof(Connection))) return hr;
	} else {
		if (m_spConnection == NULL){
			// Create m_spConnection if possible. We create it on a just in time basis (as opposed to using a
			// the base class initialiser 'm_spConnection(NULL)' in this class' constructor) since, for example,
			// in a VB-exe that references this dll may attempt to instantiate the global CSiriusComModule object
			// before VB has initialised COM (via CoInitialise).
			if (hr = spConnection.CreateInstance(__uuidof(Connection))) return hr;
			m_spConnection = spConnection;
		} else {
			spConnection = m_spConnection;
		}
		if (spConnection == NULL) return E_FAIL;
	}

	if (!(spConnection->State & adStateOpen)){
		// we need to establish a connection to a Sirius database		
		std::string		szDSN = pcp ? pcp->m_szDSN : GetSiriusSQLServerDSN();
		bool			bNTAuthenticated = GetNTAuthentication(szDSN);
		std::string		szConnect;
		
		if (bNTAuthenticated){		
			szConnect = "DSN=" + szDSN + ";Trusted_Connection=yes";
		} else {			
			szConnect = "DSN=" + szDSN + ";UID=" + (pcp ? pcp->m_szUserName : GetSiriusUserName()) + ";PWD=" + (pcp ? pcp->m_szPassword : GetSiriusPassword());
		}		
		spConnection->Open(_bstr_t(szConnect.c_str()), _bstr_t(""), _bstr_t(""), -1);		
	}

	if (spConnection->State & adStateOpen){
		return S_OK;
	} else {
		CParameterMap::ThrowComErrorR(IDS_UNHANDLED_EXCEPTION);
		return E_FAIL;
	}
}

HRESULT CSiriusComModule::DisplayExplorer(void)
{
	HRESULT								hr;
	static const CLSID					CLSID_SiriusExplorer = {0x72BD51D3,0x6047,0x4950,{0xAC,0xB6,0x14,0x1D,0x63,0x3E,0x21,0xEB}};
	CComVariant							vObject;
	
	if (!m_spExplorer && (hr = m_spExplorer.CoCreateInstance(CLSID_SiriusExplorer))) return hr;
	if (hr = CComDispatchDriverEx(m_spExplorer).Invoke0(L"GetObject", &vObject)){
		// The pointer may have been invalidated.
		long nCode = HRESULT_CODE(hr);
		if (nCode == RPC_S_SERVER_UNAVAILABLE){			
			m_spExplorer = NULL;
			if (hr = m_spExplorer.CoCreateInstance(CLSID_SiriusExplorer)) return hr;
		} else {
			return hr;
		}
	}
	if (vObject.vt == VT_EMPTY || vObject.vt == VT_DISPATCH && !vObject.pdispVal){
		// we need to set up the explorer
		if (hr = CComDispatchDriverEx(m_spExplorer).Invoke1(L"SetObject", &CComVariant(GetSiriusApplication().p))) return hr;
	}
	if (hr = CComDispatchDriverEx(m_spExplorer).Invoke0(L"Show")) return hr;
	return S_OK;
}

bool CSiriusComModule::GetEnableLoad(void) const
{
	return m_state.m_dbm != NoDatabaseConnection ? true : false;
}

HRESULT CSiriusComModule::GetExcel(CComPtr<IDispatch>& spExcel) const
{
	spExcel = m_spExcel;
	return spExcel ? S_OK : E_FAIL;
}

HRESULT CSiriusComModule::GetExcel(CComDispatchDriverEx& ddExcel) const
{
	ddExcel = m_spExcel;
	return m_spExcel ? S_OK : E_FAIL;
}

HRESULT CSiriusComModule::GetMarketDataRoot(std::string* psz) const
{
	std::string							szFileSystem = GetFileSystemRoot();
	std::string							szDSN = GetSiriusSQLServerDSN();
	FileSystemEnum						fs;

	szFileSystem += "\\MarketData\\";
	if ((fs = GetFileSystem()) == fsNTFS){
		psz->assign(szFileSystem + "Local");
	} else if (fs == fsSQLServer){
		psz->assign(szFileSystem + szDSN);
	} else {
		throw "Invalid file system value";
	}
	return S_OK;
}

std::string CSiriusComModule::GetName(void) const
{
	static std::string					szName;
	
	if (szName.size()) return szName;

	HINSTANCE							h = _Module.GetModuleInstance();
	TCHAR								szDllName[MAX_PATH + 1];

	::GetModuleFileName(h, szDllName, sizeof(szDllName));
	return szName = estring(szDllName).right("\\");
}

//	Returns true if an input DSN has the Trusted_Connection attribute set to "yes".
bool CSiriusComModule::GetNTAuthentication(const std::string& szDSN) const
{
	DWORD								dwCount = 16;	// size of buffer, 16 will do for a value that should be set to "yes" / "no"
	CRegKey								key;
	char								szValue[MAX_PATH];
		
	if (!szDSN.size()) return false;
	if (key.Open(HKEY_LOCAL_MACHINE, ("Software\\ODBC\\ODBC.INI\\" + szDSN).c_str(), KEY_READ) == ERROR_SUCCESS){
		dwCount = MAX_PATH;
		if (key.QueryValue(szValue, _T("Trusted_Connection"), &dwCount) == ERROR_SUCCESS){
			return estring("yes").CompareNoCase(szValue) ? false : true;
		}
	}
	return false;		
}

//
//	Returns the ID of the user and the time when an object was saved 
//  to the database.
//
CComVariant CSiriusComModule::GetObjectSignature(const std::string& szObjectName, const std::string& szIdentifier, const std::string& szDataSource, long nDate) const
{	
	FileSystemEnum						fs;	
	CComVariant							vUser, vTime;
	
	if (!nDate) throw "Invalid date " + estring(nDate);
	
	// load the XML for the asset	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){		
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;						// SQL query to execute
		_RecordsetPtr					prs;							// ADO recordset		
				
		ssQuery << "sp_user_get_object_signature '" << szObjectName << "', '" << szIdentifier << "', '" << szDataSource << "', " << nDate << ", '" << GetLocation() << "'";
		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
			if (prs->adoEOF) throw "No object found";
			vUser = prs->GetFields()->GetItem(CComVariant(L"User"))->GetValue();
			vTime = prs->GetFields()->GetItem(CComVariant(L"Time"))->GetValue();			
		} catch (_com_error& e){
			throw estring(e);
		}				
	} else if (fs == fsNTFS){		
		IID							iid;				
		std::string					szFileName;
		DataSourceEnum				ds;
		WIN32_FILE_ATTRIBUTE_DATA	fad;
		SYSTEMTIME					st;
		DATE						date;
		
		if (g_pApplication->GetObjectManager().NameToIID(szObjectName, &iid)) propagate_error;		
		if (CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, szDataSource, &ds)) throw "Invalid data source '" + szDataSource + "'";
		if (g_pApplication->GetObjectManager().GetFileName(iid, GetLocation(), szIdentifier, ds, nDate, &szFileName)) propagate_error;						
		if (!::GetFileAttributesEx(szFileName.c_str(), GetFileExInfoStandard, &fad)) throw "No object found";
		if (!::FileTimeToSystemTime(&fad.ftLastAccessTime, &st) || !::SystemTimeToVariantTime(&st, &date)) throw "Unhandled exception in CSiriusComModule::GetObjectSignature";
		vTime = date;
		vTime.ChangeType(VT_DATE);
		vUser = "???";		
	} else {
		throw "Unhandled exception in CSiriusComModule::GetObjectSignature";
	}	

	CParameterMap pmRet;		
	vTime.ChangeType(VT_BSTR);
	CComVariant vUserHeader("User");
	CComVariant vTimeHeader("Time");
	pmRet.SetValues(2, 2, vUserHeader, vUser, vTimeHeader, vTime);		// Can't use anonymous temporaries here due to Microsoft limitation?
	return pmRet.GetValue();
}

// get a pointer to the singleton sirius application object; creating it if necessary
CComPtr<ISiriusApplication> CSiriusComModule::GetSiriusApplication(void) const
{
	if (!m_spApplication){
		if (m_spApplication.CoCreateInstance(CLSID_SiriusApplication)){
			CParameterMap::DisplayError(L"Fatal error encountered when creating the global Sirius Application object!\n\nYou cannot use Sirius in this state.", MB_ICONSTOP);
		}
	}
	return m_spApplication;
}
CSiriusApplication* CSiriusComModule::_GetSiriusApplication(void) const
{
	return dynamic_cast<CSiriusApplication*>(GetSiriusApplication().p);
}

HRESULT CSiriusComModule::GetSiriusCaption(CComVariant* pv) const
{
	HRESULT								hr;
	estring								sz;
	if (hr = GetSiriusCaption(&sz)) return hr;
	return sz.GetValue(pv);
}

HRESULT CSiriusComModule::GetSiriusCaption(std::string* psz) const
{
	HRESULT								hr;	
	FileSystemEnum						fs;
	estring								szVersion;
	std::stringstream					ssCaption;
	CComDispatchDriverEx				ddExcel;
	CComVariant							vName;

	if (hr = GetExcel(ddExcel)) return hr;
	if (hr = ddExcel.GetPropertyByName(L"Name", &vName)) return hr;
	if (hr = g_pApplication->GetVersion(GetName(), &szVersion)) return hr;
	if ((fs = GetFileSystem()) == fsSQLServer){
		std::stringstream		ssQuery;
		_ConnectionPtr			spConnection;
		_RecordsetPtr			prs;
		
		ssQuery << "sp_user_get_caption " << "'" << GetName() << "', '" << szVersion << "'";
		try {
			GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
			if (prs->adoEOF){
				ssCaption << "Unavailable";
			} else {
				int nCount = 0;
				while (!prs->adoEOF){				
					if (nCount++) ssCaption << ", ";
					ssCaption << estring((CComVariant)prs->GetFields()->GetItem(0L)->GetValue());
					prs->MoveNext();
				}
			}
		} catch (...){
			ssCaption << "Unavailable";
		}
	} else if (fs == fsNTFS){
		ssCaption << "Local Database";
	}

	psz->assign(estring(vName) + " (" + ssCaption.str() + ")");	
	return psz->size() ? S_OK : E_FAIL;
}

/*static*/ bool	CSiriusComModule::IsReutersRunning(void)
{
	return FindWindow("DataDictPDD35", NULL) ? true : false;
}

double CSiriusComModule::LoadBwlSpot(const std::string& szCode, const std::string& szCodeType, const std::string& szDatabase, const std::string& szDataSource, const std::string& szInterpolation)
//	szDataSource - if blank then we use the default bwl spot schedule sources
{
	CheckBwlConnection();
	HRESULT							hr;
	std::vector<std::string>		aszDataSource;
	CComVariant						vProperty;	

	if (szDatabase == "Reuters"){
		// szDataSource must be blank or set to "Reuters"
		if (szDataSource.size() && szDataSource != "Reuters"){
			throw "Inconsistent database and data source names in CSiriusComModule::LoadBwlSpot. You have the database set to '" + szDatabase + "' and the data source set to '" + szDataSource + "'";
		}
	}
		
	if (szDataSource.size()){
		aszDataSource.push_back(szDataSource);
	} else {
		GetBwlSpotScheduleSource(&aszDataSource);
	}
		
	if (!aszDataSource.size()) aszDataSource.resize(1);	// delegate the defaulting to bwl
	for (std::vector<std::string>::const_iterator it = aszDataSource.begin(); it != aszDataSource.end(); it++){	
		CComVariant av[] =				{	estring::GetValue(szInterpolation)/*MissingRuleOpt*/, 
											estring::GetValue(*it)/*DataSrcOpt*/,
											estring::GetValue(szCodeType)/*CodeTypeOpt*/,
											0L/*DateOpt*/,
											L"Close"/*PropertyName*/,
											estring::GetValue(szCode)/*Code*/,
											estring::GetValue(szDatabase)/*Database*/ };
		if (!(hr = CComDispatchDriverEx(m_spBwlEngine).InvokeN(L"LoadProperty", av, 7, &vProperty))) break;
	}
	if (hr){
		// failed for all the data sources
		propagate_error;
	}
	if (vProperty.vt != VT_R8) throw "Invalid spot value for '" + szCode + "'";	
	return vProperty.dblVal;
}

/*static*/ HRESULT CSiriusComModule::MemberToRegistry(HKEY hKey, const std::string& szKey, std::string szValue, bool bSilent, bool bEncrypt)
//	bSilent - true if we suppress error messages
{
	CRegKey								key;
	LONG								n;
	std::string							sz(szValue);
	
	key.Create(hKey, s_sz_registry_path);
	if (bEncrypt) Crypt(const_cast<char*>(szValue.data()), szValue.size(), szKey.data(), szKey.size());			
	if ((n = key.SetValue(hKey, s_sz_registry_path, szValue.c_str(), szKey.c_str())) == ERROR_SUCCESS) return S_OK;
	if (!bSilent){
		void* p;
        ::FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, NULL, n, 0, (LPTSTR)&p, 0, NULL);
		MessageBox(::GetActiveWindow(), (LPTSTR)p, NULL, MB_OK | MB_ICONSTOP);		
        ::LocalFree(p);
	}
	return E_FAIL;
}

/*static*/ HRESULT CSiriusComModule::RegistryToMember(HKEY hKey, const std::string& szName, std::string* psz, bool bEncrypted)
{
	DWORD								dwCount;
	CRegKey								key;
	char								szValue[MAX_PATH];
		
	if (key.Open(hKey, s_sz_registry_path, KEY_READ) == ERROR_SUCCESS){
		dwCount = MAX_PATH;
		if (key.QueryValue(szValue, szName.c_str(), &dwCount) == ERROR_SUCCESS){
			if (bEncrypted) Crypt(szValue, std::strlen(szValue), szName.data(), szName.size());
			psz->assign(szValue);
			estring::trim(psz);
			return S_OK;
		}
	}
	return E_FAIL;
}

void CSiriusComModule::ResetAllConnections(void)
{
	ResetSiriusDatabaseConnection();		// Sirius database
	ResetBwlDatabaseConnection();			// Bwl database
	TaurusAttach(false);					// Taurus connection
}

void CSiriusComModule::ResetBwlDatabaseConnection(void)
{	
	if (!m_spBwlEngine){	
		m_spBwlEngine.CoCreateInstance(CLSID_BwlEngine);
		if (!m_spBwlEngine){
			m_BwlUsable = Unusable;
			return;
		}
	}
	if (m_BwlUsable == Unknown) return;	// no need to do anything

	// disconnect from bwl			
	CComDispatchDriverEx(m_spBwlEngine).InvokeN(L"Close", &CComVariant(VARIANT_TRUE), 1);				
	m_BwlUsable = Unknown;
}

void CSiriusComModule::ResetSiriusDatabaseConnection(void)
{		
	try {			
		m_spConnection->Close();
	} catch (_com_error& e){
		ATLTRACE("%s\n", estring(e).c_str());
	}
}

//	Set m_state to the registry settings
void CSiriusComModule::Revert(void)
{
	RegistryToMembers(&m_state, true);
	m_state_init = m_state;
}

//	Attach this process to Taurus (bouncing any existing connection)
HRESULT CSiriusComModule::TaurusAttach(bool bAlwaysAttach)
{
	HRESULT								hr = S_OK;
		
	m_spSiriusServer = NULL;
	m_spTaurus = NULL;
	
	if (bAlwaysAttach || m_state.m_dbm == UseTaurus || m_state.m_dbm == TaurusExceptPositions){
		hr = m_spTaurus.CoCreateInstance(CLSID_Taurus);
	}

	if (m_state.m_dbm == UseSiriusServer || m_state.m_dbm == DirectExceptPositions || m_state.m_dbm == TaurusExceptPositions){
		HRESULT hrSirius;
		// Try a few times but don't worry about the fact that this is a kludge since we're getting rid of Sirius Server.
		for (long nTry = 1; nTry <= 10; nTry++){
			if (!(hrSirius = m_spSiriusServer.CoCreateInstance(CLSID_SiriusDataManager, NULL, CLSCTX_LOCAL_SERVER))) break;			
			// Wait for nTry seconds
			::Sleep(nTry * 1000);
		}
		if (hrSirius) return hrSirius;
	}

	return hr;	
}

void CSiriusComModule::Term(void)
{		
	// Write display product only member of m_state to the registry. We do this
	// since this member is implemented in the Insert Product dialog. This does
	// not have 'Save' and 'Revert' buttons like the dialog(s) correponding to
	MemberToRegistry(HKEY_CURRENT_USER, s_sz_registry_product_only, m_state.m_bDisplayProductOnly ? "True" : "False", false, false);
	ImageList_Destroy(m_hImageList);
	::DeleteObject(m_hFont);
	m_spExplorer = NULL;
	m_spExcel = NULL;
	
	try {		
		m_spTaurus = NULL;
	} catch (...){
		// Probably already been deleted. Nullify the object so you don't get a crash when m_spTaurus goes out of scope.
		ATLTRACE("Exeption encountered when attempting to delete the Taurus pointer!\n");
		std::memset(&(m_spTaurus.p), NULL, sizeof(ITaurus));
	}
	try {
		m_spConnection.Detach();		// ToDo - not including this causes a crash and we shouldn't need this try catch block.
	} catch (...){
		// Do nothing
	}
	try {
		m_spSiriusServer = NULL;		// ToDo - we shouldn't need this try catch block.
	} catch (...){
		// Probably already been deleted. Nullify the object so you don't get a crash when m_spSiriusServer goes out of scope.
		ATLTRACE("Exeption encountered when attempting to delete the Sirius Server pointer!\n");
		std::memset(&(m_spSiriusServer.p), NULL, sizeof(ISiriusDataManager));
	}
	
	try {
		m_spBwlEngine = NULL;
	} catch (...){
		// Probably already been deleted. Nullify the object so you don't get a crash when m_spBwlEngine goes out of scope.
		ATLTRACE("Exeption encountered when attempting to delete the Bwl Engine pointer!\n");
		std::memset(&(m_spBwlEngine.p), NULL, sizeof(IBwlEngine));
	}

	try {
		m_spApplication = NULL;
	} catch (...){
		std::memset(&(m_spApplication.p), NULL, sizeof(ISiriusApplication));
	}
	GDA::Shutdown();	
	CComModule::Term();
}

std::string CSiriusComModule::ValidateLocation(const std::string& sz, const connection_parameters* pcp, FileSystemEnum* pfs) const
//	pcp, pfs - nullable
{
	estring								szRet(sz);
	FileSystemEnum						fs;
	
	if (pfs){
		fs = *pfs;
	} else {
		fs = GetFileSystem();
	}	
	szRet.trim();		
	szRet.ReplaceStrInStr(" ", "");	

	if (fs == fsSQLServer && szRet.size()){
		std::map<estring, estring>					asz;
		std::map<estring, estring>::const_iterator	it;
		estring										sz_lc(szRet);
		_ConnectionPtr								spConnection;
				
		sz_lc.lc();
		spConnection = g_pApplication->GetLocations(fs, pcp, NULL, &asz, true, false);
		if ((it = asz.find(sz_lc)) == asz.end()){
			// We need to add szRet to the database.
			if (::MessageBox(::GetActiveWindow(), ("Do you want to add the location '" + szRet + "' to the database?").c_str(), "Add Location", MB_YESNOCANCEL | MB_ICONQUESTION) != IDYES){
				throw "The addition of location '" + szRet + "' was cancelled";
			}
			std::stringstream ssQuery;
			ssQuery << "sp_user_add_location '" << szRet << "'";
			spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
		} else {
			// This location is already in the locations list
			szRet.assign(it->second);
		}
	}

	return szRet;
}
