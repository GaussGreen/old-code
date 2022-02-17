//	handle_assets.cpp : Implementation of CAssets
//
//	Author :			David Cuin
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionfunctions.h"
#include "handle_assets.h"
#include "handle_volatilitystructure.h"
#include "handle_volatilitystructures.h"
#include "handle_zerocurve.h"
#include "handle_correlationmatrix.h"
#include "handle_dividendschedule.h"
#include "handle_spotschedule.h"
#include "siriusapplication.h"


/*static*/ std::set<IAsset*>			CAssets::s_setDontAdd;

/*static*/ void	CAssets::DontAddToMaintainedCollection(CComPtr<IAsset> spAsset)
{	
	s_setDontAdd.insert(spAsset.p);
}

STDMETHODIMP CAssets::GetBaseUnderlyings(IAssets** pVal)
{
	begin_function
	CComPtr<IAssets>					spAssets;
		
	GetBaseUnderlyings(spAssets);
	return spAssets.CopyTo(pVal);
	end_function
}

void CAssets::GetBaseUnderlyings(CComPtr<IAssets>& spAssets)
{
	if (!spAssets) spAssets.CoCreateInstance(CLSID_Assets);
	
	for (ContainerType::const_iterator it = m_coll.begin(); it != m_coll.end(); it++){
		CComPtr<IAsset>		spAsset;	// item in spAssets
		if (spAsset = dynamic_cast<IAsset*>(it->second.pdispVal)){
			CComPtr<IAssets>	spConstituents;	// basket constituents of spAsset
			CAssets*			pConstituents;
			long				nCount = 0L;	// number of assets in the basket
			
			spAsset->get_Assets(&spConstituents);
			if (pConstituents = dynamic_cast<CAssets*>(spConstituents.p)){
				spConstituents->get_Count(&nCount);
			}
			if (nCount){
				// recall this function on this collection
				pConstituents->GetBaseUnderlyings(spAssets);
			} else {
				// spAsset is not a basket - therefore add it to spAssets
				dynamic_cast<CAssets*>(spAssets.p)->MergeSingular(spAsset);
			}
		}
	}	
}

//	Check that we can form the asset with key szKey from an existing date
//  in the collection. We do this by mapping the input key to a form without
//  composite and pay currencies defined and compare that to the list of
//  assets in memory with keys generated in a similar way.
//
HRESULT CAssets::GetClone(const std::string& szKey, CComPtr<IDispatch>& spObject)
{	
	estring											szIdentifierKey = GetIdentifierKey(szKey);
	IdentifierKeyToKeyType::const_iterator			it = m_IdentifierKeyToKey.find(szIdentifierKey);	 // Might as well use this element!

	if (it != m_IdentifierKeyToKey.end()){
		// We can clone it->second
		KeyToHandleType::const_iterator itKeyToHandle = m_KeyToHandle.find(estring::GetBSTR(it->second));
		if (itKeyToHandle == m_KeyToHandle.end()){
			// should never happen
			ATLASSERT(false);
			return E_FAIL;
		}

		// We can clone this one
		CComPtr<IAsset> spParentAsset = dynamic_cast<IAsset*>(m_coll.find(itKeyToHandle->second)->second.pdispVal);
		estring szCompositeCcy = estring::mid(&szKey, "(", ")"); szCompositeCcy.trim();
		estring szPayCcy = estring::mid(&szKey, ">", "@"); szPayCcy.trim();
		CurrencyEnum ceComposite = CEnumMap::GetEnum("CurrencyEnum", LIBID_Sirius, szCompositeCcy, NoCurrency);
		CurrencyEnum cePay = CEnumMap::GetEnum("CurrencyEnum", LIBID_Sirius, szPayCcy, NoCurrency);
		CComPtr<IAsset> spAsset;
		CAsset::CloneAsset(spParentAsset, ceComposite, cePay, spAsset);
		spObject = spAsset.p;
		return spObject ? S_OK : E_FAIL;		
	}
	return E_FAIL;
}

//	Transforms A (B) > C@Date.DS to A@Date.DS. This is used to
//  store the clone map.
std::string CAssets::GetIdentifierKey(const std::string& szKey)
{	
	estring szA = estring::left(&szKey, "@");
	estring szB = estring::right(&szKey, "@");
						
	if (szA.find("(") != szA.npos){
		szA.assign(szA.left("("));
	} else if (szA.find(">") != szA.npos){
		szA.assign(szA.left(">"));
	}
	szA.trim();
	szB.trim();

	if (!szB.size()){
		return szA;
	} else {
		return szA + "@" + szB;
	}
}

//	Returns the set of volatility structures on which the assets collection depends.
STDMETHODIMP CAssets::GetVolatilityStructures(IVolatilityStructures** pVal)
{	
	begin_function
		
	HRESULT								hr;
	CComPtr<IAssets>					spAssets;
	CComPtr<IVolatilityStructures>		spVolatilityStructures;
	CVolatilityStructures*				pVolatilityStructures;
	long								nAssets = 0L;
	
	GetBaseUnderlyings(spAssets);
	if (hr = spVolatilityStructures.CoCreateInstance(CLSID_VolatilityStructures)) return hr;
	pVolatilityStructures = dynamic_cast<CVolatilityStructures*>(spVolatilityStructures.p);
	ATLASSERT(pVolatilityStructures);	
	
	spAssets->get_Count(&nAssets);
	for (CComVariant vAsset = 1L; vAsset.lVal <= nAssets; vAsset.lVal++){
		CComPtr<IAsset> spAsset;
		CComPtr<IVolatilityStructure> spVolatilityStructure;
		if (spAssets->get_Item(vAsset, &spAsset)) continue;
		if (spAsset->get_VolatilityStructure(&spVolatilityStructure)) continue;
		pVolatilityStructures->MergeSingular(spVolatilityStructure);
	}

	return spVolatilityStructures.CopyTo(pVal);
	end_function
}

/*static*/ void CAssets::InsertCorrelations(CComPtr<IAssets> spUnderlyings)
{
	long							nUnderlyings = 0;
	std::vector<estring>			aszNames;
	DATE							date = 0.0;
	DataSourceEnum					ds = NoDataSource;
			
	if (spUnderlyings->get_Count(&nUnderlyings)) return;
	for (long n = 1; n <= nUnderlyings; n++){
		CComPtr<IAsset>		sp;
		CComBSTR			sName;		
		spUnderlyings->get_Item(CComVariant(n), &sp);
		if (sp){
			sp->get_Identifier(&sName);
			aszNames.push_back(estring(sName));

			if (!(long)date){
				sp->get_Date(&date);
			}
			if (ds == NoDataSource){
				sp->get_DataSource(&ds);
			}
		}
	}
	if (aszNames.size() == nUnderlyings){
		// don't bother doing this if the acquisition of underlyings failed
		CAssets::InsertCorrelations(aszNames, ds, date);
	}

}

/*static*/ void CAssets::InsertCorrelations(std::vector<MlEqAssetHandle> ahAssets)
{
	std::vector<estring>				aszNames;
	
	if (!ahAssets.size()) return;	
	for (long n = 0; n < ahAssets.size(); n++){
		aszNames.push_back(ahAssets[n]->GetName());
	}	
	InsertCorrelations(aszNames, ahAssets[0]->GetDataSource(), ahAssets[0]->GetDateHandle()->GetDate());	
}

/*static*/ void CAssets::InsertCorrelations(const std::vector<estring>& aszNames, DataSourceEnum ds, long nDate)
{
	for (long nRow = 0; nRow < aszNames.size(); nRow++){
		for (long nCol = nRow + 1; nCol < aszNames.size(); nCol++){
			CSiriusApplication::InsertCorrelation(aszNames[nRow], aszNames[nCol], ds, nDate, false);
			CSiriusApplication::InsertCorrelation("~" + aszNames[nRow], "~" + aszNames[nCol], ds, nDate, false);
		}
	}
}

bool CAssets::IsClonable(void) const
{
	return true;
}

void CAssets::KeyAdded(const std::string& szKey)
{
	estring szIdentifierKey = GetIdentifierKey(szKey);	
	m_IdentifierKeyToKey.insert(std::pair<std::string, std::string>(szIdentifierKey, szKey));
}

void CAssets::KeyCleared(const std::string& szKey)
{
	estring											szIdentifierKey = GetIdentifierKey(szKey);			
	IdentifierKeyToKeyType::iterator				itLower = m_IdentifierKeyToKey.find(szIdentifierKey);
	
	if (itLower != m_IdentifierKeyToKey.end()){				
		for (IdentifierKeyToKeyType::iterator it = itLower; it != m_IdentifierKeyToKey.end(); ++it){
			if (it->second == szKey){
				m_IdentifierKeyToKey.erase(it);				
				return;
			}
		}
	}

	// Should never reach this point since {szIdentifierKey, szKey} should always have existed in m_IdentifierKeyToKey
	ATLASSERT(false);
}

void CAssets::KeyRekeyed(const std::string& szKeyOld, const std::string& szKeyNew)
{	
	KeyCleared(szKeyOld);
	KeyAdded(szKeyNew);
}

void CAssets::KeysAllCleared(void)
{	
	m_IdentifierKeyToKey.clear();
}

/*static*/ HRESULT CAssets::Load(const std::string& szDummy, DataSourceEnum ds, long nDate, CComPtr<IAssets>& spAssets)
{
	HRESULT								hr;	
	FileSystemEnum						fs;
	
	if (szDummy.size()) throw "The identifier value '" + szDummy + "' is invalid for the load assets request";
	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();	
	if (hr = spAssets.CoCreateInstance(CLSID_Assets)) return hr;
	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		std::stringstream				ssQuery;
		std::string						szError = "No assets found on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
		ssQuery << "sp_user_get_assets " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "'";
		return CComObjectCollectionFunctions<IAsset>(&szError).ImplementLoadCollection_SQL(ssQuery, ds, spAssets);
	} else if (fs == fsNTFS){										
		return CComObjectCollectionFunctions<IAsset>().ImplementLoadCollection_NTFS(spAssets, 0, "", _Module.GetLocation(), "", ds, nDate, CAsset::Load);
	}
	return E_FAIL;
}

//  Check the spelling of pIndex to ensure the correct spacing for assets of
//  the form A (B) > C, A (B) or A > C.
//
//	We also spell basket assets correctly (recall that they always have an
//  explicit pay currency).
//
HRESULT CAssets::MapIndex(CComVariant* pIndex)
{	
	if (pIndex->vt == VT_BSTR){
		estring szIdentifier(*pIndex);
		// We don't correct spelling of valid spreadsheet handles.
		bool bSpreadsheetHandle;
		CParameterMap::RemoveCalculationNumberFromHandle(CLSID_Asset, &szIdentifier, &bSpreadsheetHandle);
		if (bSpreadsheetHandle) return S_OK;	// We don't repell this!									
		Respell(&szIdentifier);
	}
	return S_OK;
}

/*static*/ void CAssets::Respell(std::string* pszIndex)
{	
	estring								szIdentifier;	
	estring								szSuffix = estring::right(pszIndex, "@"); szSuffix.trim();

	if (szSuffix.size()){
		szIdentifier = estring::left(pszIndex, "@");
	} else {
		szIdentifier = *pszIndex;
	}
	szIdentifier.trim();
		
	estring								szNaturalAsset = estring::left(&szIdentifier, "(").left(">"); szNaturalAsset.trim();
	estring								szCompositeCcy = estring::mid(&szIdentifier, "(", ")");	szCompositeCcy.trim();
	estring								szPayCcy = estring::mid(&szIdentifier, ">"); szPayCcy.trim();					
	
	pszIndex->assign(szNaturalAsset);
	if (szCompositeCcy.size()) pszIndex->append(" (" + szCompositeCcy + ")");
	if (szPayCcy.size()) pszIndex->append(" > " + szPayCcy);
	if (szSuffix.size()) pszIndex->append("@" + szSuffix);	
}

STDMETHODIMP CAssets::Refresh(VARIANT_BOOL Refresh)
{
	begin_function
	CComObjectCollectionFunctions<IAsset>(&m_coll).Refresh(Refresh);
	end_function
}

/*static*/ void CAssets::RemoveAssetFromDontAddSet(IAsset* pAsset)
{	
	s_setDontAdd.erase(pAsset);
}

bool CAssets::ShouldAddToMaintainedCollection(CComPtr<IDispatch> spObject)
{	
	CComQIPtr<IAsset> spAsset = spObject;	
	ATLASSERT(spAsset);

	if (s_setDontAdd.find(spAsset) == s_setDontAdd.end()){
		return true;
	} else {
		return false;
	}	
}