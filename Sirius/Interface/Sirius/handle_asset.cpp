//	handle_asset.cpp : Implementation of CAsset
//
//	author:			   David Cuin
//
//////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionfunctions.h"
#include "handle_asset.h"
#include "handle_assets.h"
#include "siriusapplication.h"
#include "handle_spotschedule.h"
#include "handle_dividendschedule.h"
#include "MlEqCorrelationMatrix.h"
#include "handle_dateschedule.h"
#include "xmlstreamer.h"
#include "handle_zerocurve.h"
#include "smart.h"
#include "handle_volatilitystructure.h"
#include "handle_correlationmatrix.h"
#include "comobjectcollectionserialisabledefaulter.h"
#include "ComObjectCollectionSerialisableKey.h"

// The rekeying in this object is interesting. Note the following:
// (i)  -  The name is not puttable. So there are no rekeying considerations here.
// (ii) -  The date is puttable. However, the quant code reserves the right to 
//		   take a date handle from the asset (as well as a date). This uses the
//         volatility structure. Therefore, for full rigour, put_VolatilityStructure
//		   needs to call the rekey. put_Assets also needs to call rekey since in the basket
//         case, the date is taken from the 1st basket constituent's volatility
//         structure.
// (iii) - The data source does require rekeying. But we need to adjust the
//         correlation matrix. Therefore, I use the ...rekey_ro macro and 
//         manually implement the put_DataSource function.
// (iv)  - We need to rekey after put_Identifier as that adjusts the asset name.
// (v)   - We need to rekey after put_CurrencyCurve as that adjusts the asset name.
// (vi)  - We need to rekey if any of the currencies change.
// (vii) - We have enhanced the standard Rekey function so it sends the
//         appropriate correlation matrix to MlEqAsset. This needs to be done
//		   every time the date or data source could change.


implement_member_variable_rekey_ro(Asset, Assets, DATE, Date, long, InternalDate);
implement_member_variable_rekey_ro(Asset, Assets, DataSourceEnum, DataSource, DataSourceEnum, DataSource);
implement_serialisable(Asset);


///////////////////////////////////////////////////////////////////////////////
//	CloneAsset
//
//	Create a new asset that is equal to an input parent apart from the given
//  natural, composite and pay currencies.
//
/*static*/ void CAsset::CloneAsset(CComPtr<IAsset> spParentAsset, CurrencyEnum ceCloneComposite, CurrencyEnum ceClonePay, CComPtr<IAsset>& spNewAsset)
{		
	CurrencyEnum						ceParentNatural;
	CurrencyEnum						ceParentComposite;
	CurrencyEnum						ceParentPay;	
	CComPtr<IAsset>						spParentCurrencyAsset;
	CComPtr<IAsset>						spCloneCompositeCurrencyAsset;
	CComPtr<IAsset>						spClonePayCurrencyAsset;	
	CComPtr<IVolatilityStructure>		spVolatilityStructure;
	CComPtr<ISpotSchedule>				spSpotSchedule;
	CComPtr<IDividendSchedule>			spDividendSchedule;
	CComPtr<IAssets>					spAssets;
	CComBSTR							sIdentifier;
	DataSourceEnum						ds = NoDataSource;	
	DATE								date;
	CAssets*							pMaintainedAssets = g_pApplication->GetAssets();
	VARIANT_BOOL						bIsBasket = VARIANT_FALSE;
	
	spParentAsset->IsBasket(&bIsBasket);
	if (spParentAsset->get_DataSource(&ds)) propagate_error;
	if (spParentAsset->get_Date(&date)) propagate_error;
	if (spParentAsset->get_Identifier(&sIdentifier)) propagate_error;
	if (spParentAsset->get_Currency(&ceParentNatural)) propagate_error;			
	if (!bIsBasket && spParentAsset->get_CurrencyAsset(&spParentCurrencyAsset)) propagate_error;								
	if (spParentAsset->get_CompositeCurrency(&ceParentComposite)) propagate_error;
	if (spParentAsset->get_PayCurrency(&ceParentPay)) propagate_error;
	
	spParentAsset->get_VolatilityStructure(&spVolatilityStructure);
	spParentAsset->get_SpotSchedule(&spSpotSchedule);
	spParentAsset->get_DividendSchedule(&spDividendSchedule);
	spParentAsset->get_Assets(&spAssets);
	
	// Check we really need to clone.
	if (ceParentComposite == ceCloneComposite && ceParentPay == ceClonePay){
		spNewAsset = spParentAsset;
		return;
	}

	// If the asset {sIdentifier, date, ds, ceCloneComposite, ceClonePay} already exists in the maintained collection then we can return that one.	
	CComObjectCollectionSerialisableKey key(sIdentifier, ceCloneComposite, ceClonePay, date, ds);
	if (pMaintainedAssets->IsInCollection(key)){
		if (pMaintainedAssets->get_Item(key, &spNewAsset)) propagate_error;
		return;
	}

	// A new asset needs to be created if this point is reached.
	spNewAsset.CoCreateInstance(CLSID_Asset);	
	if (ceParentComposite == ceCloneComposite){
		spParentAsset->get_CompositeCurrencyAsset(&spCloneCompositeCurrencyAsset);
	} else {
		spCloneCompositeCurrencyAsset = CurrencyEnumToCurrencyAsset(ceCloneComposite, date, ds);
	}	
	if (ceParentPay == ceClonePay){
		spParentAsset->get_PayCurrencyAsset(&spClonePayCurrencyAsset);
	} else {
		spClonePayCurrencyAsset = CurrencyEnumToCurrencyAsset(ceClonePay, date, ds);
	}
		
	// Put the new asset properties in the (approximate) order of the IAsset interface (see comments in sirius.idl for reasons).	
	if (spNewAsset->put_Identifier(sIdentifier)) propagate_error;
	if (spNewAsset->put_Assets(spAssets)) propagate_error;
	if (spNewAsset->put_DataSource(ds)) propagate_error;
	if (spNewAsset->put_Date(date)) propagate_error;
	if (spNewAsset->put_DividendSchedule(spDividendSchedule)) propagate_error;
	if (spNewAsset->put_SpotSchedule(spSpotSchedule)) propagate_error;
	if (spNewAsset->put_VolatilityStructure(spVolatilityStructure)) propagate_error;		
	if (spNewAsset->put_CurrencyAsset(spParentCurrencyAsset)) propagate_error;		
	if (spNewAsset->put_CompositeCurrencyAsset(spCloneCompositeCurrencyAsset)) propagate_error;
	if (spNewAsset->put_PayCurrencyAsset(spClonePayCurrencyAsset)) propagate_error;
#	ifdef _DEBUG
		CComObjectCollectionSerialisableKey keyNew(spNewAsset);
		CComObjectCollectionSerialisableKey keyParent(spParentAsset);			
		ATLTRACE("Cloned asset '%s' to create '%s'\n", ((std::string)keyParent).c_str(), ((std::string)keyNew).c_str());	
#	endif

	// We don't ever (for as long as the spNewAsset pointer is valid) want to add this asset to the global assets collection; hence the next line.
	CAssets::DontAddToMaintainedCollection(spNewAsset);
}


///////////////////////////////////////////////////////////////////////////////
//	CreateBasket
//
//	Provides a way of creating an asset basket directly. This is useful
//	for spreadsheets.
//
/*static*/ void CAsset::CreateBasket(VARIANT Identifier, VARIANT DataSourceOpt, VARIANT DateOpt, VARIANT AssetsArr, VARIANT CompositeCurrenciesArr, VARIANT WeightsArr, VARIANT PayCurrency, CComPtr<IAsset>& spAsset)
//	spAsset - (returned) - this is the created asset
{		
	CurrencyEnum						ceComposite;					// Composite currency for the basket.
	
	map_parameter(Identifier, estring, szIdentifier);		
	map_optional_enum_parameter(DataSourceOpt, DataSourceEnum, ds, CComObjectCollectionSerialisableDefaulter::GetDataSource());
	map_optional_parameter(DateOpt, long, nDate, CComObjectCollectionSerialisableDefaulter::GetDate());	
	map_optional_enum_vector_parameter(CompositeCurrenciesArr, CurrencyEnum, aceComposite);	
	map_parameter(WeightsArr, std::vector<double>, afWeights);
	map_enum_parameter(PayCurrency, CurrencyEnum, cePay);
	pin_date(ds, nDate);												// Do this so the next line picks up the correct assets.
	map_com_object_vector_parameter(AssetsArr, Asset, aspAssets);		// Call this one last as it takes the longest.
				
	if (aspAssets.size() != afWeights.size()) throw CStringResource(IDS_VECTOR_SIZES_DIFFER);
	if (aceComposite.size() && aceComposite.size() != 1 && aspAssets.size() != aceComposite.size()) throw CStringResource(IDS_VECTOR_SIZES_DIFFER);
	if (cePay == NoCurrency) throw "Invalid pay currency";
		
	// Get aceComposite and ceComposite
	if (!aceComposite.size()){
		// Take the composite currencies from the input assets.
		aceComposite.resize(aspAssets.size());
		for (long nAsset = 0; nAsset < aspAssets.size(); nAsset++){
			CurrencyEnum ce;
			aspAssets[nAsset].m_T->get_CompositeCurrency(&ce);
			aceComposite[nAsset] = ce;
		}				
	} else if (aceComposite.size() == 1){
		aceComposite.resize(aspAssets.size(), aceComposite[0]);
	}	
	ceComposite = aceComposite[0];
	for (long n = 1; n < aceComposite.size(); n++){
		if (aceComposite[n] != ceComposite){
			// Quanto basket - indicate this by no composite currency.
			ceComposite = NoCurrency;
			break;
		}
	}
		
	// Construct two parameter maps which will form the 'Value' property of a new asset.
	// We need Identifier, Date, DateSource, CompositeCurrencyAsset, PayCurrencyAsset, {Assets}.
	CParameterMap						pmNames(aspAssets.size() + 5, 1);
	CParameterMap						pmValues(aspAssets.size() + 5, 1);

	pmNames.SetValue(0, 0, CStringResource(IDS_HEADER_IDENTIFIER));
	pmNames.SetValue(1, 0, CStringResource(IDS_HEADER_DATE));
	pmNames.SetValue(2, 0, CStringResource(IDS_HEADER_DATA_SOURCE));
	pmNames.SetValue(3, 0, CStringResource(IDS_HEADER_PAYCURRENCYASSET));
	pmNames.SetValue(4, 0, CStringResource(IDS_HEADER_COMPOSITECURRENCYASSET));
	pmValues.SetValue(0, 0, szIdentifier);
	pmValues.SetValue(1, 0, nDate);
	pmValues.SetValue(2, 0, CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds));
	if (cePay != NoCurrency){
		pmValues.SetValue(3, 0, CurrencyEnumToCurrencyAseetName(cePay));
	}
	if (ceComposite != NoCurrency){
		pmValues.SetValue(4, 0, CurrencyEnumToCurrencyAseetName(ceComposite));
	}
	for (long nAsset = 0; nAsset < aspAssets.size(); nAsset++){
		CComPtr<IAsset> spAssetAdj;
		CloneAsset(CComPtr<IAsset>(aspAssets[nAsset].m_T), aceComposite[nAsset], cePay, spAssetAdj);
		pmNames.SetValue(5 + nAsset, 0, CComObjectCollectionSerialisableKey(spAssetAdj).GetKeyAndObjectName(false, false));
		pmValues.SetValue(5 + nAsset, 0, afWeights[nAsset]);
	}
	
	// Create the basket.
	CComVariant v;
	spAsset.CoCreateInstance(CLSID_Asset);
	if (CParameterMap::VariableArgumentListToArray((VARIANT*)&v, 2, pmNames, pmValues)) propagate_error;
	if (spAsset->put_Value(v)) propagate_error;
}

/*static*/ CComPtr<IAsset> CAsset::CurrencyEnumToCurrencyAsset(CurrencyEnum ce, long nDate, DataSourceEnum ds)
{			
	if (ce == NoCurrency) return NULL;	
		
	CComPtr<IAsset>							spCurrencyAsset;
	CAssets*								pMaintainedAssets = g_pApplication->GetAssets();		
	std::string								szCurrencyAssetName = CurrencyNameToCurrencyAssetName(CEnumMap::GetString("CurrencyEnum", LIBID_Sirius, ce));
	CComObjectCollectionSerialisableKey		key(szCurrencyAssetName, nDate, ds);
	bool									bIsInCollection = pMaintainedAssets->IsInCollection(key);			
	
	if (pMaintainedAssets->get_Item(CComObjectCollectionSerialisableKey(szCurrencyAssetName, nDate, bIsInCollection ? ds : NoDataSource), &spCurrencyAsset)) propagate_error;
	return spCurrencyAsset;
}

/*static*/ std::string CAsset::CurrencyEnumToCurrencyAseetName(CurrencyEnum ce)
{
	std::string szCurrency;
	CEnumMap::GetString("CurrencyEnum", LIBID_Sirius, ce, &szCurrency);
	if (!szCurrency.size()) return "";
	return CurrencyNameToCurrencyAssetName(szCurrency);
}

/*static*/ std::string CAsset::CurrencyNameToCurrencyAssetName(const std::string& szName)
{
	// This is the only place I allow myself to explicitly append ".USD"
	return szName + ".USD";
}

HRESULT CAsset::FinalConstruct(void)
{	
	begin_function	
	SetCorrelations();	// Sends a correlation matrix down to the analytics, however incomplete!
	end_function	
}

HRESULT CAsset::FinalRelease(void)
{
	CAssets::RemoveAssetFromDontAddSet(this);
	return S_OK;
}

//
//	Assets only have names if they are valid asset types.
//
//	Eggs have the form Identifier (CompositeCurrency) > PayCurrency where we are
//  only explicit about the composite and pay currencies if they differ from the
//  natural currency.
//
//  Baskets also have the form Identifier (CompositeCurrency) > PayCurrency but
//  we only state the composite currency for composite baskets.
//
estring CAsset::GetName(void) const
{
	AssetTypeEnum						at = GetAssetType();
		
	switch (at){
	case CurrencyAsset:		
		return GetIdentifierFromCurrencyCurve();
	case QuantoBasket:
		{
			// Return Identifier > PayCurrency
			CurrencyEnum cePay = GetPayCurrency();			
			return m_h->GetName() + " > " + CEnumMap::GetString("CurrencyEnum", LIBID_Sirius, cePay);
		}
	case CompositeBasket:
		{
			// Return Identifier (CompositeCurrency) > PayCurrency
			CurrencyEnum cePay = GetPayCurrency();
			CurrencyEnum ceComposite = GetCompositeCurrency();
			return m_h->GetName() + " (" + CEnumMap::GetString("CurrencyEnum", LIBID_Sirius, ceComposite) + ") > " + CEnumMap::GetString("CurrencyEnum", LIBID_Sirius, cePay);
		}
	case PrototypeNaturalAsset: case NaturalAsset: case SingleUnderlying:
		{
			estring							szName(m_h->GetName());
			CComBSTR						sNatural, sComposite, sPay;

			if (m_spzcEggNatural) m_spzcEggNatural->get_Name(&sNatural);
			if (m_spzcEggComposite) m_spzcEggComposite->get_Name(&sComposite);
			if (m_spzcEggPay) m_spzcEggPay->get_Name(&sPay);
			
			estring							szNatural(sNatural);
			estring							szComposite(sComposite);
			estring							szPay(sPay);
			
			if (szComposite.size() && szComposite != szNatural) szName += (" (" + szComposite + ")");
			if (szPay.size() && szPay != szNatural) szName += (" > " + szPay);
			return szName;
		}
	default:
		throw "CAsset::get_Name is not supported for asset type " + CEnumMap::GetString("AssetTypeEnum", LIBID_Sirius, at);
	}
}

STDMETHODIMP CAsset::get_Name(BSTR* pVal)
{		
	begin_function
	return GetName().GetBSTR(pVal);
	end_function
}

// Rationalise without adjusting the current composite and pay currencies
void CAsset::GenerateRealConstituents(void)
{		
	GenerateRealConstituents(GetCompositeCurrency(), GetPayCurrency());
}

// Form the m_spBasketRealAssets, m_h->m_assets and m_h->m_assetWeights from m_spBasketAssets and a composite and pay currency
void CAsset::GenerateRealConstituents(CurrencyEnum ceCompositeIn, CurrencyEnum cePayIn)
{
	long								nAssets;
	CComPtr<IAssets>					spMaintainedAssets;
		
	if (!m_spBasketRealAssets){
		m_spBasketRealAssets.CoCreateInstance(CLSID_Assets);
	} else {		
		m_spBasketRealAssets->Clear();
	}
	
	SetAnalyticToEgg(false);
	if (!m_spBasketAssets) return;
	if (m_spBasketAssets->get_Count(&nAssets)) propagate_error;
	if (!nAssets) return;				// Nothing to rationalise.
	SetAnalyticToEgg(true);
	if (CSiriusApplication::GetAssets(spMaintainedAssets)) propagate_error;

	DataSourceEnum						ds;
	long								nDate;			
	
	SafeGetDateAndDataSource(&nDate, &ds);
	
	// If cePayIn is NoCurrency then set it to the pay currency of the first constituent.
	if (cePayIn == NoCurrency){
		CComPtr<IAsset> spAsset;
		m_spBasketAssets->get_Item(CComVariant(1L), &spAsset);
		spAsset->get_PayCurrency(&cePayIn);
		if (cePayIn == NoCurrency) throw "Invalid pay currency encountered in CAsset::GenerateRealConstituents";
	}

	std::string							szCompositeIn = CEnumMap::GetString("CurrencyEnum", LIBID_Sirius, ceCompositeIn);
	std::string							szPayIn = CEnumMap::GetString("CurrencyEnum", LIBID_Sirius, cePayIn);
						
	m_h->m_assets.clear();
	m_h->m_assetWeights.resize(nAssets);
	m_h->m_bIsBasket = true;
	for (long nAsset = 1; nAsset <= nAssets; nAsset++){
		CComPtr<IAsset>		spAsset;				// (nAsset)th constituent.
		double				fWeight;		
		CurrencyEnum		ceComposite;			// Composite currency of spAsset
		CurrencyEnum		cePay;					// Pay currency of spAsset
		CComBSTR			sIdentifier;			// Identifier of spAsset		
		CComPtr<IAsset>		spAssetAdj;				// Adjusted form of spAsset to match the basket characteristics.
		CurrencyEnum		ceCompositeAdj;			// Composite currency that we assign to spAssetAdj.
		CurrencyEnum		cePayAdj;				// Pay currency that we assign to spAssetAdj.
		
		if (m_spBasketAssets->get_Item(CComVariant(nAsset), &spAsset)) propagate_error;
		if (spAsset->get_CompositeCurrency(&ceComposite)) propagate_error;
		if (spAsset->get_PayCurrency(&cePay)) propagate_error;
		if (m_spBasketAssets->get_Notional(CComVariant(spAsset), &fWeight)) propagate_error;
		if (spAsset->get_Identifier(&sIdentifier)) propagate_error;
		
		if (ceCompositeIn == NoCurrency){
			// This is a quanto basket. Therefore keep the composite currency the same.
			ATLASSERT(cePayIn != NoCurrency);
			ceCompositeAdj = ceComposite;
			cePayAdj = cePayIn;
		} else {
			// This is a composite basket. Therefore set the composite currency to ceCompositeIn.
			ATLASSERT(ceCompositeIn != NoCurrency && cePayIn != NoCurrency);
			ceCompositeAdj = ceCompositeIn;
			cePayAdj = cePayIn;
		}
						
		CloneAsset(spAsset, ceCompositeAdj, cePayAdj, spAssetAdj);			
						
		// Now update m_spBasketRealAssets and MlEqAsset
		m_spBasketRealAssets->Add(CComVariant(), spAssetAdj);
		m_spBasketRealAssets->put_Notional(CComVariant(spAssetAdj), fWeight);	// Unnecessary but done for completeness' sake.		
		map_com_to_analytic(spAssetAdj, Asset, hAssetAdj);		
		m_h->m_assets.push_back(hAssetAdj);		
		m_h->m_assetWeights[nAsset - 1] = fWeight;
	}

	// Now adjust the MlEqAsset object so the natural, pay and composite currencies match the pay currencies of the
	// basket constituents (which must be all equal).
	MlEqZeroCurveHandle					hZeroCurve;
	MlEqAssetHandle						hAsset;	
	hZeroCurve = m_h->m_assets[0]->GetPayZeroCurve(true);
	hAsset = m_h->m_assets[0]->GetFxPayAsset();
	m_h->PutNaturalZeroCurve(hZeroCurve);	
	m_h->PutCompositeZeroCurve(hZeroCurve);
	m_h->PutPayZeroCurve(hZeroCurve);
	m_h->PutFxNaturalAsset(hAsset);
	m_h->PutFxCompositeAsset(hAsset);
	m_h->PutFxPayAsset(hAsset);
}

STDMETHODIMP CAsset::get_Assets(IAssets** pVal)
{	
	return m_spBasketAssets.CopyTo(pVal);
}

AssetTypeEnum CAsset::GetAssetType(void) const
{
	long								nConstituents = 0L;
	long								nRealConstituents = 0L;		// must equal nConstituents for validity.
	
	if (m_spBasketRealAssets) m_spBasketRealAssets->get_Count(&nRealConstituents);
	if (m_spBasketAssets) m_spBasketAssets->get_Count(&nConstituents);
	if (nConstituents != nRealConstituents) return InvalidAsset;
	if (m_spCurrencyCurve){
		// Candidate currency asset. There should be no basket constituents and the interface-level
		// currency asset pointers should all be NULL. We can't check the analytic level pointers
		// since these are often set in a redundant manner to facilitate quant programming.
		if (!nConstituents && !m_spEggAssetNatural && !m_spEggAssetComposite && !m_spEggAssetPay){			
			return CurrencyAsset;
		}
	} else if (m_h->IsBasket()){
		// Is the basket quanto or composite? If the basket is composite if all the composite
		// currencies of the underlyings are equal; else it is quanto.
		if (nConstituents){
			// Otherwise invalid.		
			CurrencyEnum ceTag = NoCurrency;
			for (long n = 1; n <= nConstituents; n++){
				CComPtr<IAsset> spAsset;
				m_spBasketRealAssets->get_Item(CComVariant(n), &spAsset);
				CurrencyEnum ce = NoCurrency;
				spAsset->get_CompositeCurrency(&ce);
				if (ce == NoCurrency){
					return InvalidAsset;					
				} else if (ceTag == NoCurrency){
					ceTag = ce;					
				} else if (ceTag != ce){
					return QuantoBasket;
				}
			}
			return CompositeBasket;
		}
	} else {
		// Egg case. All the quant level zero curves must be set.
		if (!!m_h->GetNaturalZeroCurve(false) && !!m_h->GetCompositeZeroCurve(false) || !!m_h->GetPayZeroCurve(false)){
			if (m_h->GetPayZeroCurve(false) == m_h->GetNaturalZeroCurve(false) && m_h->GetCompositeZeroCurve(false) == m_h->GetNaturalZeroCurve(false) && m_h->HasVolatilityStructure()){
				return NaturalAsset;
			} else {
				return SingleUnderlying;
			}
		} else {
			return PrototypeNaturalAsset;
		}
	}
	return InvalidAsset;
}

STDMETHODIMP CAsset::get_AssetType(AssetTypeEnum* pVal)
{	
	begin_function	
	*pVal = GetAssetType();
	end_function
}

STDMETHODIMP CAsset::GetBaseUnderlyings(IAssets** pVal)
{
	CComPtr<IAssets>					spAssets;
	HRESULT								hr;
	long								nCount = 0L;

	get_Assets(&spAssets);
	if (spAssets) spAssets->get_Count(&nCount);	
	if (!nCount){
		// The asset is not a basket. Therefore return this asset in the output collection
		spAssets = NULL;
		if (hr = spAssets.CoCreateInstance(CLSID_Assets)) return hr;
		spAssets->Add(CComVariant(1L), this);
		return spAssets.CopyTo(pVal);
	}
	return spAssets->GetBaseUnderlyings(pVal);
}

CurrencyEnum CAsset::GetCompositeCurrency(void) const
{	
	if (!m_h->IsBasket()){
		// Egg case
		return GetCurrency(m_h->GetCompositeZeroCurve(true));
	} else {
		// Basket case - the composite currency is only valid for a composite basket; 
		// for the quanto case, return NoCurrency.
		AssetTypeEnum at = GetAssetType();		
		if (at == CompositeBasket){
			// Return the composite currency of the first constituent.
			CComPtr<IAsset> spAsset;
			if (!m_spBasketRealAssets || m_spBasketRealAssets->get_Item(CComVariant(1L), &spAsset)){
				return NoCurrency;
			} else {
				CurrencyEnum ce;
				if (spAsset->get_CompositeCurrency(&ce)) propagate_error;
				return ce;
			}
		} else if (at == QuantoBasket){
			return NoCurrency;
		} else {
			return NoCurrency;
		}
	}	
}
STDMETHODIMP CAsset::get_CompositeCurrency(/*[out, retval]*/ CurrencyEnum* pVal)
{
	begin_function
	if (!m_h->IsBasket()){
		// Egg case
		*pVal = GetCurrency(m_h->GetCompositeZeroCurve(true));
	} else {
		// Basket case - the composite currency is only valid for a composite basket; 
		// for the quanto case, return NoCurrency.
		AssetTypeEnum at = GetAssetType();		
		if (at == CompositeBasket){
			// Return the composite currency of the first constituent.
			CComPtr<IAsset> spAsset;
			if (!m_spBasketRealAssets || m_spBasketRealAssets->get_Item(CComVariant(1L), &spAsset)){
				*pVal = NoCurrency;
			} else {						
				return spAsset->get_CompositeCurrency(pVal);
			}
		} else if (at == QuantoBasket){
			*pVal = NoCurrency;
		} else {
			*pVal = NoCurrency;
		}
	}
	end_function
}

STDMETHODIMP CAsset::get_CompositeCurrencyAsset(IAsset** pVal)
{
	begin_function
	if (!m_h->IsBasket()){
		// Egg case	
		return m_spEggAssetComposite.CopyTo(pVal);
	} else {
		// Basket case - the composite currency asset is only valid for a composite basket.
		AssetTypeEnum at = GetAssetType();		
		if (at == CompositeBasket){
			// The composite currencies of all the basket (real) assets are all equal.
			// We can therefore take the first one.
			CComPtr<IAsset>	spConstituent, spCurrency;
			m_spBasketRealAssets->get_Item(CComVariant(1L), &spConstituent);
			spConstituent->get_CompositeCurrencyAsset(&spCurrency);
			if (!spCurrency){
				// The constituent is a natural asset. Therefore take the natural currency.
				spConstituent->get_CurrencyAsset(&spCurrency);
			}
			return spCurrency.CopyTo(pVal);
		} else if (at == QuantoBasket){						
			throw "No composite currency asset is defined for '" + GetName() + "' since it is a quanto basket";			
		} else {
			throw "Unhandled exception in CAsset::get_CompositeCurrencyAsset";
		}		
	}
	end_function
}

STDMETHODIMP CAsset::get_ConstituentWeight(IAsset* Constituent, double *pVal)
{
	if (!m_spBasketAssets) return E_POINTER;
	return m_spBasketAssets->get_Notional(CComVariant(Constituent), pVal);
}

STDMETHODIMP CAsset::get_Currency(/*[out, retval]*/ CurrencyEnum* pVal)
{
	begin_function
	if (!m_h->IsBasket()){
		// Egg case
		*pVal = GetCurrency(m_h->GetNaturalZeroCurve(true));
	} else {
		// Basket case - the natural currency is conceptually meaningless so we return NoCurrency
		*pVal = NoCurrency;
	}
	end_function
}

CurrencyEnum CAsset::GetCurrency(MlEqZeroCurveHandle hZeroCurve) const
{
	if (!hZeroCurve) return NoCurrency;
	return CEnumMap::GetEnum("CurrencyEnum", LIBID_Sirius, hZeroCurve->GetName(), NoCurrency);
}

STDMETHODIMP CAsset::get_CurrencyAsset(IAsset** pVal)
{
	begin_function
	if (!m_h->IsBasket()){
		// Egg case	
		return m_spEggAssetNatural.CopyTo(pVal);
	} else {
		// Basket case - not defined
		throw "No currency asset defined for '" + GetName() + "' since it is a basket";
	}
	end_function
}

STDMETHODIMP CAsset::get_CurrencyCurve(IZeroCurve** newVal)
{
	begin_function
	return m_spCurrencyCurve.CopyTo(newVal);
	end_function
}

STDMETHODIMP CAsset::get_DividendSchedule(IDividendSchedule** pVal)
{
	return m_spDividendSchedule.CopyTo(pVal);
}

STDMETHODIMP CAsset::GetForward(DATE Today, DATE Maturity, double *pVal)
{	
	begin_function
	*pVal = m_h->GetForward(Today, Maturity, false);
	end_function
}

STDMETHODIMP CAsset::get_Identifier(BSTR* pVal)
{
	return estring(m_h->GetName()).GetBSTR(pVal);
}

estring CAsset::GetIdentifierFromCurrencyCurve(void) const
{
	CComBSTR sIdentifier;
	if (m_spCurrencyCurve->get_Name(&sIdentifier)) propagate_error;
	return CurrencyNameToCurrencyAssetName(estring(sIdentifier));
}

STDMETHODIMP CAsset::get_IsCalibrated(VARIANT_BOOL *pVal)
{
	begin_function
	*pVal = m_h->IsCalibrated() ? VARIANT_TRUE : VARIANT_FALSE;	
	end_function
}

STDMETHODIMP CAsset::get_NaturalSpot(DATE Date, double* pVal)
{
	try {
		*pVal = m_h->GetNaturalSpot(Date);
		return S_OK;
	} catch (...){
		// assume this is due to no data being defined on date Date
		std::string szDate;
		CParameterMap::DateToString(Date, &szDate);
		return CParameterMap::ReturnErrorRS(IDS_NO_VALUE_AT_DATE, szDate, IID_IAsset);
	}
}

CurrencyEnum CAsset::GetPayCurrency(void) const
{
	if (!m_h->IsBasket()){
		// Egg case
		return GetCurrency(m_h->GetPayZeroCurve(true));
	} else {
		// Basket case - the pay currency is equal to the pay currency of the constituents (which are always identical).
		CComPtr<IAsset> spAsset;
		if (!m_spBasketRealAssets || m_spBasketRealAssets->get_Item(CComVariant(1L), &spAsset)){
			return NoCurrency;
		} else {			
			CurrencyEnum ce;
			if (spAsset->get_PayCurrency(&ce)) propagate_error;
			return ce;
		}
	}
}

STDMETHODIMP CAsset::get_PayCurrency(/*[out, retval]*/ CurrencyEnum* pVal)
{
	begin_function
	*pVal = GetPayCurrency();	
	end_function
}

STDMETHODIMP CAsset::get_PayCurrencyAsset(IAsset** pVal)
{
	begin_function
	if (!m_h->IsBasket()){
		// Egg case	
		return m_spEggAssetPay.CopyTo(pVal);
	} else {
		// Basket case - the pay currency asset is equal to the pay currency asset of the constituents (which are always identical).		
		CComPtr<IAsset>	spConstituent, spCurrency;
		m_spBasketRealAssets->get_Item(CComVariant(1L), &spConstituent);
		spConstituent->get_PayCurrencyAsset(&spCurrency);
		if (!spCurrency){
			// The constituent is a natural asset. Therefore return the natural currency.
			spConstituent->get_CurrencyAsset(&spCurrency);
		}
		return spCurrency.CopyTo(pVal);
	}
	end_function
}

STDMETHODIMP CAsset::get_RealAssets(IAssets** pVal)
{
	return m_spBasketRealAssets.CopyTo(pVal);
}

STDMETHODIMP CAsset::get_Spot(DATE Date, double* pVal)
{
	try {
		*pVal = m_h->GetSpot(Date);
		return S_OK;
	} catch (...){
		// assume this is due to no data being defined on date Date
		std::string szDate;
		CParameterMap::DateToString(Date, &szDate);
		return CParameterMap::ReturnErrorRS(IDS_NO_VALUE_AT_DATE, szDate, IID_IAsset);
	}
}

STDMETHODIMP CAsset::GetSpots(DATE Today, IDateSchedule* Dates, ISpotSchedule** pVal)
{		
	HRESULT								hr;
	CComPtr<ISpotSchedule>				spSpotSchedule;		
	std::vector<long>					vectorDates;	
	CDateSchedule*						pDates = dynamic_cast<CDateSchedule*>(Dates);

	begin_function
	if (!m_spSpotSchedule) return CParameterMap::ReturnErrorRS(IDS_NO_SPOT_SCHEDULE_IN_ASSET, m_h->GetName(), IID_IAsset);
	if (!pDates) return CParameterMap::ReturnErrorR(IDS_SCHEDULE_EMPTY, IID_IAsset);
	pDates->m_h->GetDates(vectorDates);		
	if (hr = spSpotSchedule.CoCreateInstance(CLSID_SpotSchedule)) return hr;
	CSpotSchedule*						pSpotSchedule = dynamic_cast<CSpotSchedule*>(spSpotSchedule.p);		
	m_h->GetNaturalSpots(Today, vectorDates, *pSpotSchedule->m_h);		
	return spSpotSchedule.CopyTo(pVal);
	end_function	
}

STDMETHODIMP CAsset::get_SpotSchedule(ISpotSchedule** pVal)
{
	if (!m_spSpotSchedule){
		*pVal = NULL;
		return CParameterMap::ReturnErrorRS(IDS_NO_SPOT_SCHEDULE_IN_ASSET, m_h->GetName(), IID_IAsset);
	}
	return m_spSpotSchedule.CopyTo(pVal);
}

STDMETHODIMP CAsset::get_Value(VARIANT* pVal)
{							
	return g_pApplication->GetObjectManager().ImplementGetValue(this, pVal);
}

STDMETHODIMP CAsset::get_VolatilityStructure(IVolatilityStructure** pVal)
{
	return m_spVolatilityStructure.CopyTo(pVal);
}

STDMETHODIMP CAsset::GetZeroCurve(ZeroCurveTypeEnum ZeroCurveType, IZeroCurve** pVal)
{		
	begin_function
	if (!m_h->IsBasket()){	
		// Egg case
		switch (ZeroCurveType){
		case Natural:
			return m_spzcEggNatural.CopyTo(pVal);		
		case Pay:
			return m_spzcEggPay.CopyTo(pVal);
		case Composite:
			return m_spzcEggComposite.CopyTo(pVal);
		}
		return CParameterMap::ReturnErrorRN(IDS_INVALID_ZERO_CURVE_TYPE, ZeroCurveType, IID_IAsset);
	} else {
		// Basket case
		switch (ZeroCurveType){
		case Natural:
			throw "No natural zero curve is defined for asset '" + GetName() + "' since it is a basket";
		case Pay:
			// The pay zero curve asset is equal to the pay zero curve of of the constituents (which are always identical).
			{
				CComPtr<IAsset>	spAsset;
				m_spBasketRealAssets->get_Item(CComVariant(1L), &spAsset);		
				return spAsset->GetZeroCurve(ZeroCurveType, pVal);
			}
		case Composite:
			{
				AssetTypeEnum at = GetAssetType();				
				if (at == CompositeBasket){
					CComPtr<IAsset> spAsset;
					m_spBasketRealAssets->get_Item(CComVariant(1L), &spAsset);
					return spAsset->GetZeroCurve(ZeroCurveType, pVal);				
				} else if (at == QuantoBasket){						
					throw "No composite zero curve is defined for '" + GetName() + "' since it is a quanto basket";
				} else {
					throw "Unhandled exception in CAsset::GetZeroCurve";
				}
			}		
		}
	}
	end_function
}

STDMETHODIMP CAsset::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IAsset };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CAsset::IsBasket(VARIANT_BOOL* pVal)
{
	begin_function
	*pVal = m_h->IsBasket() ? VARIANT_TRUE : VARIANT_FALSE;
	end_function
}

// Loads an asset from the database, creating one on the fly if possible
/*static*/ HRESULT CAsset::Load(const std::string& szIdentifier, DataSourceEnum ds, long nDate, CComPtr<IAsset>& spAsset)
{	
	FileSystemEnum						fs;
	estring								szTrim(estring::mid(&szIdentifier, "Asset::"));

	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();
	
	// load the XML for the asset	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){		
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;						// SQL query to execute
		_RecordsetPtr					prs;							// ADO recordset
		CComVariant						vField;							// asset field value
						
		ssQuery << "sp_user_get_asset '" << (szTrim.size() ? szTrim : szIdentifier) << "', " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "'";
		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
			if (prs->adoEOF) throw "Cannot find asset '" + szIdentifier + "' on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
			vField = prs->GetFields()->GetItem(CComVariant(0, VT_I4))->GetValue();						
		} catch (_com_error& e){
			throw estring(e);
		}
		if (vField.vt != VT_BSTR) return E_FAIL;				
		return CXmlStreamer::GetObject((char*)_bstr_t(vField.bstrVal), ds, spAsset);
	} else if (fs == fsNTFS){
		CComObjectCollectionFunctions<IAsset>::ImplementLoad_NTFS(IID_IAsset, _Module.GetLocation(), "", szTrim.size() ? szTrim : szIdentifier, ds, nDate, spAsset);
		return S_OK;
	} else {
		return E_FAIL;
	}	
}

STDMETHODIMP CAsset::put_Assets(IAssets* newVal)
{						
	begin_function
	AssetTypeEnum						at = GetAssetType();	
	std::string							szOldKey = CComObjectCollectionSerialisableKey(this);
	long								nCount;
			
	if (at == CurrencyAsset){
		if (newVal){
			throw "Cannot set the assets collection for asset '" + GetName() + "' since it has a currency curve defined";
		} else {
			return S_OK;
		}
	}
	
	if (!(m_spBasketAssets = newVal) || m_spBasketAssets->get_Count(&nCount) || !nCount){
		// Non-basket base - the analytic object has already been set.		
		// It is the caller's responsibility to deal with setting asset currencies.
		SetAnalyticToEgg(false);
	} else {
		// Basket case
		SetAnalyticToEgg(true);
		GenerateRealConstituents();
	}
	Rekey(szOldKey);	// We need to rekey to (i) set up the correlations and (ii) the date might have changed (see above).
	end_function
}

STDMETHODIMP CAsset::put_CompositeCurrencyAsset(IAsset* newVal)
{
	begin_function
	AssetTypeEnum						at = GetAssetType();	
		
	if (at == CurrencyAsset){
		if (newVal){
			throw "Cannot set the composite currency for asset '" + GetName() + "' since it has a currency curve defined";
		} else {
			return S_OK;
		}
	}
	
	if (m_h->IsBasket()){
		// Basket case		
		CurrencyEnum ceComposite = NoCurrency, cePay = GetPayCurrency();		
		if (newVal) newVal->get_CompositeCurrency(&ceComposite);
		GenerateRealConstituents(ceComposite, cePay);
	} else {		
		MlEqZeroCurveHandle					hZeroCurve;
		MlEqAssetHandle						hAsset;
		std::string							szOldKey = CComObjectCollectionSerialisableKey(this);

		PutEggCurrency(newVal, m_spzcEggComposite, hZeroCurve, m_spEggAssetComposite, hAsset);
		m_h->PutCompositeZeroCurve(hZeroCurve);
		m_h->PutFxCompositeAsset(hAsset);
		SynchroniseEggCurrencies(szOldKey);
	}
	end_function
}

STDMETHODIMP CAsset::put_ConstituentWeight(IAsset* Constituent, double newVal)
{
	begin_function
	HRESULT								hr;

	if (!m_spBasketAssets) return E_POINTER;
	if (hr = m_spBasketAssets->put_Notional(CComVariant(Constituent), newVal)) return hr;
	GenerateRealConstituents();
	end_function
}

STDMETHODIMP CAsset::put_CurrencyAsset(IAsset* newVal)
{
	begin_function
	AssetTypeEnum						at = GetAssetType();	
		
	if (at == CurrencyAsset){
		if (newVal){
			throw "Cannot set the currency for asset '" + GetName() + "' since it has a currency curve defined";
		} else {
			return S_OK;
		}
	}		
	
	if (m_h->IsBasket()){
		if (newVal) throw "The currency asset cannot be set for a basket asset";
		// We don't need to do any processing for the newVal == NoCurrency case
		return S_OK;
	} else {		
		MlEqZeroCurveHandle					hZeroCurve;
		MlEqAssetHandle						hAsset;
		std::string							szOldKey = CComObjectCollectionSerialisableKey(this);

		PutEggCurrency(newVal, m_spzcEggNatural, hZeroCurve, m_spEggAssetNatural, hAsset);
		m_h->PutNaturalZeroCurve(hZeroCurve);
		m_h->PutFxNaturalAsset(hAsset);
		SynchroniseEggCurrencies(szOldKey);
	}
	end_function
}

STDMETHODIMP CAsset::put_CurrencyCurve(IZeroCurve* newVal)
{	
	begin_function
	std::string							szOldKey = CComObjectCollectionSerialisableKey(this);
	MlEqZeroCurveHandle					hCurrencyCurve;
			
	if (m_spCurrencyCurve = newVal){
		// The identifier needs to be adjusted for non-null input.
		if (put_Identifier((CComBSTR)GetIdentifierFromCurrencyCurve())) throw "Unhandled exception in CAsset::put_CurrencyCurve";
	}
	
	// Put the MlEqAsset members
	m_h->PutFxNaturalAsset(m_h.dumb());
	m_h->PutFxCompositeAsset(m_h.dumb());
	m_h->PutFxPayAsset(m_h.dumb());
	CParameterMap::ComToAnalytic(m_spCurrencyCurve, hCurrencyCurve);
	m_h->PutNaturalZeroCurve(hCurrencyCurve);
	m_h->PutCompositeZeroCurve(hCurrencyCurve);
	m_h->PutPayZeroCurve(hCurrencyCurve);

	// The object needs to be rekeyed (see notes at top of file).
	Rekey(szOldKey);
	end_function
}

STDMETHODIMP CAsset::put_Date(DATE newVal)
{	
	begin_function	
	if ((long)newVal != m_h->GetInternalDate()){
		std::string szOldKey = CComObjectCollectionSerialisableKey(this);
		m_h->PutInternalDate((long)newVal);
		Rekey(szOldKey);
	} else if (newVal){				
		// Include this case. This assists CObjectManager::Refresh() which reputs the date and data source to perform any extra cleanup.
		GenerateRealConstituents();
		SetCorrelations();	
	}
	end_function																											\
}

STDMETHODIMP CAsset::put_DataSource(DataSourceEnum newVal)
{
	begin_function
	if (newVal != m_h->GetDataSource()){
		std::string szOldKey = CComObjectCollectionSerialisableKey(this);
		m_h->PutDataSource(newVal);
		Rekey(szOldKey);
	} else if (newVal != NoDataSource){				
		// Include this case. This assists CObjectManager::Refresh() which reputs the date and data source to perform any extra cleanup.
		GenerateRealConstituents();
		SetCorrelations();
	}
	end_function																											\
}

STDMETHODIMP CAsset::put_DividendSchedule(IDividendSchedule* pVal)
{
	begin_function
	m_spDividendSchedule = pVal;
	CDividendSchedule* p = dynamic_cast<CDividendSchedule*>(m_spDividendSchedule.p);
	m_h->PutDividendSchedule(p ? p->m_h : NULL);
	end_function
}

void CAsset::PutEggCurrency(IAsset* newVal, CComPtr<IZeroCurve>& spZeroCurve, MlEqZeroCurveHandle& hZeroCurve, CComPtr<IAsset>& spAsset, MlEqAssetHandle& hAsset) const
{				
	if (!newVal){
		// Clear the zero curves.
		spZeroCurve = NULL;
		hZeroCurve = NULL;
		// Clear the currency assets.
		spAsset = NULL;
		hAsset = NULL;
	} else {				
		// You can only put a currency asset.
		AssetTypeEnum		at;
		
		newVal->get_AssetType(&at);
		if (at != CurrencyAsset) throw "The asset '" + GetName() + "' is not a currency asset";

		CComPtr<IAsset>		spAsset_New;
		MlEqAssetHandle		hAsset_New;
		CComPtr<IZeroCurve> spZeroCurve_New;
		MlEqZeroCurveHandle	hZeroCurve_New;				
		long				nDate;
		DataSourceEnum		ds;		

		SafeGetDateAndDataSource(&nDate, &ds);
		
		// Set the currency asset (spAsset, hAsset)
		spAsset_New = newVal;
		CParameterMap::ComToAnalytic(spAsset_New, hAsset_New);		
				
		// Set the currency curve (spZeroCurve, hZeroCurve)
		if (spAsset_New->get_CurrencyCurve(&spZeroCurve_New)) propagate_error;
		CParameterMap::ComToAnalytic(spZeroCurve_New, hZeroCurve_New);

		// All is well if this point is reached.
		spZeroCurve = spZeroCurve_New;
		hZeroCurve = hZeroCurve_New;
		spAsset = spAsset_New;
		hAsset = hAsset_New;		
	}
}

STDMETHODIMP CAsset::put_Identifier(BSTR newVal)
{	
	begin_function	
	AssetTypeEnum						at = GetAssetType();
	estring								szNewValue(newVal);
	
	szNewValue.trim();
	if (szNewValue == m_h->GetName()) return S_OK;

	std::string szOldKey = CComObjectCollectionSerialisableKey(this);	
	if (at == CurrencyAsset){
		// The identifier must be identical to the currency curve.
		estring szIdentifier = GetIdentifierFromCurrencyCurve();
		if (szIdentifier.CompareNoCase(szNewValue)) throw "The identifier value '" + szNewValue + "' is inconsistent with the currency asset curve";
		szNewValue.assign(szIdentifier);
	}

	if (!szNewValue.size()) return CParameterMap::ReturnErrorR(IDS_ASSET_NAME_INVALID, IID_IAsset);
	m_h->PutName(szNewValue);
	Rekey(szOldKey);	// We need to call this here; see the comment at the top of this file.
	end_function
}

STDMETHODIMP CAsset::put_NaturalSpot(DATE Date, double newVal)
{
	return m_spSpotSchedule->put_ValueAt(Date, newVal);	
}

STDMETHODIMP CAsset::put_PayCurrencyAsset(IAsset* newVal)
{
	begin_function
	AssetTypeEnum						at = GetAssetType();	
		
	if (at == CurrencyAsset){
		if (newVal){
			throw "Cannot set the pay currency for asset '" + GetName() + "' since it has a currency curve defined";
		} else {
			return S_OK;
		}
	}		
	
	if (m_h->IsBasket()){
		// Basket case.		
		CurrencyEnum ceComposite = GetCompositeCurrency(), cePay = NoCurrency;
		if (newVal) newVal->get_PayCurrency(&cePay);
		if (cePay == NoCurrency) throw "You cannot set the pay currency of '" + GetName() + "' to 'NoCurrency' since it is a basket";		
		GenerateRealConstituents(ceComposite, cePay);
	} else {		
		MlEqZeroCurveHandle					hZeroCurve;
		MlEqAssetHandle						hAsset;
		std::string							szOldKey = CComObjectCollectionSerialisableKey(this);

		PutEggCurrency(newVal, m_spzcEggPay, hZeroCurve, m_spEggAssetPay, hAsset);
		m_h->PutPayZeroCurve(hZeroCurve);
		m_h->PutFxPayAsset(hAsset);
		SynchroniseEggCurrencies(szOldKey);
	}
	end_function
}

STDMETHODIMP::CAsset::put_SpotSchedule(ISpotSchedule* pVal)
{
	if (m_spSpotSchedule = pVal){
		m_h->PutSpotSchedule(dynamic_cast<CSpotSchedule*>(m_spSpotSchedule.p)->m_h);
	} else {
		m_h->PutSpotSchedule(NULL);
	}
	return S_OK;
}

STDMETHODIMP CAsset::put_Value(VARIANT newVal)
{	
	begin_function		
	return g_pApplication->GetObjectManager().ImplementPutValue(newVal, this);	
	end_function
}

STDMETHODIMP CAsset::put_VolatilityStructure(IVolatilityStructure* pVal)
{
	std::string szOldKey = CComObjectCollectionSerialisableKey(this);
	if (m_spVolatilityStructure = pVal){		
		m_h->PutVolatilityStructure((dynamic_cast<CVolatilityStructure*>(m_spVolatilityStructure.p))->m_h);		
	} else {
		m_h->PutVolatilityStructure(NULL);
	}
	Rekey(szOldKey);	// We need to call this here; see the comment at the top of this file.
	return S_OK;
}

STDMETHODIMP CAsset::Refresh(VARIANT_BOOL Recursive)
{
	begin_function
	g_pApplication->GetObjectManager().Refresh(this, Recursive ? true : false);			
	end_function
}

void CAsset::Rekey(const std::string& szOldKey)
{
	CComPtr<IAssets>					spAssets;
	CComPtr<IAsset>						spAsset = dynamic_cast<IAsset*>(this);
	
	if (g_pApplication->GetObjectManager().IIDToCollectionPtr(IID_IAsset, spAssets)) throw "Unhandled exception in CAsset::Rekey";	
	if (!m_h->GetInternalDate()){
		// Explicitly set the internal asset date if it's zero and we can!				
		long nDate = 0L;		
		if (m_h->HasVolatilityStructure()){
			MlEqConstDateHandle hDate;
			if (!!(hDate = m_h->GetDateHandle())){
				nDate = hDate->GetDate();
			}
		}		
		m_h->PutInternalDate(nDate);
	}
	g_pApplication->GetObjectManager().Rekey(spAsset, spAssets, szOldKey, NULL);	
	SetCorrelations();	
}

void CAsset::SafeGetDateAndDataSource(long* pnDate, DataSourceEnum* pds) const
{	
	*pds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	*pnDate = CComObjectCollectionSerialisableDefaulter::GetDate();
	if (!!m_h){
		DataSourceEnum	ds;
		long			nDate;
		if ((ds = m_h->GetDataSource()) != NoDataSource) *pds = ds;
		if (nDate = m_h->GetInternalDate()) *pnDate = nDate;
	}	
}

STDMETHODIMP CAsset::Save(BSTR* pVal)
{
	HRESULT								hr;	
	xmlstreamer							ssXML;							// XML representation of vData	
	FileSystemEnum						fs;
	estring								szName(GetName());
	AssetTypeEnum						at = GetAssetType();

	begin_function
	check_publishing_enabled
	if (!m_h->GetInternalDate()) throw "You cannot save the asset '" + m_h->GetName() + "' since it does not have a date.";
	
	hr = CXmlStreamer::GetXML(CComPtr<IAsset>(this), ssXML);
	if (hr) return hr;

	// check to see if the asset is valid	
	if (at == InvalidAsset){
		std::string						szEnum;
		CEnumMap::GetString("AssetTypeEnum", LIBID_Sirius, at, &szEnum);
		return CParameterMap::ReturnErrorRS(IDS_ASSET_SAVE_INVALID, szEnum);
	}

	// save it using the appropriate stored procedure
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;						// SQL query to execute	
		_RecordsetPtr					prs;							// ADO recordset
				
		if (m_h->GetDataSource() != Last) throw "Assets with data sources other than 'Last' cannot be saved. The asset '" + GetName() + "' has the source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "'";
		ssQuery << "sp_user_add_asset '" << GetName() << "', '" << m_h->GetInternalDate() << "', '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) << "', '" << _Module.GetLocation() << "', '" << (char*)ssXML << "'";
		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
		} catch (_com_error& e){
			throw estring(e);
		}
		// set the return type as the asset name (which is the returned value of the stored procedure)
		CComVariant vValue(prs->GetFields()->GetItem(CComVariant(0, VT_I4))->GetValue());
		if (vValue.vt != VT_BSTR){
			ATLASSERT(false);
			return CParameterMap::ReturnErrorR(IDS_UNHANDLED_EXCEPTION);
		}		
		return CComBSTR(vValue.bstrVal).CopyTo(pVal);
	} else if (fs == fsNTFS){		
		// we relax the condition that only 'Last' data sources can be saved (or the UpdatePL process gets tricky)
		if (m_h->GetDataSource() == NoDataSource) throw "The asset '" + GetName() + "' has the data source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "' and so cannot be saved";
		return CComObjectCollectionFunctions<IAsset>::ImplementSave_NTFS(ssXML, IID_IAsset, _Module.GetLocation(), szName, m_h->GetDataSource(), m_h->GetInternalDate(), pVal);
	} else {
		return E_FAIL;
	}	
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	SetAnalyticToEgg
//
//	Sets the analytic level asset parameters to the Non-basket case.
//
void CAsset::SetAnalyticToEgg(bool bDeclareAsBasket)
{
	m_h->m_bIsBasket = bDeclareAsBasket;
	m_h->m_assets.clear();
	m_h->m_assetWeights.resize(1);				
	m_h->m_assets.push_back(m_h.dumb());	
	m_h->m_assetWeights[0] = 1.0;	
}


// Ensures that all the appropriate correlations are defined in the asset.
void CAsset::SetCorrelations(void)
{
	DataSourceEnum						ds;
	long								nDate;	

	// Send the appropriate correlation matrix to the analytics	
	SafeGetDateAndDataSource(&nDate, &ds);	
	m_h->PutCorrelationMatrix(CSiriusApplication::GetCorrelationMatrixHandle(ds, nDate, true));	
	if (m_h->IsBasket()){	
		// Basket case
		std::vector<estring> aszNames;
		for (long n = 0; n < m_h->m_assets.size(); n++){
			aszNames.push_back(m_h->m_assets[n]->GetName());
		}	
		for (long n1 = 0; n1 < m_h->m_assets.size(); n1++){
			for (long n2 = n1 + 1; n2 < m_h->m_assets.size(); n2++){
				CSiriusApplication::InsertCorrelation(aszNames[n1], aszNames[n2], ds, nDate, false);
			}
		}
	} else {
		// Egg case
		// FX.FX correlations
		std::string szNatural, szComposite, szPay;		
		if (!!m_h->GetNaturalZeroCurve(false)) szNatural = CurrencyNameToCurrencyAssetName(m_h->GetNaturalZeroCurve(false)->GetName());
		if (!!m_h->GetCompositeZeroCurve(false)) szComposite = CurrencyNameToCurrencyAssetName(m_h->GetCompositeZeroCurve(false)->GetName());
		if (!!m_h->GetPayZeroCurve(false)) szPay = CurrencyNameToCurrencyAssetName(m_h->GetPayZeroCurve(false)->GetName());
		if (szNatural.size() && szComposite.size()) CSiriusApplication::InsertCorrelation(szNatural, szComposite, ds, nDate, false);
		if (szNatural.size() && szPay.size()) CSiriusApplication::InsertCorrelation(szNatural, szPay, ds, nDate, false);
		if (szPay.size() && szComposite.size()) CSiriusApplication::InsertCorrelation(szPay, szComposite, ds, nDate, false);

		// Asset.FX correlations
		if (!m_h->GetName().size()) return;
		if (szNatural.size()) CSiriusApplication::InsertCorrelation(m_h->GetName(), szNatural, ds, nDate, false);
		if (szComposite.size()) CSiriusApplication::InsertCorrelation(m_h->GetName(), szComposite, ds, nDate, false);
		if (szPay.size()) CSiriusApplication::InsertCorrelation(m_h->GetName(), szPay, ds, nDate, false);
	}
}


///////////////////////////////////////////////////////////////////////////////
//	SynchroniseEggCurrencies
//
//	If the composite currency is null then we set the composite currency pointer
//  in MlEqAsset to the natural currency. If the pay currency is null then we set
//  the pay currency pointer in MlEqAsset to the composite currency.
//	We don't synchronise the COM level pointers. Doing so would result in odd
//	behaviour if the user modified the natural currency more than once without
//	changing the other currencies explicitly.
//
//	(The Rekey function will check the presence of the necessary FX.FX correlations.)
//
void CAsset::SynchroniseEggCurrencies(const std::string& szOldKey)
{
	// We have to test nullness of the interface level pointers rather than the quant level ones
	// as the quant composite and pay pointers are always not null if the natural pointers are
	// not null.

	if (m_spzcEggNatural){
		if (m_spzcEggComposite && m_spzcEggPay){
			std::string szComposite, szPay;			
			szComposite = CurrencyNameToCurrencyAssetName(m_h->GetCompositeZeroCurve(true)->GetName());
			szPay = CurrencyNameToCurrencyAssetName(m_h->GetPayZeroCurve(true)->GetName());
		}
		if (!m_spzcEggComposite){
			m_h->PutCompositeZeroCurve(m_h->GetNaturalZeroCurve(true));
		} else if (m_spzcEggNatural){
			std::string szNatural, szComposite;
			szNatural = CurrencyNameToCurrencyAssetName(m_h->GetNaturalZeroCurve(true)->GetName());
			szComposite = CurrencyNameToCurrencyAssetName(m_h->GetCompositeZeroCurve(true)->GetName());
		}
		if (!m_spzcEggPay){
			m_h->PutPayZeroCurve(m_h->GetNaturalZeroCurve(true));
		} else if (m_spzcEggNatural){
			std::string szNatural, szPay;
			szNatural = CurrencyNameToCurrencyAssetName(m_h->GetNaturalZeroCurve(true)->GetName());
			szPay = CurrencyNameToCurrencyAssetName(m_h->GetPayZeroCurve(true)->GetName());
		}
	}
	if (m_spEggAssetNatural){		
		if (!m_spEggAssetPay)				m_h->PutFxPayAsset(m_h->GetFxNaturalAsset());
		if (!m_spEggAssetComposite)			m_h->PutFxCompositeAsset(m_h->GetFxNaturalAsset());
	}
	Rekey(szOldKey);
}