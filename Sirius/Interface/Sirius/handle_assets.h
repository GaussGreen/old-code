//	handle_assets.h : Declaration of the CAssets
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __ASSETS_H_
#define __ASSETS_H_

#include "resource.h"								// main symbols
#include "comobjectcollectionserialisable.h"		// template collection class
#include "handle_asset.h"

class ATL_NO_VTABLE CAssets : public CComObjectCollectionSerialisable<IAsset, CAssets, IAssets, &CLSID_Assets, &IID_IAssets, &LIBID_Sirius, &CLSID_Asset, IDR_ASSETS>
{
protected:
	typedef std::multimap<std::string, std::string>			IdentifierKeyToKeyType;
	IdentifierKeyToKeyType									m_IdentifierKeyToKey;				// multimap of A@Date.DS to A (B) > C@Date.DS

	// Virtual overrides
	HRESULT										GetClone(const std::string& szKey, CComPtr<IDispatch>& spObject);
	bool										IsClonable(void) const;	
	void										KeyAdded(const std::string& szKey);
	void										KeyCleared(const std::string& szKey);
	void										KeyRekeyed(const std::string& szKeyOld, const std::string& szKeyNew);
	void										KeysAllCleared(void);
	HRESULT										MapIndex(CComVariant* pIndex);
	bool										ShouldAddToMaintainedCollection(CComPtr<IDispatch> spObject);

public:	
	STDMETHOD(GetVolatilityStructures)(/*[out, retval]*/ IVolatilityStructures** pVal);
	STDMETHOD(GetBaseUnderlyings)(/*[out, retval]*/ IAssets** pVal);
	STDMETHOD(Refresh)(/*[in, defaultvalue(0L)]*/ VARIANT_BOOL Recursive);
	
	static void									DontAddToMaintainedCollection(CComPtr<IAsset> spAsset);
	static void									InsertCorrelations(CComPtr<IAssets> spUnderlyings);
	static void									InsertCorrelations(std::vector<MlEqAssetHandle> ahAssets);
	static HRESULT								Load(const std::string& szDummy, DataSourceEnum ds, long nDate, CComPtr<IAssets>& spAssets);
	static void									RemoveAssetFromDontAddSet(IAsset* pAsset);
	static void									Respell(std::string* pszIndex);
	
protected:			
	void										GetBaseUnderlyings(CComPtr<IAssets>& spAssets);
	static std::string							GetIdentifierKey(const std::string& szKey);
	static void									InsertCorrelations(const std::vector<estring>& aszNames, DataSourceEnum ds, long nDate);
	
	// The set below holds assets that we don't want to add to a maintained collection.
	// (Examples include assets that have been cloned from others - we want the cloning
	// process to be repeated).
	//
	// We don't want CComPtr<> here since then the mecahnism I use in CAsset::FinalRelease
	// mechanism could not be used to remove items from CAsset. (Reason is simple - the
	// FinalRelease would never be called since a CComPtr<> set would hold reference
	// counts indefinitely.)	
	
	static std::set<IAsset*>					s_setDontAdd;
};

#endif