//	handle_asset.h : Declaration of the CAsset
//
//	Author:			 David Cuin
//
//  Despite trying, I couldn't recall or discover a sensible opposite of a basket
//  asset. 'Single Underlying', though acceptable to the COM interface is too
//  much of a mouthful, and beside, another type of non-basket is a 'Natural Asset'.
//  I have therefore decided to use the term 'egg' for an asset that is neither a
//  currency asset nor a basket. My apologies.
//
//	The three basic asset types are therefore.
//		1)			Currency Asset	(e.g. GBP.USD)
//		2)			Egg Asset	    (e.gs. .FTSE or .FTSE > CHF)
//		3)			Basket Asset	(collection of (1) and (2))
//
/////////////////////////////////////////////////////////////////////////////


#ifndef __ASSET_H_
#define __ASSET_H_

#include "resource.h"
#include "MlEqAsset.h"

class ATL_NO_VTABLE CAsset : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CAsset, &CLSID_Asset>,
	public ISupportErrorInfo,
	public IDispatchImpl<IAsset, &IID_IAsset, &LIBID_Sirius>
{
public:
	CAsset() : m_h(new MlEqAsset) {}
	HRESULT FinalConstruct(void);
	HRESULT FinalRelease(void);
	DECLARE_REGISTRY_RESOURCEID(IDR_ASSET)
	DECLARE_NOT_AGGREGATABLE(CAsset)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CAsset)
		COM_INTERFACE_ENTRY(IAsset)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()			
	associate_analytic_object(CAsset, MlEqAsset);
	declare_serialisable;
	static HRESULT						Load(const std::string& szIdentifier, DataSourceEnum ds, long nDate, CComPtr<IAsset>& spAsset);
		
	STDMETHOD(get_AssetType)(AssetTypeEnum* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);	
	STDMETHOD(get_CompositeCurrency)(/*[out, retval]*/ CurrencyEnum* pVal);	
	STDMETHOD(get_Identifier)(/*[out, retval]*/ BSTR* pVal);		
	STDMETHOD(get_Currency)(/*[out, retval]*/ CurrencyEnum* pVal);
	STDMETHOD(put_Identifier)(/*[in]*/ BSTR newVal);
	STDMETHOD(get_PayCurrency)(/*[out, retval]*/ CurrencyEnum* pVal);
	STDMETHOD(GetZeroCurve)(/*[in]*/ ZeroCurveTypeEnum ZeroCurveType, /*[out, retval]*/ IZeroCurve** pVal);			
	STDMETHOD(get_SpotSchedule)(/*[out, retval]*/ ISpotSchedule** pVal);
	STDMETHOD(put_SpotSchedule)(/*[in]*/ ISpotSchedule* pVal);
	STDMETHOD(get_DividendSchedule)(/*[out, retval]*/ IDividendSchedule** pVal);
	STDMETHOD(put_DividendSchedule)(/*[in]*/ IDividendSchedule* pVal);
	STDMETHOD(get_Assets)(/*[out, retval]*/ IAssets** pVal);
	STDMETHOD(put_Assets)(/*[in]*/ IAssets* newVal);
	STDMETHOD(get_VolatilityStructure)(/*[out, retval]*/ IVolatilityStructure** pVal);
	STDMETHOD(put_VolatilityStructure)(/*[in]*/ IVolatilityStructure* pVal);
	STDMETHOD(get_ConstituentWeight)(/*[in]*/ IAsset* Constituent, /*[out, retval]*/ double *pVal);
	STDMETHOD(put_ConstituentWeight)(/*[in]*/ IAsset* Constituent, /*[in]*/ double newVal);		
	STDMETHOD(get_Name)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_IsCalibrated)(/*[out, retval]*/ VARIANT_BOOL *pVal);
	STDMETHOD(GetSpots)(/*[in]*/ DATE Today, /*[in]*/ IDateSchedule* Dates, /*[out, retval]*/ ISpotSchedule** pVal);
	STDMETHOD(Save)(/*[out, retval]*/BSTR* pVal);
	STDMETHOD(GetForward)(/*[in]*/ DATE Today, /*[in]*/ DATE Maturity, /*[out, retval]*/ double* pVal);	
	STDMETHOD(get_PayCurrencyAsset)(/*[out, retval]*/ IAsset** pVal);
	STDMETHOD(put_PayCurrencyAsset)(/*[in]*/ IAsset* pVal);
	STDMETHOD(get_CompositeCurrencyAsset)(/*[out, retval]*/ IAsset** pVal);
	STDMETHOD(put_CompositeCurrencyAsset)(/*[in]*/ IAsset* pVal);
	STDMETHOD(get_CurrencyAsset)(/*[out, retval]*/ IAsset** pVal);
	STDMETHOD(put_CurrencyAsset)(/*[in]*/ IAsset* pVal);
	STDMETHOD(put_NaturalSpot)(/*[in]*/ DATE Date, /*[in]*/ double newVal);
	STDMETHOD(get_NaturalSpot)(/*[in]*/ DATE Date, /*[out, retval]*/ double* pVal);
	STDMETHOD(get_Spot)(/*[in]*/ DATE Date, /*[out, retval]*/ double* pVal); // The lack of a put_Spot is intentional - you don't want to modify the spot schedule of a quanto asset.
	STDMETHOD(GetBaseUnderlyings)(/*[out, retval]*/ IAssets** pVal);
	STDMETHOD(get_RealAssets)(/*[out, retval]*/ IAssets** pVal);
	STDMETHOD(Refresh)(/*[in, defaultvalue(0L)]*/ VARIANT_BOOL Recursive);
	STDMETHOD(IsBasket)(/*[out, retval]*/ VARIANT_BOOL* pVal);
	STDMETHOD(get_CurrencyCurve)(/*[out, retval]*/ IZeroCurve** pVal);
	STDMETHOD(put_CurrencyCurve)(/*[in]*/ IZeroCurve* newVal);
	
	static void										CreateBasket(VARIANT Identifier, VARIANT DataSourceOpt, VARIANT DateOpt, VARIANT AssetsArr, VARIANT CompositeCurrenciesArr, VARIANT WeightsArr, VARIANT PayCurrency, CComPtr<IAsset>& spAsset);	
	static void										CloneAsset(CComPtr<IAsset> spParentAsset, CurrencyEnum ceCloneComposite, CurrencyEnum ceClonePay, CComPtr<IAsset>& spNewAsset);

protected:
	// If we add any other items here then we need to adjust get_AssetType	
	declare_member_variable_rekey(DataSourceEnum, DataSource);
	declare_member_variable_rekey(DATE, Date);

protected:
	static CComPtr<IAsset>							CurrencyEnumToCurrencyAsset(CurrencyEnum ce, long nDate, DataSourceEnum ds);
	static std::string								CurrencyEnumToCurrencyAseetName(CurrencyEnum ce);
	static std::string								CurrencyNameToCurrencyAssetName(const std::string& szName);			
	AssetTypeEnum									GetAssetType(void) const;
	CurrencyEnum									GetCompositeCurrency(void) const;
	CurrencyEnum									GetCurrency(MlEqZeroCurveHandle hZeroCurve) const;		
	estring											GetIdentifierFromCurrencyCurve(void) const;
	CurrencyEnum									GetPayCurrency(void) const;
	void											GenerateRealConstituents(void);
	void											GenerateRealConstituents(CurrencyEnum ceCompositeIn, CurrencyEnum cePayIn);
	estring											GetName(void) const;
	void											PutEggCurrency(IAsset* newVal, CComPtr<IZeroCurve>& spZeroCurve, MlEqZeroCurveHandle& hZeroCurve, CComPtr<IAsset>& spAsset, MlEqAssetHandle& hAsset) const;
	void											Rekey(const std::string& szOldKey);
	void											SafeGetDateAndDataSource(long* pnDate, DataSourceEnum* pds) const;
	void											SetAnalyticToEgg(bool bDeclareAsBasket);
	void											SetCorrelations(void);	
	void											SynchroniseEggCurrencies(const std::string& szOldKey);

protected:	
	// Member variables for Baskets, Eggs and Currency assets.
	CComPtr<ISpotSchedule>							m_spSpotSchedule;
	CComPtr<IDividendSchedule>						m_spDividendSchedule;	
	CComPtr<IVolatilityStructure>					m_spVolatilityStructure;

	// Currency asset properties
	CComPtr<IZeroCurve>								m_spCurrencyCurve;

	// Egg properties
	CComPtr<IZeroCurve>								m_spzcEggNatural;
 	CComPtr<IZeroCurve>								m_spzcEggComposite;
	CComPtr<IZeroCurve>								m_spzcEggPay;
	CComPtr<IAsset>									m_spEggAssetNatural;
	CComPtr<IAsset>									m_spEggAssetComposite;
	CComPtr<IAsset>									m_spEggAssetPay;

	// Basket properties		
	CComPtr<IAssets>								m_spBasketAssets;
	CComPtr<IAssets>								m_spBasketRealAssets;			// m_spBasketAssets with the pay currencies adjusted - the analytic MlEqAsset uses these.
};

#endif
