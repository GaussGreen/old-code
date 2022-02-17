//	handle_volatilitystructure.h : Declaration of the CVolatilityStructure
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __VOLATILITYSTRUCTURE_H_
#define __VOLATILITYSTRUCTURE_H_

#include "resource.h"
#include "MlEqObjects.h"

class ATL_NO_VTABLE CVolatilityStructure : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CVolatilityStructure, &CLSID_VolatilityStructure>,
	public ISupportErrorInfo,
	public IDispatchImpl<IVolatilityStructure, &IID_IVolatilityStructure, &LIBID_Sirius>	
{
public:	
	CVolatilityStructure(){}
	DECLARE_REGISTRY_RESOURCEID(IDR_VOLATILITYSTRUCTURE)
	DECLARE_NOT_AGGREGATABLE(CVolatilityStructure)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CVolatilityStructure)
		COM_INTERFACE_ENTRY(IVolatilityStructure)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CVolatilityStructure, MlEqVolatilityStructure);
	declare_serialisable;
	static HRESULT							CVolatilityStructure::Load(const std::string& szIdentifier, DataSourceEnum ds, long nDate, CComPtr<IVolatilityStructure>& spVolatilityStructure);

	
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Name)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);			
	
	STDMETHOD(GetVolatility)(/*[in]*/ VARIANT Strikes, /*[in]*/ DATE Maturity, /*[in]*/ BidAskMidEnum BidAskMid, /*[out, retval]*/ VARIANT* pResult);
	STDMETHOD(get_Parameter)(/*[in]*/ BSTR Name, /*[in, optional]*/ VARIANT Index1Opt, /*[in, optional]*/ VARIANT Index2Opt, /*[in, optional]*/ VARIANT Index3Opt, /*[out, retval]*/ VARIANT *pVal);	
	STDMETHOD(Save)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(get_InterpolateType)(/*[out, retval]*/ InterpolateTypeEnum* pVal);
	STDMETHOD(get_VolatilityDataType)(/*[out, retval]*/ VolatilityDataTypeEnum* pVal);	
	STDMETHOD(get_AsymptoticVolatilityStructure)(/*[out, retval]*/ IVolatilityStructure** pVal);
	STDMETHOD(get_DateHandle)(/*[out, retval]*/ IDate** pVal);	
	STDMETHOD(get_Location)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(put_Location)(/*[in]*/ BSTR newVal);	
	STDMETHOD(get_Dates)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(get_VolatilityData)(/*[in]*/ BidAskMidEnum BidAskMid, /*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(get_Spot)(/*[out, retval]*/ double* pVal);
	STDMETHOD(get_Forwards)(/*[out, retval]*/ IArray** pVal);
	STDMETHOD(get_DiscountFactors)(/*[out, retval]*/ IArray** pVal);
	STDMETHOD(get_VolatilityDataName)(/*[in]*/ BidAskMidEnum BidAskMid, /*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(ReinitialiseStrike)(/*[in]*/ IAsset* Asset);
		
	// Scenario Functions
	STDMETHOD(Reset)(void);
	STDMETHOD(SetVol)(/*[in]*/ double vol);
	STDMETHOD(Shift)(/*[in]*/ double Amount, VARIANT_BOOL IsRelative);	
	STDMETHOD(ShiftSkew)(/*[in]*/ IStrikes* Strikes, /*[in]*/ IArray* Maturities,/*[in]*/ IMatrix* ShiftAcrossMaturity, IAsset* Asset);
	STDMETHOD(ShiftStrikeInvariantSkew)(IAsset* Asset, IInterpolator* VolShiftInterpolator, IInterpolator* VolSkewShiftInterpolator, VARIANT_BOOL IsRelative);
	STDMETHOD(ShiftStrikeSkew)(IAsset* Asset, IInterpolator* VolShiftInterpolator, IInterpolator* VolSkewShiftInterpolator);	
	STDMETHOD(Stick)(void);
	STDMETHOD(ShiftBeta)(double Amount, BSTR Tenor);
	STDMETHOD(ShiftForwardSkew)(/*[in]*/ double shiftAmount, /*[in]*/ DATE WindowStartDate, /*[in]*/ DATE WindowEndDate, /*[in]*/ long PeriodStartDays, /*[in]*/ long PeriodEndDays);
	STDMETHOD(FsaLiquidityShift)(/*[in]*/ double AverageVolume, /*[in]*/ double CoreVolatility, /*[in]*/ double NetGamma, /*[in]*/ double PositionDirection, /*[in]*/ VARIANT_BOOL IsStock);
protected:	
	declare_member_variable_rekey(DataSourceEnum, DataSource);
	declare_member_variable_rekey(DATE, Date);

protected:			
	long									FindSchedule(const estring& szName) const;
	estring									GetLocation(void) const;
	CParameterMap							GetVolatilityData(BidAskMidEnum BidAskMid) const;
	void									Make2DInterpolator(const std::map<long, std::map<long, double> > map, MlEq2DInterpolatorHandle& hOut) const;

protected:
	CComPtr<IDate>									m_spDateHandle;
	CComPtr<IInterpolator>							m_spInterpolator;
	std::vector<CAdapt<CComQIPtr<IDispatch> > >		m_aspBid;
	std::vector<CAdapt<CComQIPtr<IDispatch> > >		m_aspMid;
	std::vector<CAdapt<CComQIPtr<IDispatch> > >		m_aspAsk;
	std::vector<CAdapt<CComQIPtr<IDateSchedule> > >	m_aspSchedules;
	CComPtr<IVolatilityStructure>					m_spAsymptoticVolatilityStructure;
	estring											m_szName;
	estring											m_szLocation;
	std::vector<CParameterMap>						m_vpm;
	CComPtr<IInterpolator>							m_spExcessVolVolInterpolator;
	CComPtr<IInterpolator>							m_spElasticityInterpolator;
};

#endif
	