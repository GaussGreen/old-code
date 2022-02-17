//	handle_marketdata.h :	Declaration of the CMarketData. The purpose of this class
//							is to hold all the market data that can be used by any
//							pricing model. The application object has exactly one
//							market data object attached to it.
//							
//							The user can create as many market data objects as they
//							please and switch them into the application object.
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __MARKETDATA_H_
#define __MARKETDATA_H_

#include "resource.h"
#include "handle_assets.h"

class ATL_NO_VTABLE CMarketData : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CMarketData, &CLSID_MarketData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IMarketData, &IID_IMarketData, &LIBID_Sirius>
{
public:	
	CMarketData(){}
	HRESULT								FinalConstruct(void);	

	DECLARE_REGISTRY_RESOURCEID(IDR_MARKETDATA)
	DECLARE_NOT_AGGREGATABLE(CMarketData)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CMarketData)
		COM_INTERFACE_ENTRY(IMarketData)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);

protected:
	DECLARE_MEMBER_COM_OBJECT(IAssets, m_spAssets, Assets);
	DECLARE_MEMBER_COM_OBJECT_RO(ICorrelationMatrices, m_spCorrelationMatrices, CorrelationMatrices);
	DECLARE_MEMBER_COM_OBJECT_RO(IDividendSchedules, m_spDividendSchedules, DividendSchedules);
	DECLARE_MEMBER_COM_OBJECT_RO(ISpotSchedules, m_spSpotSchedules, SpotSchedules);
	DECLARE_MEMBER_COM_OBJECT_RO(IVolatilityStructures, m_spVolatilityStructures, VolatilityStructures);
	DECLARE_MEMBER_COM_OBJECT_RO(IZeroCurves, m_spZeroCurves, ZeroCurves);				
};

#endif