// product_basketcliquet.h : Declaration of the CBasketCliquet

#ifndef __BASKETCLIQUET_H_
#define __BASKETCLIQUET_H_

#include "resource.h"       // main symbols
#include "basketcliquetpricer.h"
#include "MlEqObjects.h"

/////////////////////////////////////////////////////////////////////////////
// CBasketCliquet
class ATL_NO_VTABLE CBasketCliquet : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CBasketCliquet, &CLSID_BasketCliquet>,
	public IDispatchImpl<IBasketCliquet, &IID_IBasketCliquet, &LIBID_Products>,
	public BasketCliquetPricer
{
public:
	HRESULT								FinalConstruct(void);
	CBasketCliquet(){}

	DECLARE_REGISTRY_RESOURCEID(IDR_BASKETCLIQUET)
	DECLARE_NOT_AGGREGATABLE(CBasketCliquet)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CBasketCliquet)
		COM_INTERFACE_ENTRY(IBasketCliquet)
		COM_INTERFACE_ENTRY(IDispatch)
	END_COM_MAP()
	
	DECLARE_MEMBER_VARIABLE(double, m_fGlobalCap, GlobalCap)
	DECLARE_MEMBER_VARIABLE(double, m_fGlobalFloor, GlobalFloor)
	DECLARE_MEMBER_VARIABLE(double, m_fGlobalRedemption, GlobalRedemption)
	DECLARE_MEMBER_VARIABLE(double, m_fGlobalGearing, GlobalGearing)
	DECLARE_MEMBER_VARIABLE(double, m_fNotional, Notional)
	DECLARE_MEMBER_OBJECT(Array, m_hCallPuts, CallPutArray);			
	DECLARE_MEMBER_OBJECT(Array, m_hLocalCaps, LocalCaps);
	DECLARE_MEMBER_OBJECT(Array, m_hLocalFloors, LocalFloors);
	DECLARE_MEMBER_OBJECT(Array, m_hStrikes, Strikes);	
	DECLARE_MEMBER_OBJECT(Array, m_hRainbowWeights, RainbowWeights);
	DECLARE_MEMBER_OBJECT(ParameterList, m_hParameterList, ModelInfo);
	DECLARE_MEMBER_OBJECT(MonteCarlo, m_hMonteCarlo, ForwardSkewHandle);
	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying);			
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hFixingSchedule, DateAndFixingSchedule);
	
public:		
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(Evaluate)(IResult* pVal);
};

#endif //__BASKETCLIQUET_H_
