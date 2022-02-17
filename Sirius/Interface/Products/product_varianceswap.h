// product_varianceswap.h : Declaration of the CVarianceSwapPricer

#ifndef __VARIANCESWAPPRODUCT_H_
#define __VARIANCESWAPPRODUCT_H_

#include "resource.h"     
#include "mleqobjects.h"


/////////////////////////////////////////////////////////////////////////////
// CVarianceSwapPricer
class ATL_NO_VTABLE CVarianceSwapPricer : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CVarianceSwapPricer, &CLSID_varianceswap>,
	public IDispatchImpl<IVarianceSwap, &IID_IVarianceSwap, &LIBID_Products>
{
public:
	CVarianceSwapPricer(){}
	HRESULT FinalConstruct(void);
	DECLARE_REGISTRY_RESOURCEID(IDR_VARIANCESWAP)
	DECLARE_NOT_AGGREGATABLE(CVarianceSwapPricer)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CVarianceSwapPricer)
		COM_INTERFACE_ENTRY(IVarianceSwap)
		COM_INTERFACE_ENTRY(IDispatch)
	END_COM_MAP()

	DECLARE_MEMBER_VARIABLE(double, m_upCutoff, UpCutoff)
	DECLARE_MEMBER_VARIABLE(double, m_downCutoff, DownCutoff)
	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying)
	DECLARE_MEMBER_VARIABLE(DATE, m_datePay, PayDate)
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hFixingSchedule, DateAndFixingSchedule)
	DECLARE_MEMBER_VARIABLE(PerformanceTypeEnum, m_perfType, PerformanceType)
	DECLARE_MEMBER_VARIABLE(StrikesTypeEnum, m_boundType, CutoffType)
	DECLARE_MEMBER_VARIABLE(YesNoEnum, m_useCurrentSpot, UseCurrentSpot)
	DECLARE_MEMBER_STRING(m_szCalendar, Calendar)
	DECLARE_MEMBER_VARIABLE(long, m_nDaysPerYear, NumberBusinessDaysPerYear)

public:	
	STDMETHOD(Evaluate)(IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);

//	STDMETHOD(EvaluateGreeks)(IResult* pVal);
};

#endif //__VARIANCESWAPPRODUCT_H_
