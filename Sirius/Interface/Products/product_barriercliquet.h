//	product_barriercliquet.h : Declaration of the CBarrierCliquet
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __BARRIERCLIQUET_H_
#define __BARRIERCLIQUET_H_

#include "resource.h"
#include "MlEqObjects.h"
#include "MlEqBarrierCliquetPricer.h"

class ATL_NO_VTABLE CBarrierCliquet : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CBarrierCliquet, &CLSID_BarrierCliquet>,
	public ISupportErrorInfo,
	public IDispatchImpl<IBarrierCliquet, &IID_IBarrierCliquet, &LIBID_Products>,
	public MlEqBarrierCliquetPricer
{
public:	
	HRESULT FinalConstruct();
	CBarrierCliquet(){}
	DECLARE_REGISTRY_RESOURCEID(IDR_BARRIERCLIQUET)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CBarrierCliquet)
		COM_INTERFACE_ENTRY(IBarrierCliquet)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()

	DECLARE_MEMBER_OBJECT(DateSchedule, m_hBarrierStartDates,BarrierStartDates )
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hBarrierEndDates, BarrierEndDates)
	DECLARE_MEMBER_OBJECT(Array, m_hBarrierPayDates, BarrierPayDatesArr)
	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying)
	DECLARE_MEMBER_BOOL(m_bDelayPaymentsToEnd, DelayPaymentsToEnd)
	DECLARE_MEMBER_OBJECT(ParameterList, m_hParameterList, ModelInfo)
	DECLARE_MEMBER_OBJECT(Matrix, m_hBarriers, Barriers)
	DECLARE_MEMBER_OBJECT(Matrix, m_hRebate, Rebates)
	DECLARE_MEMBER_OBJECT(Matrix, m_hCallsPuts, CallsPuts)
	DECLARE_MEMBER_OBJECT(Matrix, m_hStrikes, Strikes)
	DECLARE_MEMBER_OBJECT(Matrix, m_hBarrierHasHit, BarrierHasHit)
	DECLARE_MEMBER_OBJECT(Array, m_hBarrierFixings, BarrierFixings)
	DECLARE_MEMBER_BOOL(m_bKnockIn, KnockIn)

public:
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(Evaluate)(IResult* pVal);
};

#endif
