// product_iguana.h : Declaration of the CIguana

#ifndef __IGUANA_H_
#define __IGUANA_H_

#include "resource.h"       // main symbols
#include "mleqobjects.h"
#include "IguanaProduct.h"


// CIguana
class ATL_NO_VTABLE CIguana : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CIguana, &CLSID_Iguana>,
	public ISupportErrorInfo,
	public IDispatchImpl<IIguana, &IID_IIguana, &LIBID_TemporaryProductsLib>
	,	public IguanaPricer
{
public:
	CIguana(){}
	HRESULT FinalConstruct();
	DECLARE_REGISTRY_RESOURCEID(IDR_IGUANA)
	DECLARE_NOT_AGGREGATABLE(CIguana)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CIguana)
		COM_INTERFACE_ENTRY(IIguana)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()

// ISupportsErrorInfo
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);

// IIguana
public:
	
	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying)
	DECLARE_MEMBER_VARIABLE(DATE, m_dateMaturity, Maturity)
	DECLARE_MEMBER_VARIABLE(double, m_fHighStrike, HighStrike)
	DECLARE_MEMBER_VARIABLE(double, m_fLowStrike, LowStrike)
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hFixingSchedule, DateAndFixingSchedule )
	DECLARE_MEMBER_OBJECT(MonteCarlo, m_hMonteCarlo, ForwardSkewHandle)
	DECLARE_MEMBER_OBJECT(ParameterList, m_hParameterList, ModelInfo)
		
	STDMETHOD(Evaluate)(IResult* pVal);

protected:
	double BSFormula(double k, double cp);
	double impliedVol(double k);
};



#endif //__IGUANA_H_
