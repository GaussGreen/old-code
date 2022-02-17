//	product_europeanoption.h : Declaration of the CEuropeanOption
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __PRODUCTEUROPEANOPTION_H_
#define __PRODUCTEUROPEANOPTION_H_

#include "resource.h"
#include "mleqobjects.h"


class EuropeanBasketProduct	:	public product
{
protected:
	double		m_fStrike;
	double		m_cp;
	int			m_nMaturity ;
	CVector		m_weight;

public:
	EuropeanBasketProduct(){};
	~EuropeanBasketProduct(){};

	void initialize(DATE dateStart, DATE dateMaturity, double fStrike, double cp, NotionalTypeEnum nte, MlEqAssetHandle udly, CMatrix& fixings);
					
	void payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc);
	void setUp(CMatrix& value,MlEqMonteCarlo& mc);
};


class ATL_NO_VTABLE CEuropeanOption : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CEuropeanOption, &CLSID_EuropeanOption>,
	public ISupportErrorInfo,
	public IDispatchImpl<IEuropeanOption, &IID_IEuropeanOption, &LIBID_Products>,
	public EuropeanBasketProduct
{
public:
	HRESULT FinalConstruct(void);
	DECLARE_REGISTRY_RESOURCEID(IDR_EUROPEANOPTION)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CEuropeanOption)
		COM_INTERFACE_ENTRY(IEuropeanOption)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	CEuropeanOption(){}

	// member variables with COM wrappers	
	DECLARE_MEMBER_VARIABLE(double, m_fStrike, Strike)
	DECLARE_MEMBER_VARIABLE(DATE, m_dateMaturity, Maturity)
	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying)
	DECLARE_MEMBER_VARIABLE(DATE, m_datePay, PayDate)
	DECLARE_MEMBER_VARIABLE(DATE, m_dateStart, StartDate)
	DECLARE_MEMBER_VARIABLE(double, m_fNotional, Notional)
	DECLARE_MEMBER_VARIABLE(NotionalTypeEnum, m_nte, NotionalType)
	DECLARE_MEMBER_VARIABLE(PayoffTypeEnum, m_cpe, CallOrPut)
	DECLARE_MEMBER_OBJECT(ParameterList, m_hParameterList, ModelInfo)

public:	
	STDMETHOD(Evaluate)(IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);

protected:
	HRESULT			EvaluateBlackScholes(long nToday, MlEqDateHandle hDate, double fDiscountFactor, NotionalTypeEnum nte, double fCallPut, double fStrike, double fNotional, IResult** ppVal);
	HRESULT			EvaluateMonteCarlo(long nToday, MlEqDateHandle hDate, double fDiscountFactor, NotionalTypeEnum nte, double fCallPut, double fStrike, double fNotional, IResult** ppVal);
	HRESULT			EvaluateMonteCarlo(IResult** ppVal);
};




#endif
