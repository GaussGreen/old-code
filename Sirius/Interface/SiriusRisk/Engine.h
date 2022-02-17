//	Engine.h : Declaration of the CEngine
//
//////////////////////////////////////////////////////////////////////

#ifndef __ENGINE_H_
#define __ENGINE_H_

#include "resource.h"
#include "JetAggregator.h"

class ATL_NO_VTABLE CEngine : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CEngine, &CLSID_Engine>,
	public ISupportErrorInfo,
	public IDispatchImpl<IEngine, &IID_IEngine, &LIBID_SiriusRisk>
{
private:
	static const std::string			s_szRiskLibraries;				// This is a comma separated list of risk scenario library names. ToDo - remove once we use COM categories

protected:	
	enum portfolio_type {
		unknown = 0,
		product = 1,
		position = 2,
		deal = 3,
		products = 4,
		positions = 5,
		deals = 6
	};

	struct scenario_details {	
		std::string						m_szName;
		CComPtr<IDispatch>				m_spObject;			// may not be defined if the object is created at run-time		
		CLSID							m_clsid;			// CLSID of an object if we need to create it at run-time
	};

	struct portfolio {
		CComPtr<IEvaluatable>			m_spPortfolio;
		portfolio_type					pt;
		union {
			// Obviously we can't have CComPtr<> members in this union
			// since the former has a copy constructor.
			// m_spPortfolio will perform the necessary reference counting.
			// Alternatively, we could avoid the union, BUT given that the 
			// members are mutually exclusive, the union is the most natural
			// way of representing this.
			IDispatch*					m_p;
			IProduct*					m_pProduct;
			IPosition*					m_pPosition;
			IDeal*						m_pDeal;
			IProducts*					m_pProducts;
			IPositions*					m_pPositions;
			IDeals*						m_pDeals;
		};								
	};

public:
	CEngine(){}

	HRESULT								FinalConstruct();
	DECLARE_REGISTRY_RESOURCEID(IDR_ENGINE)
	DECLARE_NOT_AGGREGATABLE(CEngine)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CEngine)
		COM_INTERFACE_ENTRY(IEngine)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Scenario)(/*[out, retval]*/ IScenario** pVal);
	STDMETHOD(put_Scenario)(/*[in]*/ IScenario* newVal);	
	STDMETHOD(get_Portfolio)(/*[out, retval]*/ IEvaluatable** pVal);
	STDMETHOD(put_Portfolio)(/*[in]*/ IEvaluatable* newVal);		
	STDMETHOD(get_Aggregators)(/*[out, retval]*/ IAggregators** pVal);
	STDMETHOD(put_Aggregators)(/*[in]*/ IAggregators* newVal);	
	STDMETHOD(Run)(/*[out, retval]*/ IResultsCollection** pVal);	

protected:
	CComPtr<IScenario>					m_spScenario;
	portfolio							m_portfolio;							// this is the product-derived data on which the risk engine operates	
	CComPtr<IAggregators>				m_spAggregators;
				
	estring								GetProgID(const CLSID& clsid) const;
	void								RunProduct(const std::vector<CAdapt< CComPtr<IParameter> > >& aspParameters, const std::vector<scenario_details>& asdScenarios, double fNotional, CComPtr<IDeal> spDeal, CComPtr<IPosition> spPosition, CComPtr<IProduct> spProduct, CComPtr<IResultsCollection> spResultsCollection) const;
	void								RunProducts(const std::vector<CAdapt< CComPtr<IParameter> > >& aspParameters, const std::vector<scenario_details>& asdScenarios, double fNotional, CComPtr<IDeal> spDeal, CComPtr<IPosition> spPosition, CComPtr<IProducts> spProducts, CComPtr<IResultsCollection> spResultsCollection) const;
	void								RunPosition(const std::vector<CAdapt< CComPtr<IParameter> > >& aspParameters, const std::vector<scenario_details>& asdScenarios, double fNotional, CComPtr<IDeal> spDeal, CComPtr<IPosition> spPosition, CComPtr<IResultsCollection> spResultsCollection) const;
	void								RunPositions(const std::vector<CAdapt< CComPtr<IParameter> > >& aspParameters, const std::vector<scenario_details>& asdScenarios, double fNotional, CComPtr<IDeal> spDeal, CComPtr<IPositions> spPositions, CComPtr<IResultsCollection> spResultsCollection) const;
	void								RunDeal(const std::vector<CAdapt<CComPtr<IParameter> > >&, const std::vector<scenario_details>& asdScenarios, double fNotional, CComPtr<IDeal> spDeal, CComPtr<IResultsCollection> spResultsCollection) const;
	void								RunDeals(const std::vector<CAdapt<CComPtr<IParameter> > >&, const std::vector<scenario_details>& asdScenarios, double fNotional, CComPtr<IDeals> spDeals, CComPtr<IResultsCollection> spResultsCollection) const;	
};

#endif
