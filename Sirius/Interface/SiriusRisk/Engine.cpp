//	engine.cpp : Implementation of CEngine.
//				 This is the risk controller object.
//
//	Author :	 David Cuin
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SiriusRisk.h"
#include "engine.h"

/*static*/ const std::string CEngine::s_szRiskLibraries = "BasicRisk,SiriusDice";	// This is a comma separated list of risk scenario library names. ToDo - remove once we use COM categories

HRESULT CEngine::FinalConstruct(void)
{	
	return S_OK;
}

STDMETHODIMP CEngine::get_Aggregators(IAggregators** pVal)
{
	return m_spAggregators.CopyTo(pVal);	
}

estring CEngine::GetProgID(const CLSID& clsid) const
{
	CComBSTR							sProgID;
	CParameterMap::ProgIDFromCLSID(clsid, sProgID);
	return sProgID;
}

STDMETHODIMP CEngine::get_Portfolio(IEvaluatable** pVal)
{
	return m_portfolio.m_spPortfolio.CopyTo(pVal);
}

STDMETHODIMP CEngine::get_Scenario(IScenario** pVal)
{
	return m_spScenario.CopyTo(pVal);
}

STDMETHODIMP CEngine::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IEngine };
	for (int i = 0; i < sizeof(arr) / sizeof(arr[0]); i++) {
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CEngine::put_Aggregators(IAggregators* newVal)
{
	m_spAggregators = newVal;
	return S_OK;
}

STDMETHODIMP CEngine::put_Portfolio(IEvaluatable* newVal)
{
	begin_function
	portfolio							p;

	if (!(p.m_p = newVal)){
		// assignment to NULL is OK (corresponds to "SiriusRisk.Engine.Portfolio = Nothing" in VB)
		p.pt = unknown;
	} else if (p.m_pProduct = dynamic_cast<IProduct*>(newVal)){
		p.pt = product;
	} else if (p.m_pPosition = dynamic_cast<IPosition*>(newVal)){
		p.pt = position;
	} else if (p.m_pDeal = dynamic_cast<IDeal*>(newVal)){
		p.pt = deal;
	} else if (p.m_pProducts = dynamic_cast<IProducts*>(newVal)){
		p.pt = products;
	} else if (p.m_pPositions = dynamic_cast<IPositions*>(newVal)){
		p.pt = positions;
	} else if (p.m_pDeals = dynamic_cast<IDeals*>(newVal)){
		p.pt = deals;
	} else {
		std::string szObject;
		CParameterMap::GetObjectName(CComPtr<IDispatch>(newVal), &szObject);
		throw "'" + szObject + "' is an invalid portfolio object";
	}
	
	p.m_spPortfolio = newVal;
	m_portfolio = p;
	end_function
}

STDMETHODIMP CEngine::put_Scenario(IScenario* newVal)
{
	m_spScenario = newVal;
	return S_OK;
}


/////////////////////////////////////////////////////////////////////////////
//	Run
//
//	Implementation of the risk engine Run method.
//
STDMETHODIMP CEngine::Run(IResultsCollection** pVal)
{
	begin_function
	std::vector<scenario_details>				asdScenarios;			// Scenario object or creation details for a scenario object.
	std::vector<CAdapt<CComPtr<IParameter> > >	aspParameters;			// Set of parameters for each scenario.
					
	if (m_spScenario){
		// Resolve the scenario into the scenarios and parmaeters. Note that
		// the Sirius scenarios are a bit like sheep in that 'scenario' (as
		// in Sirius.Scenario is the collective name for scenario). Each line 
		// in the scenario object has an object name or an instantiated object
		// in the first column. The rest of the columns in the scenario (if
		// populated) comprise the arguments of the "Setup" function (if any)
		// contained within each scenario object.
		CComVariant		vScenario;
		CParameterMap	pmScenario;		

		m_spScenario->get_Value(&vScenario);
		pmScenario.SetValue(vScenario);
		aspParameters.resize(pmScenario.GetRows());

		for (long nRow = 0; nRow < pmScenario.GetRows(); nRow++){
			scenario_details	sd;
			if (!pmScenario.GetValue(nRow, 0, &sd.m_szName)){
				// ToDo - use COM categories				
				std::vector<std::string> asz;
				estring::Split(s_szRiskLibraries, ",", &asz);
				for (long n = 0; n < asz.size(); n++){
					CParameterMap::CLSIDFromProgID(estring(asz[n] + "." + sd.m_szName), sd.m_clsid);
					if (sd.m_clsid != CLSID_NULL) break;
				}		
				if (sd.m_clsid == CLSID_NULL) throw "The Sirius Risk Engine could not find scenario '" + sd.m_szName + "'";
			} else if (!pmScenario.GetValue(nRow, 0, sd.m_spObject)){
				CParameterMap::GetObjectName(sd.m_spObject, &sd.m_szName);
				sd.m_clsid = CLSID_NULL;
			}		
			asdScenarios.push_back(sd);
			
			// Assemble the arguments for the setup function.
			CParameterMap pmArguments;
			pmScenario.GetRow(nRow, &pmArguments);
			pmArguments.RemoveColumn(0);
			aspParameters[nRow].m_T.CoCreateInstance(CLSID_Parameter);
			aspParameters[nRow].m_T->put_Value(pmArguments.GetValue());
		}	
	} else {
		throw "No scenario object has been defined";
	}

	// Resolve the portfolio and run the scenario.
	if (!m_portfolio.m_spPortfolio){
		throw "You have not passed in a portfolio of products, positions or deals into the risk engine";
	} else {						
		// Create the vector of result collections (one per scenario).
		CComPtr<IResultsCollection> spResultsCollection;
		if (spResultsCollection.CoCreateInstance(CLSID_ResultsCollection)) throw "Failure when attempting to create the output results collection";		
		spResultsCollection->put_Size(aspParameters.size());
		switch (m_portfolio.pt){
		case (product):
			RunProduct(aspParameters, asdScenarios, 1.0, NULL, NULL, m_portfolio.m_pProduct, spResultsCollection);
			break;
		case (products):
			RunProducts(aspParameters, asdScenarios, 1.0, NULL, NULL, m_portfolio.m_pProducts, spResultsCollection);
			break;
		case (position):
			RunPosition(aspParameters, asdScenarios, 1.0, NULL, m_portfolio.m_pPosition, spResultsCollection);
			break;
		case (positions):
			RunPositions(aspParameters, asdScenarios, 1.0, NULL, m_portfolio.m_pPositions, spResultsCollection);
			break;
		case (deal):
			RunDeal(aspParameters, asdScenarios, 1.0, m_portfolio.m_pDeal, spResultsCollection);
			break;
		case (deals):
			RunDeals(aspParameters, asdScenarios, 1.0, m_portfolio.m_pDeals, spResultsCollection);
			break;
		default:
			// Should have already been disallowed.
			throw "Unhandled exception in CEngine::Run";
		}
		
		// Process Aggregation
		std::vector<scenario_details>::const_iterator itScenarios;
		CComVariant vResults = 0L;
		for (itScenarios = asdScenarios.begin(); itScenarios != asdScenarios.end(); itScenarios++){
			CComPtr<IResults>	spResults;
			HRESULT				hr;
			vResults.lVal++;
			spResultsCollection->get_Item(vResults, &spResults);			
			if (hr = CComQIPtr<IAggregatable>(m_spAggregators)->Aggregate(estring::GetBSTR(itScenarios->m_szName), spResults)) return hr;
		}

		// Return the results collection
		return spResultsCollection.CopyTo(pVal);		
	}
			
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	Run[...]
//
//	These functions deal with different collection inputs.
//  They all ultimately call into RunProduct.
//	It is my intention that these functions are NOT distributable.
//
void CEngine::RunProducts(const std::vector<CAdapt<CComPtr<IParameter> > >& aspParameters, const std::vector<scenario_details>& asdScenarios, double fNotional, CComPtr<IDeal> spDeal, CComPtr<IPosition> spPosition, CComPtr<IProducts> spProducts, CComPtr<IResultsCollection> spResultsCollection) const
{
	long nProducts;
	spProducts->get_Count(&nProducts);
	for (CComVariant vIndex = 1L; vIndex.lVal <= nProducts; vIndex.lVal++){
		CComPtr<IProduct>	spProduct;
		double				f;
		spProducts->get_Item(vIndex, &spProduct);
		spProducts->get_Notional(vIndex, &f);
		RunProduct(aspParameters, asdScenarios, fNotional * f, spDeal, spPosition, spProduct, spResultsCollection);
	}
}

void CEngine::RunPosition(const std::vector<CAdapt<CComPtr<IParameter> > >& aspParameters, const std::vector<scenario_details>& asdScenarios, double fNotional, CComPtr<IDeal> spDeal, CComPtr<IPosition> spPosition, CComPtr<IResultsCollection> spResultsCollection) const
{
	CComPtr<IProducts> spProducts;
	double f;	
	
	spPosition->get_Products(&spProducts);
	spPosition->get_Notional(&f);	
	RunProducts(aspParameters, asdScenarios, fNotional * f, spDeal, spPosition, spProducts, spResultsCollection);
}

void CEngine::RunPositions(const std::vector<CAdapt<CComPtr<IParameter> > >& aspParameters, const std::vector<scenario_details>& asdScenarios, double fNotional, CComPtr<IDeal> spDeal, CComPtr<IPositions> spPositions, CComPtr<IResultsCollection> spResultsCollection) const
{
	long nPositions;
	spPositions->get_Count(&nPositions);
	for (CComVariant vIndex = 1L; vIndex.lVal <= nPositions; vIndex.lVal++){
		CComPtr<IPosition>	spPosition;
		double				f;		
		spPositions->get_Item(vIndex, &spPosition);
		spPositions->get_Notional(vIndex, &f);		
		RunPosition(aspParameters, asdScenarios, fNotional * f, spDeal, spPosition, spResultsCollection);
	}
}

void CEngine::RunDeal(const std::vector<CAdapt<CComPtr<IParameter> > >& aspParameters, const std::vector<scenario_details>& asdScenarios, double fNotional, CComPtr<IDeal> spDeal, CComPtr<IResultsCollection> spResultsCollection) const
{
	CComPtr<IPositions> spPositions;
	double f;	
	
	spDeal->get_Positions(&spPositions);
	spDeal->get_Notional(&f);	
	RunPositions(aspParameters, asdScenarios, fNotional * f, spDeal, spPositions, spResultsCollection);
}

void CEngine::RunDeals(const std::vector<CAdapt<CComPtr<IParameter> > >& aspParameters, const std::vector<scenario_details>& asdScenarios, double fNotional, CComPtr<IDeals> spDeals, CComPtr<IResultsCollection> spResultsCollection) const
{
	long nDeals;
	spDeals->get_Count(&nDeals);
	for (CComVariant vIndex = 1L; vIndex.lVal <= nDeals; vIndex.lVal++){
		CComPtr<IDeal>		spDeal;
		spDeals->get_Item(vIndex, &spDeal);
		RunDeal(aspParameters, asdScenarios, fNotional, spDeal, spResultsCollection);
	}
}


/////////////////////////////////////////////////////////////////////////////
//	RunProduct
//
//	This is the only method that makes a call to the relevant scenario(s).
//  It is my intention that this function is distributable.
//
void CEngine::RunProduct(const std::vector<CAdapt<CComPtr<IParameter> > >& aspParameters, const std::vector<scenario_details>& asdScenarios, double fNotional, CComPtr<IDeal> spDeal, CComPtr<IPosition> spPosition, CComPtr<IProduct> spProduct, CComPtr<IResultsCollection> spResultsCollection) const
//	May have to pass aspParameters and vclsid by value for distribution (or this function is a stub).
{
	if (aspParameters.size() != asdScenarios.size()) throw "Error in CEngine::RunProduct. The parameter and scenario vectors are not the same size.";	

	std::vector<CAdapt<CComPtr<IParameter> > >::const_iterator		itParameter;	
	std::vector<scenario_details>::const_iterator					itScenario;
	long															nResults = 0L;
	long															nResult = 0L;

	spResultsCollection->get_Count(&nResults);	
	for (itScenario = asdScenarios.begin(), itParameter = aspParameters.begin(); itScenario != asdScenarios.end() && itParameter != aspParameters.end(); itScenario++, itParameter++){
		// Instantiate or get a scenario object
		CComPtr<IDispatch>	spScenario;
		
		if (itScenario->m_spObject){
			// Object already exists. THIS IS NON-DISTRIBUTABLE.
			spScenario = itScenario->m_spObject;
		} else {
			spScenario.CoCreateInstance(itScenario->m_clsid);
			if (!spScenario){
				throw "The Sirius Risk Engine could not create the scenario with name '" + itScenario->m_szName + "'";
			}
		}

		// Get the appropriate results object
		CComPtr<IResults>	spResults;
		spResultsCollection->get_Item(CComVariant(++nResult), &spResults);
		
		// Invoke the scenario's Setup method (if any).
		DISPID dispidSetup = -1L;
		if (CComDispatchDriverEx(spScenario).GetIDOfName(L"Setup", &dispidSetup)){
			// No setup method. The parameter must be blank in this case. Hopefully the parameter is blank.
			VARIANT_BOOL b = FALSE;
			if ((itParameter->m_T)->IsBlank(&b)) throw "Unhandled exception in CEngine::RunProduct";
			if (!b){
				throw "The scenario '" + itScenario->m_szName + "' does not require any input parameters";
			}
		} else {
			// Setup method found. In this case, it is the Setup method's responibility to detect whether
			// a blank aspParameters element is meaningful. I want to support blanks in general to enable
			// default arguments.			
			CComVariant av[1] = {itParameter->m_T.p};			
			if (CComDispatchDriverEx(spScenario).InvokeN(dispidSetup, av, 1, NULL)){
				propagate_error;
			}
		}

		// Invoke the object's Evaluate method to return the results collection for the object.
		DISPID dispidEvaluate = -1L;
		CComPtr<IResults> spCalcResults;
		if (CComDispatchDriverEx(spScenario).GetIDOfName(L"Evaluate", &dispidEvaluate)){
			throw "The scenario '" + itScenario->m_szName + "' does not implement the Evaluate method";
		} else {			
			if (spCalcResults.CoCreateInstance(CLSID_Results)) throw "Failure when attempting to create the results collection";
			CComVariant av[5] = {spCalcResults.p, spProduct.p, spPosition.p, spDeal.p, fNotional};
			if (CComDispatchDriverEx(spScenario).InvokeN(dispidEvaluate, av, 5, NULL)){
				propagate_error;
			}
			
			// Add spCalcResults to spResults
			long nCount = 0L;
			spCalcResults->get_Count(&nCount);
			for (CComVariant vItem = 1L; vItem.lVal <= nCount; vItem.lVal++){
				CComPtr<IResult> spResult;			
				if (spCalcResults->get_Item(vItem, &spResult)) ATLASSERT(false);
				if (spResults->Add(CComVariant(), spResult)) ATLASSERT(false);
			}
		}
	}
}
