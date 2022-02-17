// product_basketcliquet.cpp : Implementation of CBasketCliquet
#include "stdafx.h"
#include "Products.h"
#include "product_basketcliquet.h"
#include "MlEqObjects.h"
#include "MlEqParameterList.h"
#include "MonteCarlo.h"
#include "mleqpde.h"

#undef max
#undef min


HRESULT CBasketCliquet::FinalConstruct(void)
{	
	m_fGlobalCap = std::numeric_limits<double>::max();
    m_fGlobalFloor = std::numeric_limits<double>::min();
    m_fGlobalRedemption = 0.0;
    m_fGlobalGearing = 1.0;
	m_fNotional = 1.0;
	return S_OK;
}

STDMETHODIMP CBasketCliquet::Evaluate(IResult* pVal)
{
	begin_function
			
	if (!m_hLocalCaps) throw "No local caps defined";
	if (!m_hLocalFloors) throw "No local floors defined";
	if (!m_hStrikes) throw "No strikes defined";
	if (!m_hCallPuts) throw "No call / put schedule defined";
	if (!m_hRainbowWeights) throw "No weights defined";
	if (!m_hFixingSchedule) throw "No fixing schedule defined";

	if( m_hRainbowWeights->getsize() != m_hUnderlying->m_assets.size() )throw "Wrong weights size";

	HRESULT								hr;
	MlEqMonteCarloHandle				hMonteCarlo;
	CMatrix								mResults;


	// We deal with two cases:
	// (1) If the forward skew handle is given then the fixing schedule must contain fixing values. 
	//	   The asset handle must not be given.
	//	   The model info parameter must not be given.
	// (2) If the asset handle is given, then the fixing schedule can contain the fixing values, but might not (in which case we get them from the aseet).
	//	   The forward skew handle should not be given.
	//	   The model info parameter must be given.	
	if (!m_hUnderlying && !m_hMonteCarlo){
		throw "Either a value for parameter 'ForwardSkewHandle' or 'Underlying' must be supplied";
	} else if (!m_hUnderlying){
		// m_hMonteCarlo passed in. The current date is set from MlEqDate (since it cannot be obtained from an asset).
		GVector<long>					afDates;
		CMatrix							mFixings;		
		hMonteCarlo = m_hMonteCarlo;					
		if (m_hFixingSchedule->GetWidth() != 1){			
			// Error if there are dates in the past
			if (m_hFixingSchedule->GetFirstDate() < MlEqDate::GetCurrentDate()){			
				throw "Invalid fixing schedule. It should be a two column matrix of dates and values.";
			}
		}
		m_hFixingSchedule->GetDates(afDates);
		m_hFixingSchedule->GetColumnBeforeToday(MlEqDate::GetCurrentDate() /*This is OK in this instance*/, 0, mFixings);	
	
		initialize(*m_hLocalCaps, *m_hLocalFloors, *m_hStrikes, *m_hCallPuts, m_fGlobalCap, m_fGlobalFloor, m_fGlobalRedemption, m_fGlobalGearing, m_fNotional, afDates, *m_hRainbowWeights );
		
		hMonteCarlo->initializeMcHelper(*this, mFixings);
		hMonteCarlo->simulate(mResults, *this);
	} else if (!m_hMonteCarlo) {
		// an asset handle has been passed in
		CVector										modelParameters(16);
		GVector<long>								afDates;
				
		if (!m_hParameterList) throw "No model info parameter list has been defined";
		m_hParameterList->BeginExhaust();
		m_hParameterList->GetValue("cL, cR, addTanhWing, yPower, seed, npaths, calibflag, numberVolStates, localVolFlag, savevolvolinfo, numberGridPoints, saveVolGridFlag, randomNumberFlag, ControlVariate, GlobalCalib, UseHermite = 0", ",", "=", &modelParameters[0]);
		m_hParameterList->CheckUnexhausted();
		m_hFixingSchedule->GetDates(afDates);
								
		// there are two cases here - (i) basket case, (ii) non-basket case
		if (m_hUnderlying->IsBasket()){
			// basket case
			CMatrix									mFixings;
			long									nWidth = m_hFixingSchedule->GetWidth();

			// initialise the CappedFlooredCliquetPricer
			initialize(*m_hLocalCaps, *m_hLocalFloors, *m_hStrikes, *m_hCallPuts, m_fGlobalCap, m_fGlobalFloor, m_fGlobalRedemption, m_fGlobalGearing, m_fNotional, afDates, *m_hRainbowWeights );
		
			// set up the fixings
			if (!nWidth){
				// fixing schedule contains just dates - we get the values from the basket constituent spot schedules
				std::vector<long>			afDates;
				m_hFixingSchedule->GetDates(afDates);												
				mFixings.resize(m_hUnderlying->m_assets.size(), m_hFixingSchedule->size(), 0.0);				
				for (long nAsset = 0; nAsset < m_hUnderlying->m_assets.size(); nAsset++){			
					std::vector<double>	vSpots;
					m_hUnderlying->m_assets[nAsset]->GetSpotsBeforeToday(m_hUnderlying->GetDateHandle()->GetDate(), afDates, &vSpots);
					for (long nSpot = 0; nSpot < vSpots.size(); nSpot++){
						mFixings[(int)nAsset][(int)nSpot] = vSpots[nSpot];
					}
				}
			} else if (nWidth == m_hUnderlying->m_assets.size()){
				// fixing schedule contains both the dates and the fixings
				m_hFixingSchedule->GetColumnsBeforeToday(m_hUnderlying->GetDateHandle()->GetDate(), mFixings);
			} else {
				throw "The number of columns in the fixing schedule is invalid - you need to pass as many data columns as there are assets in the basket";
			}
			
			// Initialise Monte Carlo
			CMultiAssetForwardSkewMC* pMC;			
			hMonteCarlo = pMC = new CMultiAssetForwardSkewMC(m_hUnderlying->GetDateHandle());			
			pMC->Initialize(m_hUnderlying->m_assets, *this, mFixings, modelParameters, CMatrix(), CMatrix(), std::vector<std::vector<MlEqStrikeHandle> >());
			hMonteCarlo->generatePaths();
			hMonteCarlo->simulate(mResults, *this);
			
			// return the handle back
			begin_calculate(pVal, L"MonteCarloHandle")
				map_analytic_to_com(hMonteCarlo, MonteCarlo, spMonteCarlo);
				if (hr = pVal->AddValue(L"MonteCarloHandle", CComVariant(spMonteCarlo))) return hr;
			end_calculate
		} else {
			// non-basket case where a forward skew handle has NOT been passed in
			CMatrix										mFixings;
			long										nWidth = m_hFixingSchedule->GetWidth();			
			MlEqStochBetaVolatilityStructureHandle		hBetaVolStruct = dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*m_hUnderlying->GetVolatilityStructure());

			if (!hBetaVolStruct) throw "You need to supply a stochastic beta volatility matrix.";				
						
			// Initialise CappedFlooredCliquetPricer			
			if (!nWidth){
				// fixing schedule contains just dates - we get the values from the asset spot schedule
				m_hUnderlying->GetSpotsBeforeToday(m_hUnderlying->GetDateHandle()->GetDate(), afDates, &mFixings);
			} else if (nWidth == 1){
				// fixing schedule contains both the dates and the fixings			
				m_hFixingSchedule->GetColumnBeforeToday(m_hUnderlying->GetDateHandle()->GetDate(), 0, mFixings);
			} else {
				throw "The fixing schedule contains too many columns";
			}
			
			initialize(*m_hLocalCaps, *m_hLocalFloors, *m_hStrikes, *m_hCallPuts, m_fGlobalCap, m_fGlobalFloor, m_fGlobalRedemption, m_fGlobalGearing, m_fNotional, afDates, *m_hRainbowWeights );
		
			// Initialise Monte Carlo.			
			CForwardSkewMC* pMC;			
			hMonteCarlo = pMC = new CForwardSkewMC(hBetaVolStruct->getDateToDouble());						
			pMC->Initialize(*m_hUnderlying, *this, *hBetaVolStruct, mFixings, modelParameters, CMatrix(), CMatrix(), vector<vector<MlEqStrikeHandle> >());
			hMonteCarlo->generatePaths();
			
			// Return the forward skew handle in the result handle.
			begin_calculate(pVal, L"MonteCarloHandle")
				map_analytic_to_com(hMonteCarlo, MonteCarlo, spMonteCarlo);
				if (hr = pVal->AddValue(L"MonteCarloHandle", CComVariant(spMonteCarlo))) return hr;			
			end_calculate

			// simulate
			hMonteCarlo->simulate(mResults, *this);
		}
	} else {
		// too many data have been passed in.
		throw "You should pass in either an asset handle or a forward skew monte-carlo handle into this product - not both.";
	}
	
	begin_calculate(pVal, L"Price")
		if (hr = pVal->AddValue(L"Price", CComVariant(mResults[0][0]))) return hr;
	end_calculate

	
	end_function
}


STDMETHODIMP CBasketCliquet::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IARPropFitVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}