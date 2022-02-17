//	product_himalayan.cpp : Implementation of CHimalayan
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Products.h"
#include "product_himalayan.h"
#include "mleqobjects.h"
#include "MlEqParameterList.h"
#include "MonteCarlo.h"
#include "mleqpde.h"

STDMETHODIMP CHimalayan::Evaluate(IResult* pVal)
{
	begin_function

	if (!m_hFixingSchedule) throw "No fixing schedule defined";	
	if (!m_hIsHimalayaArr) throw "No Himalaya schedule defined";
	if (!m_hRainbowWeights) throw "No rainbow weights defined";	
	if (!m_hCallLevelsArr) throw "No call levels defined";
	if (!m_hRebates) throw "No rebates defined";
		
	HRESULT											hr;	
	double											fCallPut = (m_CallOrPut == Call ? 1.0 : -1.0);	
	MlEqMonteCarloHandle							hMonteCarlo;
	CMatrix											mFixings;			
	CVector											vModelParameters(16);
	CMatrix											results;
	CVector											vSpotFixings;
	long											nWidth = m_hFixingSchedule->GetWidth();
	
	// set up vSpotFixings from m_hAssets, m_dateFixing and m_hFixingValues
	if (!m_hFixingValues){
		// derive the fixings from the asset handle	
		vSpotFixings.resize(m_hAssets.size());
		for (long nAsset = 0; nAsset < m_hAssets.size(); nAsset++){
			vSpotFixings[(int)nAsset] = m_hAssets[nAsset]->GetSpot(m_dateFixing);
		}
	} else {
		if (m_hFixingValues->getsize() != m_hAssets.size()) throw "The number of elements in the fixing values array must equal the number of assets";
		vSpotFixings = *m_hFixingValues;
	}

	// initialise the Himalaya
	initialize(m_hAssets.size(),
			   vSpotFixings,
			   m_hFixingSchedule->GetDates(),
			   *m_hIsHimalayaArr,
			   *m_hRainbowWeights,
			   m_fNotional, 
			   m_fStrike, 
			   fCallPut,			   
			   *m_hCallLevelsArr, 
			   *m_hRebates,
			   m_bDelayRebateToEnd,
			   m_bAsianFromStart);


	// set up the fixings
	if (!nWidth){		
		// fixing schedule contains just dates - we get the values from the assets spot schedules		
		std::vector<long>			afDates;
		m_hFixingSchedule->GetDates(afDates);
		mFixings.resize(m_hAssets.size(), m_hFixingSchedule->size(), 0.0);
		for (long nAsset = 0; nAsset < m_hAssets.size(); nAsset++){			
			std::vector<double>	vSpots;
			m_hAssets[nAsset]->GetSpotsBeforeToday(m_hAssets[nAsset]->GetDateHandle()->GetDate(), afDates, &vSpots);
			for (long nSpot = 0; nSpot < vSpots.size(); nSpot++){
				mFixings[(int)nAsset][(int)nSpot] = vSpots[nSpot];
			}
		}
	} else if (nWidth == m_hAssets.size()){
		// fixing schedule contains both the dates and the fixings		
		m_hFixingSchedule->GetColumnsBeforeToday(m_hAssets[0]->GetDateHandle()->GetDate(), mFixings);
	} else {
		throw "The number of columns in the fixing schedule is invalid";
	}
	
	// set up the Monte Carlo
	CMultiAssetForwardSkewMC* pMC;	
	hMonteCarlo = pMC = new CMultiAssetForwardSkewMC(m_hAssets[0]->GetDateHandle());
	
	m_hParameterList->BeginExhaust();
	m_hParameterList->GetValue("cL, cR, addTanhWing, yPower, seed, npaths, calibflag, numberVolStates, localVolFlag, savevolvolinfo, numberGridPoints, saveVolGridFlag, randomNumberFlag, ControlVariate, GlobalCalib, UseHermite", std::string(","), &vModelParameters[0]);
	m_hParameterList->CheckUnexhausted();
	pMC->Initialize(m_hAssets, *this, mFixings, vModelParameters, CMatrix(), CMatrix(), std::vector<std::vector<MlEqStrikeHandle> >());
	hMonteCarlo->generatePaths();
		
	// return the handle back
	begin_calculate(pVal, L"MonteCarloHandle")
		map_analytic_to_com(hMonteCarlo, MonteCarlo, spMonteCarlo);
		if (hr = pVal->AddValue(L"MonteCarloHandle", CComVariant(spMonteCarlo))) return hr;
	end_calculate
	
	// get the price
	hMonteCarlo->simulate(results, *this);
	begin_calculate(pVal, L"Price")
		if (hr = pVal->AddValue(L"Price", CComVariant(results[0][0]))) return hr;
	end_calculate

	end_function
}

HRESULT CHimalayan::FinalConstruct(void)
{	
	m_fStrike = 1.0;	
	m_fNotional = 1.0;
	m_bAsianFromStart = false;
	m_bDelayRebateToEnd = true;
	m_CallOrPut = Call;
	m_dateFixing = 0.0;
	m_hIsHimalayaArr = new MlEqArray;
	m_hRainbowWeights = new MlEqArray;
	m_hCallLevelsArr = new MlEqArray; 
	return S_OK;
}

STDMETHODIMP CHimalayan::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = {&IID_IHimalayan};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}