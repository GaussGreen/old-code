//	product_barriercliquet.cpp : Implementation of CBarrierCliquet
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "products.h"
#include "product_barriercliquet.h"
#include "MlEqObjects.h"
#include "MlEqParameterList.h"
#include "MonteCarlo.h"
#include "mleqpde.h"

HRESULT CBarrierCliquet::FinalConstruct(void)
{	
	m_bDelayPaymentsToEnd = false;
//	m_datePay = 0.0;
	return S_OK;
}


STDMETHODIMP CBarrierCliquet::Evaluate(IResult* pVal)
{
	begin_function

	if(!m_bDelayPaymentsToEnd) throw "Delay Payments To End cannot be set to NO, because using black scholes approximation in CBarrierCliquet::Evaluate";
	if (!m_hBarrierStartDates) throw "No StartDates schedule defined CBarrierCliquet::Evaluate";
	if (!m_hBarrierEndDates) throw "No EndDates barrier schedule defined CBarrierCliquet::Evaluate";
	if (!m_hBarrierPayDates) throw "No PayDates schedule defined CBarrierCliquet::Evaluate";
	if (!m_hUnderlying) throw "No underlying defined";
	if (!m_hBarriers) throw "No barriers defined";
	if (!m_hRebate) throw "No coupons defined";
	if (!m_hCallsPuts) throw "No callsputs defined";
	if (!m_hStrikes) throw "No strikes defined";
	if (!m_hBarrierHasHit) throw "No events defined";
	if (!m_hBarrierFixings) throw "No barrier fixings defined";

//	if (!m_bKnockIn) throw "No KnockIn defined";


	HRESULT										hr;

	MlEqMonteCarloHandle						hMonteCarlo;			
	MlEqStochBetaVolatilityStructureHandle		hVolatility = dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*m_hUnderlying->GetVolatilityStructure());	// ToDo - this will break for barrier cliquet on basket case
	CMatrix										mResults;
	CVector										vModelParameters(16);
	
	double										fVolatilityStrikePlace;
	double										fBlend;
	long										nPricingMethod;
	double										stressForwardBarrierVol;
	double										stressForwardBarrierVolSlope;
	long										nPoints;

	m_hParameterList->BeginExhaust();
	m_hParameterList->GetValue("cL, cR, addTanhWing, yPower, seed, npaths, calibflag, numberVolStates, localVolFlag, savevolvolinfo, numberGridPoints, saveVolGridFlag, randomNumberFlag, ControlVariate, GlobalCalib, UseHermite = 0", ",", "=", &vModelParameters[0]);
	m_hParameterList->GetValue("VolatilityStrikePlace", &fVolatilityStrikePlace);	
	m_hParameterList->GetValue("Blend", &fBlend);
	m_hParameterList->GetValue("nPoints", &nPoints);
	m_hParameterList->GetValue("StressForwardBarrierVol", &stressForwardBarrierVol);
	m_hParameterList->GetValue("StressForwardBarrierVolSlope", &stressForwardBarrierVolSlope);
	m_hParameterList->GetValue("BarrierPricingMethod", &nPricingMethod);
	m_hParameterList->CheckUnexhausted();

	GMatrix < bool > barrierHasHit(m_hBarrierHasHit->rows(),m_hBarrierHasHit->cols());
	for ( int i = 0 ; i < m_hBarrierHasHit->rows(); i++ ){
		for ( int j = 0 ; j < m_hBarrierHasHit->cols(); j++ ){
			barrierHasHit[i][j] = (*m_hBarrierHasHit)[i][j] ? true : false;
		}
	}

	GVector< long> BarrierStartDates;
	BarrierStartDates = m_hBarrierStartDates->GetDates();

	GVector< long> BarrierEndDates;
	BarrierEndDates = m_hBarrierEndDates->GetDates();

	GVector< long> BarrierPayDates(m_hBarrierPayDates->getsize());

	for ( int i = 0 ; i < m_hBarrierPayDates->getsize(); i++ ){
		BarrierPayDates[i] = (*m_hBarrierPayDates)[i];
	}


	int nToday = m_hUnderlying->GetDateHandle()->GetDate();
	double spot = m_hUnderlying->GetSpot(nToday);

	MlEqZeroCurveHandle hZCurve = m_hUnderlying->GetPayZeroCurve(true);

	initialize(		hZCurve,
					*m_hBarriers,
					*m_hRebate,
					barrierHasHit,
					BarrierStartDates,
					BarrierEndDates,
					BarrierPayDates,
					nPricingMethod,
					stressForwardBarrierVol,
					stressForwardBarrierVolSlope,
					fVolatilityStrikePlace,
					fBlend,
					nPoints,
					m_hUnderlying->GetDateHandle(),
					*m_hBarrierFixings,
					spot,
					m_bKnockIn
					); 


	// set up the Monte Carlo
	CForwardSkewMC* pMC;
	hMonteCarlo = pMC = new CForwardSkewMC(hVolatility->getDateToDouble());

	CMatrix fixings;
	CMatrix calibTimeIndex,calibVols;
	vector<vector<MlEqStrikeHandle> > pCalibStrikes;

	const std::vector<long> additionalSimDates;	
	pMC->Initialize(*m_hUnderlying, *this, additionalSimDates, *hVolatility, fixings, vModelParameters, calibTimeIndex,calibVols,pCalibStrikes);
	hMonteCarlo->generatePaths();


	// get the price	
	hMonteCarlo->simulate(mResults, *this);
	begin_calculate(pVal, L"Price")
		if (hr = pVal->AddValue(L"Price", CComVariant(mResults[0][0]))) return hr;
	end_calculate

	end_function

}

STDMETHODIMP CBarrierCliquet::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IARPropFitVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}