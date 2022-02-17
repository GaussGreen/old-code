//	product_napoleon.cpp : Implementation of CNapoleon
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "products.h"
#include "product_napoleon.h"
#include "MlEqObjects.h"
#include "MlEqParameterList.h"
#include "MonteCarlo.h"
#include "MLeqpde.h"


HRESULT CNapoleon::FinalConstruct(void)
{
	m_bDelayPaymentsToEnd = false;
	return S_OK;
}

STDMETHODIMP CNapoleon::Evaluate(IResult* pVal)
{			
	begin_function
				
	if (!m_hFixingSchedule) throw "No fixing schedule defined";
	if (!m_hStartsPeriodArr) throw "No start of period array defined";
	if (!m_hCallPutArr) throw "No call / put array defined";
	if (!m_hStrikesArr) throw "No strikes array defined";
	if (!m_hLocalFloorsArr) throw "No local floors array defined";
	if (!m_hLocalCapsArr) throw "No local caps array defined";
	if (!m_hWeightsArr) throw "No weights array defined";	
	if (!m_hCouponsArr) throw "No coupons array defined";
	if (!m_hGearingsArr) throw "No gearings array defined";	
	if (!m_hPeriodFloorsArr) throw "No period floors array defined";
	if (!m_hPeriodCapsArr) throw "No period caps array defined";
	if (!m_hParameterList) throw "No model info parameter list has been defined";
	if (!m_hUnderlying) throw "No underlying defined";
		
	HRESULT											hr;		
	GVector<long>									anPayoffDates;		
	GVector<long>									anPayDates;
	MlEqMonteCarloHandle							hMonteCarlo;
	CMatrix											mFixings;			
	CVector											vModelParameters(16);
	CMatrix											mResults;
	MlEqStochBetaVolatilityStructureHandle			hVolatility = dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*m_hUnderlying->GetVolatilityStructure());
	long											nWidth = m_hFixingSchedule->GetWidth();
		
	if (!hVolatility) throw "You need to supply a stochastic beta volatility matrix";

	// initialise the MlEqHimalaya
	anPayoffDates = m_hFixingSchedule->GetDates();
	if (!!m_hPayDateScheduleOpt){
		if (m_hPayDateScheduleOpt->GetWidth()) throw "The pay date array should not contain any data other than the dates";
		anPayDates = m_hPayDateScheduleOpt->GetDates();
	}
	initialize(m_hUnderlying, anPayoffDates, *m_hStartsPeriodArr, *m_hCallPutArr, *m_hStrikesArr, *m_hLocalFloorsArr, *m_hLocalCapsArr, *m_hWeightsArr, *m_hCouponsArr, *m_hGearingsArr, *m_hPeriodFloorsArr, *m_hPeriodCapsArr, anPayDates, m_bDelayPaymentsToEnd);
		
	// set up the spot fixings (mFixings)
	if (!nWidth){		
		// fixing schedule contains just dates - we get the values from the asset spot schedule
		std::vector<long>			afDates;
		m_hFixingSchedule->GetDates(afDates);
		m_hUnderlying->GetSpotsBeforeToday(m_hUnderlying->GetDateHandle()->GetDate(), afDates, &mFixings);
	} else if (nWidth == 1){
		// fixing schedule contains both the dates and the fixings		
		m_hFixingSchedule->GetColumnBeforeToday(m_hUnderlying->GetDateHandle()->GetDate(), 0, mFixings);
	} else {
		throw "The fixing schedule contains too many columns";
	}
			
	// set up the Monte Carlo
	m_hParameterList->BeginExhaust();
	m_hParameterList->GetValue("cL, cR, addTanhWing, yPower, seed, npaths, calibflag, numberVolStates, localVolFlag, savevolvolinfo, numberGridPoints, saveVolGridFlag, randomNumberFlag, ControlVariate, GlobalCalib, UseHermite = 0", ",", "=", &vModelParameters[0]);
	m_hParameterList->CheckUnexhausted();

	CForwardSkewMC* pMC;	
	hMonteCarlo = pMC = new CForwardSkewMC(hVolatility->getDateToDouble());
	pMC->Initialize(*m_hUnderlying, *this, *hVolatility, mFixings, vModelParameters, CMatrix(), CMatrix(), vector<vector<MlEqStrikeHandle> >());
	hMonteCarlo->generatePaths();

	// get the price
	hMonteCarlo->simulate(mResults, *this);
	begin_calculate(pVal, L"Price")
		if (hr = pVal->AddValue(L"Price", CComVariant(mResults[0][0]))) return hr;
	end_calculate

	end_function
}

STDMETHODIMP CNapoleon::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = {&IID_IHimalayan};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}