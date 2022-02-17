// product_callableCliquetSwap.cpp : Implementation of CCallableCliquetSwap


/////////////////////////////////////////////////////////////////////////////
// CCallableCliquetSwap



#include "stdafx.h"
#include "products.h"
#include "product_callablecliquetswap.h"
#include "MlEqObjects.h"
#include "MlEqParameterList.h"
#include "MonteCarlo.h"
//#include "MlEqSwap.h"
#include "mleqpde.h"

#undef max
#undef min


HRESULT CCallableCliquetSwap::FinalConstruct(void)
{	
//	m_fFunding = 
	m_fCallableLevel = 0.0;
	return S_OK;
}

STDMETHODIMP CCallableCliquetSwap::Evaluate(IResult* pVal)
{
	begin_function
			
	if (!m_hLocalCaps) throw "No local caps defined";
	if (!m_hLocalFloors) throw "No local floors defined";
	if (!m_hStrikes) throw "No strikes defined";
	if (!m_hCallPuts) throw "No call / put schedule defined";
	if (!m_hWeights) throw "No weights defined";
	if (!m_hFixingSchedule) throw "No fixing schedule defined";
	if (!m_hPeriodIdentifier) throw "No period identifiers defined";

	if (!m_hCouponSchedule) throw "No coupon dates schedule defined";
	if (!m_hGlobalCap)		throw "No global cap defined";
	if (!m_hGlobalFloor)	throw "No global floor defined";
	if (!m_hFixedCoupon)	throw "No fixed coupon defined";

//	if (!m_hZeroCurve)		throw "No zero curve defined...";
//	if (!m_hFundingSchedule) throw "No funding schedule defined";
//	if (!m_hSwap) throw "No swap component defined";
	

	HRESULT								hr;
	MlEqMonteCarloHandle				hMonteCarlo;
	CMatrix								mResults;

	GVector<long>	couponDates;
	GVector<long>	fundingDates;
	m_hCouponSchedule->GetDates(couponDates);
//	m_hFundingSchedule->GetDates(fundingDates);


	if (!m_hUnderlying && !m_hMonteCarlo)
		throw "Either a value for parameter 'ForwardSkewHandle' or 'Underlying' must be supplied";
	
	else if (!m_hUnderlying)
	{
		// m_hMonteCarlo passed in. The current date is set from MlEqDate (since it cannot be obtained from an asset).
		GVector<long>					afDates;
		CMatrix							mFixings;
		
		hMonteCarlo = m_hMonteCarlo;								
		if (m_hFixingSchedule->GetWidth() != 1) throw "Invalid fixing schedule. It should be a two column matrix of dates and values.";				
		m_hFixingSchedule->GetDates(afDates);
		m_hFixingSchedule->GetColumnBeforeToday(MlEqDate::GetCurrentDate() /*ToDo - ROMAIN this must go before it's finished.*/, 0, mFixings);	
		
		initialize( *m_hLocalCaps,
					*m_hLocalFloors,
					*m_hStrikes,
					*m_hCallPuts,
					*m_hWeights,
					*m_hGlobalCap,
					*m_hGlobalFloor,
					*m_hFixedCoupon,
					m_fCallableLevel,
					afDates,
					*m_hPeriodIdentifier,
					couponDates,
		//			m_hZeroCurve,
		//			fundingDates,
		//			m_fFunding,
//					m_hSwap,
					mFixings);
		
		
		hMonteCarlo->initializeMcHelper(*this, mFixings);
		hMonteCarlo->simulate(mResults, *this);
	} 
	else if (!m_hMonteCarlo)
	{
		// an asset handle has been passed in
		CVector										modelParameters(16);
		GVector<long>								afDates;
				
		if (!m_hParameterList) throw "No model info parameter list has been defined";
		m_hParameterList->BeginExhaust();
		m_hParameterList->GetValue("cL, cR, addTanhWing, yPower, seed, npaths, calibflag, numberVolStates, localVolFlag, savevolvolinfo, numberGridPoints, saveVolGridFlag, randomNumberFlag, ControlVariate, GlobalCalib, UseHermite = 0", ",", "=", &modelParameters[0]);
		m_hParameterList->CheckUnexhausted();
		m_hFixingSchedule->GetDates(afDates);
								
		CMatrix										mFixings;
		long										nWidth = m_hFixingSchedule->GetWidth();			
		MlEqStochBetaVolatilityStructureHandle		hBetaVolStruct = dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*m_hUnderlying->GetVolatilityStructure());

		if (!hBetaVolStruct) throw "You need to supply a stochastic beta volatility matrix.";				
					
				
		if (!nWidth)
			m_hUnderlying->GetSpotsBeforeToday(m_hUnderlying->GetDateHandle()->GetDate(), afDates, &mFixings);
		else if (nWidth == 1)	
			m_hFixingSchedule->GetColumnBeforeToday(m_hUnderlying->GetDateHandle()->GetDate(), 0, mFixings);
		else 
			throw "The fixing schedule contains too many columns";


		
		initialize( *m_hLocalCaps,
					*m_hLocalFloors,
					*m_hStrikes,
					*m_hCallPuts,
					*m_hWeights,
					*m_hGlobalCap,
					*m_hGlobalFloor,
					*m_hFixedCoupon,
					m_fCallableLevel,
					afDates,
					*m_hPeriodIdentifier,
					couponDates,			
					mFixings);
		
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
	} else {
		throw "You should pass in either an asset handle or a forward skew monte-carlo handle into this product - not both.";
	}
	

	begin_calculate(pVal, L"Price")
		if (hr = pVal->AddValue(L"Price", CComVariant(mResults[0][0]))) return hr;
	end_calculate

	end_function
}




STDMETHODIMP CCallableCliquetSwap::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_ICallableCliquetSwap
	};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}
