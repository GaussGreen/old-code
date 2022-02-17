// product_momentumcliquet.cpp : Implementation of CMomentumCliquet
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "products.h"
#include "product_momentumcliquet.h"
#include "MlEqObjects.h"
#include "MlEqParameterList.h"
#include "MonteCarlo.h"
#include "Mleqpde.h"

#undef max

HRESULT CMomentumCliquet::FinalConstruct(void)
{
	m_MinMaxType = Minimum;
	m_PayoffType = Straddle;
	return S_OK;
}

STDMETHODIMP CMomentumCliquet::Evaluate(IResult* pVal)
{			
	begin_function
				
	if (!m_hFixingSchedule) throw "No fixing schedule defined";
	if (!m_hStartsPeriodArr) throw "No start of period array defined";
	if (!m_hStrikesArr) throw "No strikes array defined";
	if (!m_hLocalFloorsArr) throw "No local floors array defined";
	if (!m_hLocalCapsArr) throw "No local caps array defined";
	if (!m_hWeightsArr) throw "No weights array defined";	
	if (!m_hCouponsArr) throw "No coupons array defined";
	if (!m_hGearingsArr) throw "No gearings array defined";	
	if (!m_hPeriodFloorsArr) throw "No period floors array defined";
	if (!m_hPeriodCapsArr) throw "No period caps array defined";

	if (!!m_hMonteCarlo){
		if (!m_hUnderlying){
			m_hUnderlying = m_hMonteCarlo->getAssets()[0];
		} else {
			throw "You passed an underlying AND a Monte Carlo handle";
		}
	}
	else
	{
		if (!m_hUnderlying)			throw "No underlying has been passed in";
		else if (!m_hParameterList) throw "No model info parameter list has been defined";
	}

		
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
	if (!!m_hPayDateSchedule)
	{
		if (m_hPayDateSchedule->GetWidth())
			throw "The pay date array should not contain any data other than the dates";

		anPayDates = m_hPayDateSchedule->GetDates();
	}
		
	// set up the spot fixings (mFixings)
	if (!nWidth)
	{		
		// fixing schedule contains just dates - we get the values from the asset spot schedule
		std::vector<long>			afDates;
		m_hFixingSchedule->GetDates(afDates);
		m_hUnderlying->GetSpotsBeforeToday(m_hUnderlying->GetDateHandle()->GetDate(), afDates, &mFixings);
	}
	else if (nWidth == 1)
	{
		// fixing schedule contains both the dates and the fixings		
		m_hFixingSchedule->GetColumnsBeforeToday(m_hUnderlying->GetDateHandle()->GetDate(), mFixings);
	}
	else
	{ 
		throw "The fixing schedule contains too many columns";
	}

	// useful for further changes in the GeneralCliquerPricer class...
	if( m_MinMaxType != Minimum && m_MinMaxType != Maximum )	throw "Wrong enum type";
	if( m_PayoffType != Straddle)	throw"You can price only straddle cliquets";
	
	
	initialize(m_hUnderlying, 
				anPayoffDates, 
				*m_hStartsPeriodArr, 
				*m_hStrikesArr, 
				*m_hLocalFloorsArr, 
				*m_hLocalCapsArr, 
				*m_hWeightsArr, 
				*m_hCouponsArr, 
				*m_hGearingsArr, 
				*m_hPeriodFloorsArr, 
				*m_hPeriodCapsArr, 
				anPayDates, 
				m_MinMaxType, 
				m_PayoffType,
				mFixings
				);


	if (!!m_hMonteCarlo)	// only for tests
	{
		hMonteCarlo = m_hMonteCarlo;
	}
	else
	{
		// set up the Monte Carlo
		m_hParameterList->BeginExhaust();
		m_hParameterList->GetValue("cL, cR, addTanhWing, yPower, seed, npaths, calibflag, numberVolStates, localVolFlag, savevolvolinfo, numberGridPoints, saveVolGridFlag, randomNumberFlag, ControlVariate, GlobalCalib, UseHermite = 0", ",", "=", &vModelParameters[0]);
		m_hParameterList->CheckUnexhausted();

		CForwardSkewMC* pMC;		
		hMonteCarlo = pMC = new CForwardSkewMC(hVolatility->getDateToDouble());
		pMC->Initialize(*m_hUnderlying, *this, *hVolatility, mFixings, vModelParameters, CMatrix(), CMatrix(), vector<vector<MlEqStrikeHandle> >());
		hMonteCarlo->generatePaths();
	}
	
	
	hMonteCarlo->simulate(mResults, *this);

	// get the price
	begin_calculate(pVal, L"Price")
		if (hr = pVal->AddValue(L"Price", CComVariant(mResults[0][0]))) return hr;
	end_calculate

	end_function
}





STDMETHODIMP CMomentumCliquet::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = {&IID_IHimalayan};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}