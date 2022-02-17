
#include "stdafx.h"
#include "Products.h"
#include "product_whale.h"
#include "MlEqObjects.h"
#include "MlEqParameterList.h"
#include "MonteCarlo.h"
#include "MlEqPde.h"



/////////////////////////////////////////////////////////////////////////////
// CWhale

////////////////////////////////////////////////////////////////////////////


HRESULT CWhale::FinalConstruct(void)
{
	m_fStrike = 0.0;
	m_cpe = Call;

	return S_OK;
}

STDMETHODIMP CWhale::Evaluate(IResult* pVal)
{			
	begin_function
				
	if (!m_hFixingSchedule) throw "No fixing schedule defined";
	if (!m_hWeights) throw "No weights array defined";

	if (!!m_hMonteCarlo)
	{
		if (!m_hUnderlying)
			m_hUnderlying = m_hMonteCarlo->getAssets()[0];
		else 
			throw"You passed an underlying AND a Monte Carlo handle";
	}
	else
	{
		if (!m_hUnderlying)			throw "No underlying has been passed in";
		else if (!m_hParameterList) throw "No model info parameter list has been defined";
	}

		
	HRESULT											hr;	
	MlEqMonteCarloHandle							hMonteCarlo;
	CMatrix											mFixings;			
	CVector											vModelParameters(16);
	CMatrix											mResults;
	MlEqStochBetaVolatilityStructureHandle			hVolatility = dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*m_hUnderlying->GetVolatilityStructure());
	long											nWidth = m_hFixingSchedule->GetWidth();
		
	if (!hVolatility) throw "You need to supply a stochastic beta volatility matrix";

	GVector< long > anPayoffDates;
	anPayoffDates = m_hFixingSchedule->GetDates();

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

	long nToday = m_hUnderlying->GetDateHandle()->GetDate();

	double cp;
	if(m_cpe == Call)		cp = -1.;
	else if(m_cpe == Put)	cp =  1.;
	else	throw"This payout is not supported";

	double	df	= m_hUnderlying->GetPayZeroCurve(true)->GetDiscountFactor(nToday, m_datePay);
	
	initialize(*m_hWeights, m_fStrike, cp, nToday, anPayoffDates, mFixings);
	

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

	mResults[0][0] *= df ;

	// get the price
	begin_calculate(pVal, L"Price")
		if (hr = pVal->AddValue(L"Price", CComVariant(mResults[0][0]))) return hr;
	end_calculate

	end_function
}





STDMETHODIMP CWhale::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = {&IID_IHimalayan};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}