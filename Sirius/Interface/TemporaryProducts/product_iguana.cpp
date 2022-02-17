// product_iguana.cpp : Implementation of CIguana
#include "stdafx.h"
#include "TemporaryProducts.h"
#include "MlEqObjects.h"
#include "product_iguana.h" 
#include "MonteCarlo.h"
#include "MlEqParameterList.h"


/*
STDMETHODIMP CIguana::Evaluate(IResult *pVal)
{
	begin_function

	HRESULT								hr;	
								
	if (!m_hUnderlying) throw "No underlying defined";	

					
	begin_calculate(pVal, L"Price")


		double price = BSFormula(m_fLowStrike,1.0);
		double midstrike = 0.5*(m_fHighStrike+m_fLowStrike);
		price -= 2*BSFormula( midstrike,1.0);
		price += BSFormula(m_fHighStrike,1.0);

		double deltaStrike = 0.5*(m_fHighStrike-m_fLowStrike);
		price /= deltaStrike*deltaStrike ;

//		double price = impliedVol(m_fLowStrike);
	
		if (hr = pVal->AddValue(L"Price", CComVariant(price))) return hr;
	end_calculate


	end_function
}
/*/


STDMETHODIMP CIguana::Evaluate(IResult *pVal)
{
	begin_function
			
	if (!m_hFixingSchedule) throw "No fixing schedule defined";

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


	if (!m_hUnderlying && !m_hMonteCarlo)
		throw "Either a value for parameter 'ForwardSkewHandle' or 'Underlying' must be supplied";

	else if (!m_hUnderlying)	// m_hMonteCarlo passed in. The current date is set from MlEqDate (since it cannot be obtained from an asset).
	{
		
		hMonteCarlo = m_hMonteCarlo;								
		if (m_hFixingSchedule->GetWidth() != 1)
			throw "Invalid fixing schedule. It should be a two column matrix of dates and values.";	
		
		GVector<long>		afDates;
		CMatrix				mFixings;
		m_hFixingSchedule->GetDates(afDates);
		m_hFixingSchedule->GetColumnBeforeToday(MlEqDate(MlEqDate::GetCurrentDate()).GetDate(), 0, mFixings);	
		
		initialize(m_fHighStrike, m_fLowStrike, afDates, mFixings[0]);//////////
	
		hMonteCarlo->initializeMcHelper(*this, mFixings);
		hMonteCarlo->simulate(mResults, *this);
	}
	else if (!m_hMonteCarlo)	// an asset handle has been passed in
	{		
		
		MlEqStochBetaVolatilityStructureHandle		
			hBetaVolStruct = dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*m_hUnderlying->GetVolatilityStructure());

		if (!hBetaVolStruct) throw "You need to supply a stochastic beta volatility matrix.";
		
		GVector<long>	afDates;				
		m_hFixingSchedule->GetDates(afDates);
	
		CMatrix			mFixings;
		long			nWidth = m_hFixingSchedule->GetWidth();	
			
		if (!nWidth)			// fixing schedule contains just dates - we get the values from the asset spot schedule
			m_hUnderlying->GetSpotsBeforeToday(m_hUnderlying->GetDateHandle()->GetDate(), afDates, &mFixings);
			
		else if (nWidth == 1)	// fixing schedule contains both the dates and the fixings	
			m_hFixingSchedule->GetColumnBeforeToday(m_hUnderlying->GetDateHandle()->GetDate(), 0, mFixings);
		
		else 
			throw "The fixing schedule contains too many columns";
		
		initialize(m_fHighStrike, m_fLowStrike, afDates, mFixings[0]);
		

		CVector		modelParameters(16);
				
		if (!m_hParameterList) 
			throw "No model info parameter list has been defined";

		m_hParameterList->BeginExhaust();
		m_hParameterList->GetValue("cL, cR, addTanhWing, yPower, seed, npaths, calibflag, numberVolStates, localVolFlag, savevolvolinfo, numberGridPoints, saveVolGridFlag, randomNumberFlag, ControlVariate, GlobalCalib, UseHermite = 0", ",", "=", &modelParameters[0]);
		m_hParameterList->CheckUnexhausted();
	


		// Initialise Monte Carlo.			
		CForwardSkewMC* pMC;		
		hMonteCarlo = pMC = new CForwardSkewMC(hBetaVolStruct->getDateToDouble());						
		pMC->Initialize(*m_hUnderlying, *this, *hBetaVolStruct, mFixings, modelParameters, CMatrix(), CMatrix(), std::vector<std::vector<MlEqStrikeHandle> >());
		hMonteCarlo->generatePaths();		
		
		// Return the forward skew handle in the result handle.
		begin_calculate(pVal, L"MonteCarloHandle")
			map_analytic_to_com(hMonteCarlo, MonteCarlo, spMonteCarlo);
			if (hr = pVal->AddValue(L"MonteCarloHandle", CComVariant(spMonteCarlo))) return hr;			
		end_calculate

		// simulate
		hMonteCarlo->simulate(mResults, *this);		
	}
	 else 
		throw "You should pass in either an asset handle or a forward skew monte-carlo handle into this product - not both.";
	
	
	begin_calculate(pVal, L"Price")
		if (hr = pVal->AddValue(L"Price", CComVariant(mResults[0][0]))) return hr;
	end_calculate

	end_function
}



HRESULT CIguana::FinalConstruct()
{
	m_fLowStrike = m_fHighStrike = 0.0;
	m_dateMaturity = MlEqDate(MlEqDate::GetCurrentDate()).AddTenor("1y")->GetDate();
	return S_OK;
}


STDMETHODIMP CIguana::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IIguana
	};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}


double CIguana::BSFormula(double strike, double cp)
{
	long nToday = m_hUnderlying->GetDateHandle()->GetDate();
	double fForward = m_hUnderlying->GetQuantoForward(nToday, m_dateMaturity, false);
	double fVolatility = m_hUnderlying->GetCompositeVolatility(MlEqStrike(strike), nToday, m_dateMaturity);	
	double fTime = m_hUnderlying->GetVolatilityStructure()->getDateToDouble()->GetYearFraction(m_dateMaturity);				
	double d1 = (log(fForward / strike) + 0.5 * fVolatility * fVolatility * fTime) / (fVolatility * sqrt(fTime));	
	double d2 = (log(fForward / strike) - 0.5 * fVolatility * fVolatility * fTime) / (fVolatility * sqrt(fTime));
	double fDiscountFactor = m_hUnderlying->GetPayZeroCurve(true)->GetDiscountFactor(nToday, m_dateMaturity);

	return fDiscountFactor * ( cp*fForward * ::normal(cp*d1) - cp*strike * ::normal(cp*d2) );
}


double CIguana::impliedVol(double strike)
{
	long nToday = m_hUnderlying->GetDateHandle()->GetDate();
	
	return m_hUnderlying->GetCompositeVolatility(MlEqStrike(strike), nToday, m_dateMaturity);	
}




/////////////////////////////////


