//	product_europeanoption.cpp : Implementation of CEuropeanOption
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "products.h"
#include "product_europeanoption.h"
#include "MlEqParameterList.h"
#include "MonteCarlo.h"
#include "mleqpde.h"


STDMETHODIMP CEuropeanOption::Evaluate(IResult* pVal)
{
	begin_function	
	if (!m_hUnderlying) throw "No underlying defined";

	HRESULT								hr;
	long								nToday = m_hUnderlying->GetDateHandle()->GetDate();
	double								fStrike(m_fStrike);
	double								fNotional(m_fNotional);
	NotionalTypeEnum					nte(m_nte);
	MlEqDateHandle						hDate = new MlEqDate(*m_hUnderlying->GetDateHandle());	
	double								fCallPut = m_cpe == Call ? 1.0 : -1.0;
	double								fDiscountFactor	= m_hUnderlying->GetPayZeroCurve(true)->GetDiscountFactor(nToday, m_datePay);
	
	if ( fabs(m_dateStart) < 1e-6 ){
		m_dateStart = nToday;
	}

	if (nToday >= m_datePay){
		// We discount nToday == m_datePay since we only really care about precise valuation in the afternoon after which
		// back office have closed out the posisiton.
		begin_calculate(pVal, L"Price")
			if (hr = pVal->AddValue(L"Price", CComVariant(0.0))) return hr;
		end_calculate
	} else {
		ConvertConstantNotionalRamStrike(&nte, m_dateStart, nToday, &fNotional, &fStrike);
		if (m_dateStart > nToday) hDate->PutDate(m_dateStart);
		
		// Use Monte-Carlo if the option is forward starting and m_hParameterList has been passed.
		if ( m_hUnderlying->IsBasket()){
			hr = EvaluateMonteCarlo(&pVal);
		}  else if(m_dateStart > nToday && !!m_hParameterList){
			hr = EvaluateMonteCarlo(nToday, hDate, fDiscountFactor, nte, fCallPut, fStrike, fNotional, &pVal);
		} else {
			hr = EvaluateBlackScholes(nToday, hDate, fDiscountFactor, nte, fCallPut, fStrike, fNotional, &pVal);
		}
	}
	return hr;
	end_function
}

HRESULT CEuropeanOption::EvaluateBlackScholes(long nToday, MlEqDateHandle asOfDate, double fDiscountFactor, NotionalTypeEnum nte, double fCallPut, double fStrike, double fNotional, IResult** ppVal)
{	
	HRESULT								hr;
	double								fPrice;
	

	if (nToday > m_datePay)
	{
		// option has expired and has been paid
		throw "The option has already been paid";
	}
	else if (nToday >= m_dateMaturity)
	{
		// option has expired but not yet paid
		double		fSpot = m_hUnderlying->GetSpot(m_dateMaturity);
					
		fPrice = MlEqMaths::Max(fCallPut * (fSpot - fStrike), 0.0);
		fPrice *= m_hUnderlying->GetPayZeroCurve(true)->GetDiscountFactor(nToday, m_datePay) * fNotional;
		
		begin_calculate(*ppVal, L"Price")
			if (hr = (*ppVal)->AddValue(L"Price", CComVariant(fPrice))) return hr;
		end_calculate
	} 
	else if (nToday >= m_dateStart)
	{
		// option has started but has not yet expired
		MlEqVolatilityStructureHandle		hVolatilityStructure = m_hUnderlying->GetVolatilityStructure();
		
		double		fSpot = m_hUnderlying->GetSpot(nToday);			
		double		fwd	= m_hUnderlying->GetQuantoForward(nToday, m_dateMaturity, false);
		double		fTime = asOfDate->GetYearFraction(m_dateMaturity);			
		double		volS;

		double fFixing;
		if (nte == ConstantNotional || nte == RAMConstantNotional)
		{
			fFixing = m_hUnderlying->GetSpot(m_dateStart);	
			fStrike *= fFixing;
			fNotional /= fFixing;
		}

		volS = hVolatilityStructure->MlEqVolatilityStructure::getFutureVol(MlEqStrike(fStrike), asOfDate->GetDate(), m_dateMaturity, fSpot);

		double		fVolatility = m_hUnderlying->GetCompositeVolatility(nToday, m_dateMaturity, volS);

		fPrice = ::Bs(fwd, fVolatility, fTime, fStrike, fDiscountFactor, fCallPut);
		fPrice *= fNotional;
	
		begin_calculate(*ppVal, L"Price")
			if (hr = (*ppVal)->AddValue(L"Price", CComVariant(fPrice))) return hr;	
		end_calculate	
	} 
	else
	{
		// option has not yet started		
		MlEqVolatilityStructureHandle		hVolatilityStructure = m_hUnderlying->GetVolatilityStructure();
			
		double		fForward = m_hUnderlying->GetQuantoForward(asOfDate->GetDate(), m_dateMaturity, false);
		double		fwdLong = m_hUnderlying->GetQuantoForward(nToday, m_dateMaturity, false);
		double		fwdShort = m_hUnderlying->GetQuantoForward(nToday, m_dateStart, false);
		double		fTime = asOfDate->GetYearFraction(m_dateMaturity);		
		double		volS = hVolatilityStructure->MlEqVolatilityStructure::getFutureVol(MlEqStrike(fStrike * fwdShort), asOfDate->GetDate(), m_dateMaturity, fwdShort);
		double		fVolatility = m_hUnderlying->GetCompositeVolatility(asOfDate->GetDate(), m_dateMaturity, volS);
		
		if (nte == ConstantNotional || nte == RAMConstantNotional){
			fPrice = ::Bs(fwdLong / fwdShort, fVolatility, fTime, fStrike, fDiscountFactor, fCallPut);
		} else { 
			fPrice = ::Bs(fwdLong, fVolatility, fTime, fStrike * fwdShort, fDiscountFactor, fCallPut);
		}

		fPrice *= fNotional;
		begin_calculate(*ppVal, L"Price")
			if (hr = (*ppVal)->AddValue(L"Price", CComVariant(fPrice))) return hr;
		end_calculate
	}

	return S_OK;
}

HRESULT CEuropeanOption::EvaluateMonteCarlo(long nToday, MlEqDateHandle hDate, double fDiscountFactor, NotionalTypeEnum nte, double fCallPut, double fStrike, double fNotional, IResult** ppVal)
{
	HRESULT								hr;		
	double								fPrice;		
	CVector								vModelParameters(16);
	MlEqMonteCarloHandle				hMonteCarlo;
	std::vector<long>					anPayoffDates;
			
	m_hParameterList->BeginExhaust();
	m_hParameterList->GetValue("cL = 0.01, cR = 3.0, addTanhWing = 1.0, yPower, seed, npaths, calibflag, numberVolStates, localVolFlag, savevolvolinfo, numberGridPoints, saveVolGridFlag, randomNumberFlag, ControlVariate, GlobalCalib, UseHermite = 0", ",", "=", &vModelParameters[0]);
	m_hParameterList->CheckUnexhausted();

	if (m_dateStart < nToday){		// not in use as bs formula --> this case actually crashes		
		anPayoffDates.resize(2);
		anPayoffDates[0] = m_dateStart;
		anPayoffDates[1] = m_dateMaturity;
	} else /*if(m_dateStart>=nToday)*/{
		int interval = m_dateMaturity - m_dateStart;
		int NumberOfDates = floor((m_dateMaturity - nToday) / interval) + 1;//+2; // +1

		if ((double)((m_dateMaturity - nToday) / interval) - floor((m_dateMaturity - nToday) / interval) > 0.0 ){
			NumberOfDates++;
		}

		anPayoffDates.resize(NumberOfDates);//+1?
		anPayoffDates[NumberOfDates - 1] = m_dateMaturity;
		anPayoffDates[NumberOfDates - 2] = m_dateStart;
		int nDate=NumberOfDates-3;				
		int	currentDate(m_dateStart);

		while (nDate > 0){// >=0?	
			currentDate -= interval;
			anPayoffDates[nDate--] = currentDate;
		}

		anPayoffDates[0]=nToday;
	}
	
	// set up the Monte Carlo
	CForwardSkewMC* pMC;	
	hMonteCarlo = pMC = new CForwardSkewMC(m_hUnderlying->GetDateHandle());
	CMatrix calibTimeIndex, calibVols;
	vector<vector<MlEqStrikeHandle> > pCalibStrikes;
	pMC->Initialize(*m_hUnderlying, anPayoffDates, vModelParameters, calibTimeIndex, calibVols, pCalibStrikes);
	hMonteCarlo->generatePaths();

	// compute the payoff
	CMatrix& firstPath = hMonteCarlo->GetPathArray(0);
	int nDate = firstPath.cols() - 1;
	if (nDate == 0){
		throw("incorret set up in cliquet montecarlo");
	}

	double fPayoff = 0.0;
	for (int nPath = 0; nPath < hMonteCarlo->m_nPaths; nPath++){
		double ST = hMonteCarlo->GetPathValue(nPath, nDate);
		double St = hMonteCarlo->GetPathValue(nPath, nDate - 1);

		if (nte == ConstantNotional || nte == RAMConstantNotional){
			fPayoff += MlEqMaths::Max(fCallPut * (ST - fStrike * St) / St, 0.0);
		} else {
			fPayoff += MlEqMaths::Max(fCallPut * (ST - fStrike * St), 0.0);
		}
	}

	hMonteCarlo->averageMonteCarloResults(fPayoff);
	fPrice = fPayoff * fNotional * fDiscountFactor;	
	begin_calculate(*ppVal, L"Price")
		if (hr = (*ppVal)->AddValue(L"Price", CComVariant(fPrice))) return hr;			
	end_calculate	
	return S_OK;
}



HRESULT CEuropeanOption::FinalConstruct(void)
{	
	m_fStrike = 0.0;
	m_dateStart = 0.0;
	m_dateMaturity = 0.0;
	m_fNotional = 1.0;
	m_nte = ConstantShare;
	m_cpe = Call;
	return S_OK;
}

STDMETHODIMP CEuropeanOption::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IARPropFitVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}



HRESULT CEuropeanOption::EvaluateMonteCarlo(IResult** ppVal)
{
	HRESULT								hr;		
	CVector								vModelParameters(16);
	MlEqMonteCarloHandle				hMonteCarlo;
	CMatrix								mResults;
	CMatrix								fixings;
			
	if (!m_hParameterList) throw "You have not supplied a model information parameter list";
	m_hParameterList->BeginExhaust();
	m_hParameterList->GetValue("cL, cR, addTanhWing, yPower, seed, npaths, calibflag, numberVolStates, localVolFlag, savevolvolinfo, numberGridPoints, saveVolGridFlag, randomNumberFlag, ControlVariate, GlobalCalib, UseHermite = 0", ",", "=", &vModelParameters[0]);
	m_hParameterList->CheckUnexhausted();

	double	fCallPut = m_cpe == Call ? 1.0 : -1.0;
	long nToday = m_hUnderlying->GetDateHandle()->GetDate();

	ConvertConstantNotionalRamStrike(&m_nte, m_dateStart, nToday, &m_fNotional, &m_fStrike);	// ugly but necessary for now (as of 05/11/04).

	initialize(m_dateStart, m_dateMaturity, m_fStrike, fCallPut, m_nte, m_hUnderlying, fixings);		   

	if (m_hUnderlying->IsBasket()){		
		CMultiAssetForwardSkewMC* pMC;
		long nDate1 = m_hUnderlying->m_assets[0]->GetDateHandle()->GetDate();
		long nDate2 = m_hUnderlying->m_assets[0]->GetCorrelationMatrix()->GetDate();
		ATLASSERT(nDate1 == nDate2);		
		hMonteCarlo = pMC = new CMultiAssetForwardSkewMC(m_hUnderlying->GetDateHandle());
		pMC->Initialize(m_hUnderlying->m_assets, *this, fixings, vModelParameters, CMatrix(), CMatrix(), std::vector<std::vector<MlEqStrikeHandle> >());
		hMonteCarlo->generatePaths();
		hMonteCarlo->simulate(mResults, *this);
	} else {
		throw "CEuropeanOption::EvaluateMonteCarlo should not be called for a single underlying asset";
	}

	double fDiscountFactor	= m_hUnderlying->GetPayZeroCurve(true)->GetDiscountFactor(nToday, m_datePay);
	double fPrice = mResults[0][0] * m_fNotional * fDiscountFactor;
	
	begin_calculate(*ppVal, L"Price")
		if (hr = (*ppVal)->AddValue(L"Price", CComVariant(fPrice))) return hr;
	end_calculate
	return S_OK;
}



void EuropeanBasketProduct::initialize(DATE dateStart, 
									   DATE dateMaturity, 
									   double fStrike, 
									   double cp, 
									   NotionalTypeEnum nte,
									   MlEqAssetHandle udly, 
									   CMatrix& fixings)
{
	m_fStrike = fStrike ;
	m_cp = cp;

	int nAsset = udly->m_assets.size();
	m_weight.resize(nAsset, 0.0);
	double sumWeight = 0.0;
	for(int i=0; i<nAsset; i++)
	{
		m_weight[i] = udly->m_assetWeights[i];
		sumWeight += m_weight[i];
	}
/*
	if( nte == ConstantNotional )
	{
		for(int i=0; i<nAsset; i++)
			m_weight[i] /= sumWeight;
	}
	else if( nte == ConstantShare )
	{
		if( dateStart >= udly->GetDateHandle()->GetDate() )
			m_fStrike *= sumWeight;
	}
	else
		throw " wrong notional type ";

*/
/*	if( nte == ConstantShare )
	{
		m_fStrike *= sumWeight;
	}
	else
		throw " wrong notional type ";
*/

	// build time subdivision in the forward starting case...

	long nToday = udly->GetDateHandle()->GetDate();
			

	if (dateStart < nToday)
	{			
		m_payoffDates.resize(2);
		m_payoffDates[0] = dateStart;
		m_payoffDates[1] = dateMaturity;
		m_nMaturity = 1;
		int pastSize = (dateMaturity < nToday)?2:1;

		// get past spots in a matrix...
		fixings.resize(nAsset, pastSize, 0.0);				
		for (int n = 0; n<nAsset ; n++)
		{
			CMatrix vSpots;
			udly->m_assets[n]->GetSpotsBeforeToday(nToday, m_payoffDates, &vSpots);
			for(int j=0; j<vSpots.cols();j++)	// should be pastSize
				fixings[n][j] = vSpots[0][j];			
		}
	}
	else 
	{
		int interval = dateMaturity - dateStart;
		int NumberOfDates = floor((dateMaturity - nToday) / double(interval)) + 1;

		if ((double)((dateMaturity - nToday) / interval) - floor((dateMaturity - nToday) / interval) > 0.0 ){
			NumberOfDates++;
		}

		m_payoffDates.resize(NumberOfDates);//+1?
		m_payoffDates[NumberOfDates - 1] = dateMaturity;
		m_payoffDates[NumberOfDates - 2] = dateStart;
		
		m_nMaturity = NumberOfDates - 1;

		int nDate = NumberOfDates-3;				
		int	currentDate(dateStart);

		while (nDate > 0)
		{	
			currentDate -= interval;
			m_payoffDates[nDate--] = currentDate;
		}

		m_payoffDates[0] = nToday;
		fixings.resize(nAsset, 0, 0.0);
	}
	
}


void EuropeanBasketProduct::payout(CMatrix& values, CMatrix& pathArray, CVector& discountFactors, MlEqMonteCarlo& mc)
{			
	double basket = 0.0;
	for(int i=0; i<pathArray.rows(); i++)
	{
		basket += m_weight[i] * pathArray[i][m_nMaturity] / pathArray[i][m_nMaturity-1];
	}

	double fPayoff = std::max( m_cp*(basket - m_fStrike), 0.0) ;
	
	values[0][0] = fPayoff;
}

void EuropeanBasketProduct::setUp(CMatrix& value, MlEqMonteCarlo& mc)
{
	value.resize(1, 2);
	m_results.resize(1, 2);
	product::setUp(value,mc);
};

