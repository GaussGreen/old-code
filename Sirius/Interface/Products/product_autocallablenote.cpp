// product_autocallableNote.cpp : Implementation of CAutocallableNote
#include "stdafx.h"
#include "Products.h"
#include "product_autocallableNote.h"
#include "MlEqParameterList.h"
#include "MonteCarlo.h"


STDMETHODIMP CAutocallableNote::Evaluate(IResult* pVal)
{
	begin_function	

	HRESULT hr;

	if (!m_hUnderlying) throw "No underlying defined";
	if (!m_hZeroCurve)	throw "No zero curve defined";
	if (!m_hSpread)	throw "No spreads defined";
	if (!m_hWeights)	throw "No weights defined";
	if (!m_hResetSchedule)
		throw "No reset dates schedule defined";

	if (!m_hPaySchedule)
		throw "No pay dates schedule defined";

	
	GVector<long>	resetDates;			
	m_hResetSchedule->GetDates(resetDates);
	int ndates = resetDates.getsize();

	if( m_hSpread->getsize() != ndates-1 )	throw "wrong spread array size";
	if( m_hWeights->getsize() != ndates-1 )	throw "wrong weights array size";


	GVector<long>	payDates;			
	m_hPaySchedule->GetDates(payDates);
	int nPayDates = payDates.getsize();
	if( nPayDates != ndates-1 )	throw "wrong pay dates array size";

	long nToday = m_hUnderlying->GetDateHandle()->GetDate();

	MlEqDateHandle hDate = new MlEqDate( nToday, m_dcc );

	CVector libor_coupons(ndates-1);

	int npast = (!!m_hLiborFixings)? m_hLiborFixings->getsize() : 0 ;

	for(int i=1; i<ndates; ++i)
	{
		double dt = hDate->GetYearFraction(resetDates[i]) - hDate->GetYearFraction(resetDates[i-1]) ;

		if( i <= npast && (*m_hLiborFixings)[i-1] )
		{
			double libbborrr = (*m_hLiborFixings)[i-1];
			libor_coupons[i-1] = ( (*m_hWeights)[i-1] * libbborrr + (*m_hSpread)[i-1] ) * dt ;
		}
		else
		{
			if( resetDates[i-1] < nToday ){
				throw "You must fill the past fixings in" ;
			}
			double dff = m_hZeroCurve->GetDiscountFactor( resetDates[i-1], resetDates[i] );
			double libbborrr = (1./dff - 1.0) ;
			libor_coupons[i-1] = (*m_hWeights)[i-1] * libbborrr + (*m_hSpread)[i-1] * dt;
		}
	}

	double cp;
	if(m_cpe == Call){
		cp = 1.;
	}
	else if (m_cpe == Put){
		cp = -1.;
	}
	else	throw "Wrong payoff type";



	// start the monte carlo stuff...

	GVector<long>		barrierDates;		
	CMatrix				mFixings;
	CMatrix				mResults;

	long	nWidth = m_hBarrierSchedule->GetWidth();	

	barrierDates = m_hBarrierSchedule->GetDates();
	int nbarrierdates = barrierDates.getsize();

	if(!m_hThresholds)	throw "No thresholds defined";
	if( m_hThresholds->getsize() !=  nbarrierdates-1 )	throw "wrong thresholds array size";

	if(!m_hCouponUp)	throw "No coupon up defined";
	if( m_hCouponUp->getsize() !=  nbarrierdates-1 )	throw "wrong coupon up array size";

	if(!m_hCouponLow)	throw "No coupon low defined";
	if( m_hCouponLow->getsize() !=  nbarrierdates-1 )	throw "wrong coupon low array size";


	if( !m_hUnderlying->IsBasket() )
	{
		// set up the spot fixings (mFixings)
		if (!nWidth)
		{	
			std::vector<long>			afDates;
			m_hBarrierSchedule->GetDates(afDates);
			m_hUnderlying->GetSpotsBeforeToday(nToday + 1, afDates, &mFixings);

		}
		else if (nWidth == 1)	{
			m_hBarrierSchedule->GetColumnBeforeToday(nToday + 1, 0, mFixings);
		} 
		else {
			throw "The fixing schedule contains too many columns";
		}
	}
	else
	{
		if (!nWidth)
		{
			// fixing schedule contains just dates - we get the values from the basket constituent spot schedules
			std::vector<long>			afDates;
			m_hBarrierSchedule->GetDates(afDates);	
			
			std::vector<double>	vSpots;
			m_hUnderlying->m_assets[0]->GetSpotsBeforeToday(nToday + 1, afDates, &vSpots);

			mFixings.resize( m_hUnderlying->m_assets.size(), vSpots.size(), 0.0);	
			
			for (long nAsset = 0; nAsset < m_hUnderlying->m_assets.size(); nAsset++)
			{
				m_hUnderlying->m_assets[nAsset]->GetSpotsBeforeToday(nToday + 1, afDates, &vSpots);

				for (long nSpot = 0; nSpot < vSpots.size(); nSpot++){
					mFixings[(int)nAsset][(int)nSpot] = vSpots[nSpot];
				}
			}
		} 
		else if (nWidth == m_hUnderlying->m_assets.size())		{
				// fixing schedule contains both the dates and the fixings
				m_hBarrierSchedule->GetColumnsBeforeToday(nToday + 1, mFixings);
		}
		else {
			throw "The number of columns in the fixing schedule is invalid - you need to pass as many data columns as there are assets in the basket";
		}
	}



	if(m_hUnderlying->IsBasket())
	{
		if( !m_hRainbowWeights ){
			throw "No rainbow weighting defined";
		}
		if( m_hRainbowWeights->getsize() != m_hUnderlying->GetAssets().size() ){
			throw "The rainbow weights array has the wrong size";
		}

	}
	else
	{
		m_hRainbowWeights = new MlEqArray;
		m_hRainbowWeights->resize(1, 1.0);
	}

	initialize( *m_hRainbowWeights,
				*m_hThresholds,
				*m_hCouponUp,
				*m_hCouponLow,
				*m_hCallLevels, 
				*m_hRebates, 
				barrierDates,  
				m_fStrike, cp, m_fGearing, m_fFixedCoupon, 
				payDates, 
				resetDates, 
				libor_coupons, 
				m_hZeroCurve, 
				hDate, 
				mFixings);

	
	// now choose the model...
	// set up the Monte Carlo

	CVector						vModelInfo;
	MlEqMonteCarloHandle		hMonteCarlo;

	int nModelInfoSize = m_hParameterList->GetCount();

	if( nModelInfoSize == 16 )
	{

		vModelInfo.resize(16);
		m_hParameterList->BeginExhaust();
		m_hParameterList->GetValue("cL, cR, addTanhWing, yPower, seed, npaths, calibflag, numberVolStates, localVolFlag, savevolvolinfo, numberGridPoints, saveVolGridFlag, randomNumberFlag, ControlVariate, GlobalCalib, UseHermite = 0", ",", "=", &vModelInfo[0]);
		m_hParameterList->CheckUnexhausted();

		if(m_hUnderlying->IsBasket())
		{
			CMultiAssetForwardSkewMC* pMC;
			
			hMonteCarlo = pMC = new CMultiAssetForwardSkewMC(m_hUnderlying->GetDateHandle());
			pMC->Initialize(m_hUnderlying->m_assets, *this, mFixings, vModelInfo, CMatrix(), CMatrix(), std::vector<std::vector<MlEqStrikeHandle> >());
	
			hMonteCarlo->generatePaths();
			hMonteCarlo->simulate(mResults, *this);
		}
		else
		{
			MlEqStochBetaVolatilityStructureHandle	hVolatility = dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*m_hUnderlying->GetVolatilityStructure());
			if (!hVolatility) throw "You need to supply a stochastic beta volatility matrix";

			CForwardSkewMC* pMC;


			hMonteCarlo = pMC = new CForwardSkewMC(hVolatility->getDateToDouble());
			pMC->Initialize(*m_hUnderlying, *this, *hVolatility, mFixings, vModelInfo, CMatrix(), CMatrix(), vector<vector<MlEqStrikeHandle> >());

			hMonteCarlo->generatePaths();
			hMonteCarlo->simulate(mResults, *this);
		}
	}
	else if( nModelInfoSize == 3 )
	{

		CVector vModelInfo(3);
		m_hParameterList->BeginExhaust();
		m_hParameterList->GetValue("npaths = 50000; numberStepsAYear = 50; randomNumberFlag = 6", ";", "=", &vModelInfo[0]);
		m_hParameterList->CheckUnexhausted();

		int npaths	= vModelInfo[0];
		int nSteps	= vModelInfo[1];
		int rngFlag = vModelInfo[2];


		if(m_hUnderlying->IsBasket())
		{
			CMultiLocalVolMC* pMC;			
			hMonteCarlo = pMC = new CMultiLocalVolMC(m_hUnderlying->GetDateHandle(), m_hUnderlying);
			pMC->CLocalVolMC::initialize(*this, mFixings, npaths, rngFlag, nSteps);
			hMonteCarlo->simulate(mResults, *this);
		}
		else
		{
			CLocalVolMC* pMC;
			
			hMonteCarlo = pMC = new CLocalVolMC(m_hUnderlying->GetDateHandle(), m_hUnderlying);
			pMC->initialize(*this, mFixings, npaths, rngFlag, nSteps);
			hMonteCarlo->simulate(mResults, *this);
		}
	}
	else{
		throw "Wrong Model";
	}

	double price = mResults[0][0] ;

		// add the current coupon to pay if necessary	

	for(int t=0; t<nPayDates; ++t)
	{
		if( resetDates[t+1] < nToday )
		{
			long payDate = payDates[t];

			if(  payDate >= nToday ){
				price += libor_coupons[t] * m_hZeroCurve->GetDiscountFactor( payDate );
			}
		}
		else break;
	}


	begin_calculate(pVal, L"Price")
		if (hr = pVal->AddValue(L"Price", CComVariant( price ))) return hr;
	end_calculate


	end_function
}




HRESULT CAutocallableNote::FinalConstruct(void)
{	
	m_fStrike = 0.0;
	m_fGearing = 1.0;
	m_cpe = Call;
	m_dcc = Actual365;
	m_fFixedCoupon = 0.0;

	return S_OK;
}

STDMETHODIMP CAutocallableNote::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IARPropFitVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

