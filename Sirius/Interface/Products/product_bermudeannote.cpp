// product_bermudeanNote.cpp : Implementation of CBermudeanNote
#include "stdafx.h"
#include "Products.h"
#include "product_bermudeanNote.h"
#include "MlEqParameterList.h"
#include "pdeproducts.h"
#include "mleqpde.h"
/////////////////////////////////////////////////////////////////////////////
// CBermudeanNote


STDMETHODIMP CBermudeanNote::Evaluate(IResult* pVal)
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
	double spot = m_hUnderlying->GetSpot(nToday);

	MlEqDateHandle hDate = new MlEqDate( nToday, m_dcc );

	CVector libor_coupons(ndates-1);

	int npast = (!!m_hLiborFixings)? m_hLiborFixings->getsize() : 0 ;

	for(int i=1; i<ndates; ++i)
	{
		double dt = hDate->GetYearFraction(resetDates[i]) - hDate->GetYearFraction(resetDates[i-1]) ;

		if( i <= npast && (*m_hLiborFixings)[i-1] )
		{
			double libbborrr = (*m_hLiborFixings)[i-1];
			libor_coupons[i-1] = ( (*m_hWeights)[i-1] * libbborrr + (*m_hSpread)[i-1] ) * dt;
		}
		else
		{
			if( resetDates[i-1] < nToday ){
				throw "You musy fill the past fixings in" ;
			}
			double dff = m_hZeroCurve->GetDiscountFactor( resetDates[i-1], resetDates[i] );
			
			double libbborrr = (1./dff - 1.0) ;
			libor_coupons[i-1] = (*m_hWeights)[i-1] * libbborrr + (*m_hSpread)[i-1] * dt;
		}
	}


	double volbump = 0.01;
	DupireLocalVolHandle lv = new DupireLocalVol;

	double mat = hDate->GetYearFraction(m_dateMaturity);
	double fwd = m_hUnderlying->GetForward(m_dateMaturity, false);
	double vol = m_hUnderlying->GetVolatility(MlEqStrike(fwd),m_dateMaturity);

	double fnstdev  = (mat < 1.)?6.0 - 2.0 * sqrt(mat) : 4.0;
	double nsteps = (mat < 1.5)?250 - 200 * sqrt(mat/1.5) : 50;
	int nt = int(nsteps * mat);
	int nx = nt ;

	double highSpot = spot*exp(fnstdev*vol*sqrt(mat));
	double lowSpot  = spot*exp(- 1.5*fnstdev*vol*sqrt(mat));


	// build the whole payoff dates set
	GVector<long> tmpDates = resetDates;

	GVector<long> vmat(1);
	vmat[0] = m_dateMaturity ;

	if( !m_hCallDates )		throw "No call dates defined";
	if( !m_hRebates )		throw "No rebates defined";

	GVector<long>	callDates;			
	m_hCallDates->GetDates(callDates);

	merge( tmpDates, resetDates, vmat );
	merge( tmpDates, tmpDates, callDates );

	lv->initialize(m_hUnderlying, lowSpot, highSpot, tmpDates, nt, nx, volbump );


	RCPtr<callableNote01Helper> callableNote = new callableNote01Helper();
	GVector<long> pdeDates = lv->getDates();

	double cp;
	if(m_cpe == Call){
		cp = 1.;
	}
	else if (m_cpe == Put){
		cp = -1.;
	}
	else	throw "Wrong payoff type";

	callableNote->initialize( spot, 
							m_fRefSpot, 
							hDate, 
							libor_coupons, 
							resetDates, 
							*m_hRebates,
							callDates,
							pdeDates,
							cp, m_fStrike, m_fGearing, m_fFixedCoupon);


	
	pdeLocalVol pde;
	pde.init( lv,&(*callableNote));

	CVector result(4);
	CVector tmp(1);

	pde.pde_integrate(tmp);
	double price = tmp[0];
	result[0] = price;

	double delta, gamma ;
	callableNote->greeks( delta, gamma, 0, 0, 0, 0, &pde );
	result[1] = delta ;
	result[2] = gamma ;
	result[2] *= 0.01 * spot ;

	lv->parallelVegaBump();			// vega
	pde.init( lv,&(*callableNote));
	pde.pde_integrate(tmp);
	result[3] = tmp[0] - price ;


	// add the current coupon to pay if necessary	

	for(int t=0; t<nPayDates; ++t)
	{
		if( resetDates[t+1] < nToday )
		{
			long payDate = payDates[t];

			if(  payDate >= nToday ){
				price += libor_coupons[t+1] * m_hZeroCurve->GetDiscountFactor( payDate );
			}
		}
		else break;
	}


	begin_calculate(pVal, L"Price")
		if (hr = pVal->AddValue(L"Price", CComVariant(price))) return hr;
	end_calculate

	begin_calculate(pVal, L"Delta")
		if (hr = pVal->AddValue(L"Delta", CComVariant(delta))) return hr;
	end_calculate

	begin_calculate(pVal, L"Gamma")
		if (hr = pVal->AddValue(L"Gamma", CComVariant(result[2]))) return hr;
	end_calculate

	begin_calculate(pVal, L"Vega")
		if (hr = pVal->AddValue(L"Vega", CComVariant(result[3]))) return hr;
	end_calculate

	end_function
}




HRESULT CBermudeanNote::FinalConstruct(void)
{	
	m_fStrike = 0.0;
	m_dateMaturity = 0.0;
	m_fGearing = 1.0;
	m_cpe = Call;
	m_dcc = Actual365;
	m_fFixedCoupon = 0.0;
	m_fRefSpot = 0.0;

	return S_OK;
}

STDMETHODIMP CBermudeanNote::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IARPropFitVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}



