
#include "stdafx.h"
#include "CVarianceSwap.h"


CVarianceSwap::CVarianceSwap(	MlEqAssetHandle			hUdly,
								DATE					dateStart,
								DATE					dateMaturity,
								MlEqDateScheduleHandle	hFixingSchedule,
								std::string				szCalendar,
								PerformanceTypeEnum		perfType,
								bool					useCurrentSpot,
								StrikesTypeEnum			boundType,
								double					cutUp,
								double					cutDown,
								long					nDaysAYear):
CStaticReplication( hUdly, dateMaturity ),
m_dateStart(dateStart),
m_hFixingSchedule(hFixingSchedule),
m_szCalendar(szCalendar),
m_perfType(perfType),
m_boundType(boundType),
m_strikeUp(cutUp),
m_strikeDown(cutDown),
m_useCurrentSpot(useCurrentSpot),
m_nDaysAYear(nDaysAYear)
{

	MlEqDateHandle hDate = new MlEqDate( m_dateStart );
	m_nDaysStartMat	= hDate->GetNumBusinessDays(m_dateStart, long(m_dateMaturity), m_szCalendar ) ; 
	m_nDaysTodayMat = hDate->GetNumBusinessDays(long(m_nToday), long(m_dateMaturity), m_szCalendar ) ;
	m_nDaysTodayStart = fabs(m_nDaysStartMat - m_nDaysTodayMat) ;

	m_adjTodayStart	 = (m_nDaysTodayStart < 1e-6)? 1.: 2.* m_nDaysAYear/m_nDaysTodayStart;
	m_adjTodayMat	 = (m_nDaysTodayMat < 1e-6)? 1.: 2.* m_nDaysAYear/m_nDaysTodayMat ;
	m_adjStartMat	 = (m_nDaysStartMat < 1e-6)? 1.: 2.* m_nDaysAYear/m_nDaysStartMat ;

	m_forwardStart = m_spot;
	if( m_dateStart > m_nToday ){
		m_forwardStart = m_hUnderlying->GetQuantoForward(m_nToday, m_dateStart, false);	
	}

	switch( m_boundType )
	{
	case NoStrikesType:
		{
			m_strikeUp = 1e30;
			m_strikeDown = 0.0;
			return;
		}
	case Fixed:
		break;
	case ForwardBased:
		{
			m_strikeDown	*= m_forward;
			m_strikeUp		*= m_forward;	
		}
		break;
	case Normalised:
		{
			double vol	= m_hUnderlying->GetVolatility(MlEqStrike(m_forward),m_hUnderlying->GetDateHandle()->GetDate(),m_dateMaturity); 

			m_strikeUp		= m_forward*exp(vol*sqrt(m_mat) * m_strikeUp );
			m_strikeDown	= m_forward*exp(vol*sqrt(m_mat) * m_strikeDown );
		}
		break;
	default:
		m_strikeUp = m_strikeDown = m_forward = 0.0;
		break;
	}

	m_integrationLimits = CVector(3, 0.0);
	m_integrationLimits[0] = m_strikeDown ;
	m_integrationLimits[1] = m_splitLimit = m_forward ;
	m_integrationLimits[2] = m_strikeUp	;

	m_infiniteIntegrate = false;
}



double CVarianceSwap::f(double x)
{ 
	return -log(x/m_forward);
}

double CVarianceSwap::df(double x)
{
	return -1/x ;
}

double CVarianceSwap::ddf(double x)
{
	return 1/(x*x);
}




void CVarianceSwap::fComputePrice(CMatrix& results)
{	
	double startVar = 0.0, futVar = 0.0, value = 0.0;

	// first get the future variance --> maturity

	futVar = CStaticReplication::fComputePrice();

	value += futVar;

	if( m_dateStart > m_nToday )	// compute another variance swap
	{
		double fwd = m_forward ;
		long maturity = m_dateMaturity;
		updateMaturity( m_dateStart, m_forwardStart );
		startVar = CStaticReplication::fComputePrice();	
		value -= startVar ;
		updateMaturity( maturity, fwd );
	}
	else if( m_dateStart < m_nToday )
	{
		startVar = fComputePastVariance();
		value += startVar ;
	}
	else
	{
		startVar = 0.0;
	}

	startVar *= m_adjTodayStart ;
	futVar	 *= m_adjTodayMat ;
	value    *= m_adjStartMat ;

// fill results

	results.resize(4, 2, 0.0);

	results[0][0] = value;		results[0][1] = m_nDaysStartMat;
	results[1][0] = startVar;	results[1][1] = m_nDaysTodayStart ;
	results[2][0] = futVar;		results[2][1] = m_nDaysTodayMat ;

	results[3][0] = m_strikeUp; results[3][1] = m_strikeDown;
}


void CVarianceSwap::updateMaturity( DATE maturity, double forward )
{
	m_dateMaturity			= maturity;

	MlEqDateHandle	hDate	= new MlEqDate(*m_hUnderlying->GetDateHandle());
	m_mat					= hDate->GetYearFraction(m_dateMaturity) ;

	m_forward				= forward;	
	m_integrationLimits[1]	= m_splitLimit = m_forward;
}


double CVarianceSwap::fComputePastVariance()
{
	std::map< long, std::vector<double> >::const_iterator it = (*m_hFixingSchedule).begin();
	
	double pastVar = 0.0, spotPrev = 0.0;

	for(; it  != (*m_hFixingSchedule).end(); it++)
	{
		long date = it->first ;

		if((it->second).size() == 0 && date < long(m_nToday) )
				throw"The fixings are not properly populated";

		if( date > long(m_dateStart) && date <= long(m_nToday) )
		{
			double spot;
			if( date == long(m_nToday) && ((it->second).size() == 0 || m_useCurrentSpot) )
				spot = m_spot;
			else if( (it->second).size() > 0 )
				spot = (it->second)[0];


			double perf = spot/spotPrev ;

			if( m_perfType == LogarithmicType ){
				perf = log(perf);			
			}
			else{
				perf = perf - 1.0;
			}
			
			pastVar += perf*perf ;
			spotPrev = spot;
		}
		else if( date == long(m_dateStart) && date < long(m_nToday))
		{
			spotPrev = (it->second)[0] ;
		}
		else break;
	}

	return pastVar * 0.5 ;	// because of time adjustment
}






















/*


void CVarianceSwap::shiftSpot( double factor )
{
	m_spot			*= factor;
	m_forward		*= factor;
	m_forwardStart	*= factor;

	m_integrationLimits[1]	= m_splitLimit = m_forward;
}



double CVarianceSwap::fComputePrice(double& startVar, double& futVar)
{	
	double value = 0.0;

	// first get the future variance --> maturity

	futVar = CStaticReplication::fComputePrice();

	value += futVar;

	if( m_dateStart > m_nToday )	// compute another variance swap
	{
		double fwd = m_forward ;
		long maturity = m_dateMaturity;
		updateMaturity( m_dateStart, m_forwardStart );
		startVar = CStaticReplication::fComputePrice();	
		value -= startVar ;
		updateMaturity( maturity, fwd );
	}
	else if( m_dateStart < m_nToday )
	{
		startVar = fComputePastVariance();
		value += startVar ;
	}
	else
	{
		startVar = 0.0;
	}


	startVar *= m_adjTodayStart ;
	futVar	 *= m_adjTodayMat ;
	value    *= m_adjStartMat ;

	return value ;
}


void CVarianceSwap::fComputePriceAndGreeks(CMatrix& result)
{
	result.resize(5, 3, 0.0);
	double eps = 0.005;
	double startVar, futVar;

	double price = fComputePrice( startVar, futVar );

//	double fixVol = sqrt(futVar);	

	result[0][0]	 = price;
	result[0][1]	 = startVar;
	result[0][2]	 = futVar;

	shiftSpot( 1+eps );

	price = fComputePrice( startVar, futVar );

	result[1][0]	 = price;
	result[1][1]	 = startVar;
	result[1][2]	 = futVar;

	shiftSpot( (1-eps)/(1+eps) );

	price = fComputePrice( startVar, futVar );

	shiftSpot( 1/(1-eps) );

	eps				*= m_spot ;

	// calculate gamma

	result[2][0]	 = (result[1][0] + price	-2.* result[0][0]) / (eps*eps) ;
	result[2][1]	 = (result[1][1] + startVar -2.* result[0][1]) / (eps*eps) ;
	result[2][2]	 = (result[1][2] + futVar	-2.* result[0][2]) / (eps*eps) ;

	// calculate delta

	result[1][0]	 = (result[1][0] - price )		/ (2.*eps) ;
	result[1][1]	 = (result[1][1] - startVar )	/ (2.*eps) ;
	result[1][2]	 = (result[1][2] - futVar )		/ (2.*eps) ;

	// room for other greeks...


	result[3][0] = m_nDaysStartMat;
	result[3][1] = m_nDaysTodayStart ;
	result[3][2] = m_nDaysTodayMat ;
	result[4][0] = m_nDaysAYear ;

	result[4][2] = m_strikeUp;
	result[4][1] = m_strikeDown;
}



double CVarianceSwap::fComputeDiscreteVariance(long& nbBusDays, long dateEnd)
{
	long dateStart = std::max( m_dateStart, m_nToday );
	dateEnd   = std::min( m_dateMaturity, dateEnd );

	GVector< long > fixingDates;
	m_hFixingSchedule->GetDates(fixingDates);
	MlEqDateHandle	hDate = new MlEqDate(*m_hUnderlying->GetDateHandle());	

	RCPtr<CVarianceSwapDiscreteSlice> slice = new CVarianceSwapDiscreteSlice(m_hUnderlying, dateEnd);

	double discreteVar = 0.0;
	double fwdNext = m_hUnderlying->GetQuantoForward(hDate->GetDate(), dateStart, false);


	discreteVar -= slice->fComputeSquareTerm( dateStart, fwdNext );

	nbBusDays = 0;

	for(int k=0; k<fixingDates.getsize(); k++)	// cross terms
	{
		if( fixingDates[k+1] <= dateEnd && fixingDates[k] >= dateStart )
		{
			double fwdStart	= fwdNext;
			fwdNext	= m_hUnderlying->GetQuantoForward(hDate->GetDate(), fixingDates[k+1], false);

			discreteVar += slice->fComputeCrossTerm(fixingDates[k], fixingDates[k+1], fwdStart, fwdNext) ;
			nbBusDays++ ;
		}
		else if( fixingDates[k] >= dateEnd)	break;
	}
	

	discreteVar += slice->fComputeSquareTerm( dateEnd, fwdNext );

	return discreteVar ;
}



void CVarianceSwap::fComputeDiscretePrice(CMatrix& results, long discreteDate)
{	
// first case : t < td < ts < T

	if( discreteDate <= m_dateStart )
	{
		fComputePrice(results);
		results[3][0] = 0.0;;	
		return ;
	}

	long nbBusDays;
	double discreteVar = fComputeDiscreteVariance(nbBusDays, discreteDate);
	double adjDiscrete = m_nDaysAYear / double(nbBusDays) ;

	double futVar = 0.0, pastVar = 0.0;

// case 2: td > T
	if( discreteDate >= m_dateMaturity )
	{
		results.resize(4, 2, 0.0);
		results[0][0] = discreteVar * adjDiscrete ;
		return;
	}

// case 2:	 ts < td < T

	futVar = CStaticReplication::fComputePrice();

	if( discreteDate < m_dateMaturity )
	{
		double fwd = m_forward ;
		long maturity = m_dateMaturity;
		MlEqDateHandle	hDate	= new MlEqDate(*m_hUnderlying->GetDateHandle());
		double forwardStart =  m_hUnderlying->GetQuantoForward(hDate->GetDate(), discreteDate, false);
		updateMaturity( discreteDate, forwardStart );
		double startVar = CStaticReplication::fComputePrice();	
		futVar -= startVar ;
		updateMaturity( maturity, fwd );
		m_adjTodayMat = 2.*m_nDaysAYear/double(m_nDaysTodayMat - nbBusDays) ;

		// now calculate past variance:
		if( m_dateStart < m_nToday )
		{
			pastVar = fComputePastVariance();
		}

		double totalVar = pastVar + 0.5*discreteVar + futVar ;

		pastVar		 *= m_adjTodayStart ;
		futVar		 *= m_adjTodayMat ;
		discreteVar	 *= adjDiscrete ;
		totalVar	 *= m_adjStartMat ;

		results.resize(4, 2, 0.0);

		results[0][0] = totalVar;	
		results[1][0] = pastVar;	
		results[2][0] = futVar;	
		results[3][0] = discreteVar;	
	}

}





void CVarianceSwap::CalibrateVanilla(const CVector& strikes, const CVector& cp, CVector& weights)
{

	int nVanilla = strikes.getsize() ;
	int nCall = 0;
	int i;
	double one_over_forward = 1./m_forward;

	for(i=0; i<nVanilla; i++){
		if( cp[i] > 0.0 ) nCall++;
	}
	int nPut = nVanilla - nCall ;

	// assumes strikes are properly ordered...
	weights.resize(nVanilla, 0.0);
	double runningSum_c = 0.0;

	for(i=nPut; i<nVanilla-1; i++)
	{
		weights[i] = one_over_forward-log(strikes[i+1]/strikes[i])/(strikes[i+1] - strikes[i]) - runningSum_c ;
		runningSum_c += weights[i];
	}
//	weights[nVanilla-1] = one_over_forward - 1./strikes[nVanilla-1] ;
//	weights[nVanilla-1] = one_over_forward - 0.5*log(2.)/strikes[nVanilla-1] - runningSum_c ;

	double runningSum_p = 0.0;
	for(i=nPut-1; i>0; i--)
	{
		weights[i] = - one_over_forward + log(strikes[i+1]/strikes[i])/(strikes[i+1] - strikes[i]) - runningSum_p ;
		runningSum_p += weights[i];
	}
//	weights[0] = -one_over_forward + 1./strikes[0];
//	weights[0] = - one_over_forward + 2.*log(2.)/strikes[0] - runningSum_p ;


	double strikeUp, strikeDown;
	strikeUp = strikes[nVanilla-1];
	strikeUp *= strikeUp * one_over_forward ;

	double vol	= m_hUnderlying->GetVolatility(MlEqStrike(m_forward),m_hUnderlying->GetDateHandle()->GetDate(),m_dateMaturity); 
	strikeDown = strikes[0] * exp(-vol*sqrt(m_mat)) ;
	strikeDown = one_over_forward * strikes[0]*strikes[0];

	weights[nVanilla-1] = one_over_forward - log(strikeUp/strikes[nVanilla-1])/(strikeUp - strikes[nVanilla-1]) - runningSum_c ;
	weights[0] = - one_over_forward + log(strikes[0]/strikeDown)/(strikes[0] - strikeDown) - runningSum_p;


	double price = CalculateVanillaPortfolio( strikes, cp, weights );

	return;
}



double CVarianceSwap::CalculateVanillaPortfolio(const CVector& strikes, const CVector& cp, const CVector& weights)
{
	double price = 0.0;

	int nVanilla = strikes.getsize() ;

	for(int i=0; i<nVanilla; i++)	
	{
		price += weights[i] * vanillaPrice( strikes[i], cp[i] );
	}

	return price * m_adjTodayMat ;
}


















double	CVarianceSwapDiscreteSlice::fComputeSquareTerm( DATE dateMaturity, double forward )
{
	m_dateMaturity			= dateMaturity;

	MlEqDateHandle	hDate	= new MlEqDate(*m_hUnderlying->GetDateHandle());
	m_mat					= hDate->GetYearFraction(m_dateMaturity) ;

	if( m_mat < 1e-4 )		return log(m_spot)*log(m_spot) ;


	m_forward				= forward;	
	m_integrationLimits[1]	= m_splitLimit = m_forward;

	m_isCrossTerm = false;

	double sqterm = CStaticReplication::fComputePrice() ;
	return sqterm ;
}




double	CVarianceSwapDiscreteSlice::fComputeCrossTerm( DATE dateStart, DATE dateMaturity, double forwardStart, double forwardMat )
{
	m_dateMaturity			= dateStart;

	MlEqDateHandle	hDate	= new MlEqDate(*m_hUnderlying->GetDateHandle());
	m_mat					= hDate->GetYearFraction(m_dateMaturity) ;

	if( m_mat < 1e-4 )	// from today to next step
	{
		double localvol  = m_hVolatilityStructure->getFutureVol(MlEqStrike(m_spot), m_nToday, dateStart, m_spot);
		double dt		 = hDate->GetYearFraction(dateMaturity);

		double temp = log(m_spot) * ( localvol*localvol*dt - 2.*log(forwardMat/forwardStart) );
		return temp ;
	}

	m_forward				= forwardStart;	
	m_integrationLimits[1]	= m_splitLimit = m_forward;

	m_isCrossTerm = true;

	m_integrateVolTerm = false;
	double crossterm =	CStaticReplication::fComputePrice() * (-2.*log(forwardMat/forwardStart) ) ;

	m_integrateVolTerm = true;

	double ddt = - CStaticReplication::fComputePrice() ;
//	crossterm -= CStaticReplication::fComputePrice() ;

//	m_dateMaturity = dateMaturity ;
//	m_mat = hDate->GetYearFraction(dateMaturity);

	double dt_ = hDate->GetYearFraction(dateMaturity) - m_mat;
	m_mat = hDate->GetYearFraction(dateStart+1);
	m_dateMaturity = dateStart+1 ;
	
//	crossterm += CStaticReplication::fComputePrice() ;
	ddt +=  CStaticReplication::fComputePrice() ;
	ddt *= dt_*365. ;

	crossterm += ddt;

	return crossterm ;
}


double CVarianceSwapDiscreteSlice::f( double x )
{
	if( m_isCrossTerm )
	{
		if( m_integrateVolTerm )	return 0.0;
		return log( x );
	}

	double temp = log( x );
	return temp*temp ;
}


double CVarianceSwapDiscreteSlice::df( double x )
{
	if( m_isCrossTerm )
	{
		if( m_integrateVolTerm )	return 0.0;
		return 1./x;
	}

	double temp = 2.*log( x );
	return temp*x ;
}


double CVarianceSwapDiscreteSlice::ddf( double x )
{
	double temp = 0.0;

	if( m_isCrossTerm )
	{
		if( m_integrateVolTerm )
		{
			temp = 2.*log(x)/(x*x) ;
			return temp;
		}
		return -1./(x*x) ;			
	}
	else
	{
		temp = 1.- log(x);
		return 2*temp/(x*x) ;
	}

	return 0.0;
}










/*
double	CVarianceSwapDiscreteSlice::fDupireLocalVariance( double strike, double maturity)
{
	// just a basic Dupire formula implementation...

	double y = log(strike/m_forward) ;
	
	double vol  = m_hVolatilityStructure->getFutureVol(MlEqStrike(strike), m_nToday, m_dateMaturity, m_spot);

	double w = vol*vol*m_mat ;

	double eps = 0.005;
	double dwdy = w;
	double ddwddy = -2*w ;

	vol = m_hVolatilityStructure->getFutureVol(MlEqStrike(strike*exp(eps)), m_nToday, m_dateMaturity, m_spot);

	double dw = vol*vol*m_mat ;
	dwdy = dw ;
	ddwddy += dw ;

	vol = m_hVolatilityStructure->getFutureVol(MlEqStrike(strike*exp(-eps)), m_nToday, m_dateMaturity, m_spot);

	dw = vol*vol*m_mat ;
	dwdy -= dw ;
	ddwddy += dw ;

	eps *= y;
	dwdy /= 2.* eps;
	ddwddy /= eps*eps ;

	vol = m_hVolatilityStructure->getFutureVol(MlEqStrike(strike), m_nToday, m_dateMaturity+1, m_spot);

	double dwdt = vol*vol*(m_mat+1./365.) - w ;		// check this out !!!
	dwdt /= 1./365.;

	double num = w*w * dwdt ;

	double denom = w* ( - y*dwdy + w*( 1. + 0.5*ddwddy ) ) ;
	denom -= 0.25 * dwdy*dwdy * ( y*y + w + 0.25*w*w ) ;

	if( num < 1e-6 || denom < 1e-6 )	return 1e-4;

	double vl2 = num / denom ;

	vl2 = std::max(1e-4, std::min( vl2, 1.)) ;

	return vl2;
}
*/




