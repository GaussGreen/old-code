
#include "stdafx.h"
#include "MlEqObjects.h"
#include "MomentumCliquetPricer.h"
#include "MonteCarlo.h"
#include "MlEqPde.h"


#undef max
#undef min


void GeneralCliquetPricer::initialize(MlEqAssetHandle				hUnderlying,
									const GVector<long>&		anPayoffDates,
									const GVector<bool>&		abStartsPeriod,			// [date]
									const CVector&				afStrikes,				// [date]
									const CVector&				afLocalFloors,			// [date]
									const CVector&				afLocalCaps,			// [date]
									const CVector&				afWeights,				// [date]
									const CVector&				afCoupons,				// [period]
									const CVector&				afGearings,				// [period]									
									const CVector&				afPeriodFloors,			// [period]
									const CVector&				afPeriodCaps,			// [period]					
									GVector<long>				anPayDatesArr,			// [period]	
									MinMaxEnum					minmaxtype,
									PayoffTypeEnum					potype,
									const CMatrix&				pastFixings)	

{
	m_payoffDates	=	anPayoffDates;

	int ndates = anPayoffDates.getsize()-1;
	if(afLocalCaps.getsize() != ndates)		throw "Local cap array has the wrong size";
	if(afLocalFloors.getsize() != ndates)	throw "Local floor array has the wrong size";
	if(afStrikes.getsize() != ndates)		throw "Strikes array has the wrong size";
	if(afWeights.getsize() != ndates)		throw "Weights array has the wrong size";

	m_afLocalCaps	=	afLocalCaps;
	m_afLocalFloors	=	afLocalFloors;
	m_afStrikes		=	afStrikes;
	m_afWeights		=	afWeights;

	m_nPeriods = anPayDatesArr.getsize();
	if(afPeriodCaps.getsize() != m_nPeriods)	throw "Global Cap array has the wrong size	";
	if(afPeriodFloors.getsize() != m_nPeriods)	throw "Global Floor array has the wrong size ";
	if(afCoupons.getsize() != m_nPeriods)		throw "Fixed Coupon array has the wrong size ";
	if(afGearings.getsize() != m_nPeriods)		throw "Gearings array has the wrong size ";

	m_period.resize(m_nPeriods);

	long nToday = hUnderlying->GetDateHandle()->GetDate();
	MlEqZeroCurveHandle hZeroCurve =  hUnderlying->GetPayZeroCurve(true); 

	for(int period = 0; period < m_nPeriods; period++)
	{
		m_period[period].mCap			 = afPeriodCaps[period];
		m_period[period].mFloor			 = afPeriodFloors[period];
		m_period[period].mFixedCoupon	 = afCoupons[period];
		m_period[period].mGearing		 = afGearings[period];

		long payDate = anPayDatesArr[period];
		m_period[period].mPayDate		 = payDate;
		m_period[period].mDiscountFactor 
			= (nToday >= payDate)?	0.0:hZeroCurve->GetDiscountFactor(nToday, payDate);	
	}

	// Set  periods
	int nPeriods = 0;
	m_period[0].mStart = 1;
	for (int n = 0; n < abStartsPeriod.getsize(); n++)
	{
		if (abStartsPeriod[n])
		{		
			if(nPeriods >= m_nPeriods)	throw"The fixing and coupon schedules sizes mismatch";

			m_period[nPeriods].mStart = n+1;
			if( nPeriods > 0)
				m_period[nPeriods-1].mEnd = n;

			++nPeriods;
		}
	}
	m_period[nPeriods-1].mEnd = abStartsPeriod.getsize(); // ndates-1

	if (nPeriods != m_nPeriods) throw "The number of pay dates must match the number of periods";

	switch( potype )
	{
	case NoPayoffType:
		m_localPayoff = &GeneralCliquetPricer::fForward;
		break;
	case Call:
		m_localPayoff = &GeneralCliquetPricer::fCall;
		break;
	case Put:
		m_localPayoff = &GeneralCliquetPricer::fPut;
		break;
	case Straddle:
		m_localPayoff = &GeneralCliquetPricer::fStraddle;
		break;
	}

	switch( minmaxtype )
	{
	case Minimum:
		m_periodPayoff = &GeneralCliquetPricer::fMinimum;
		break;
	case Maximum:
		m_periodPayoff = &GeneralCliquetPricer::fMaximum;
		break;
	default:
	//	m_periodPayoff = &GeneralCliquetPricer::fProduct;
	//	m_periodPayoff = &GeneralCliquetPricer::fSum;
		fPointToSpecificPayoff();
		break;
	}

	init_payout(pastFixings);
}

void GeneralCliquetPricer::init_payout(const CMatrix& pastPathArray)
{
	m_fPastPayoff = 0.0;
	m_nStartPeriod = 0;
	
	int nPastFixingsIndex = pastPathArray[0].getsize()-1;
	if( nPastFixingsIndex <= 0)	return;
							
	m_nStartPeriod = m_nPeriods;
	CVector perfArray;

	for(int nPeriod=0; nPeriod<m_nPeriods; nPeriod++)
	{
		long start = m_period[nPeriod].mStart;
		long end   = m_period[nPeriod].mEnd;
		if( nPastFixingsIndex < end )	// account only for fully past periods
		{
			m_nStartPeriod = nPeriod ;
			return;
		}

		perfArray.resize(end-start+1,0.0);
		int array_index = 0;

		for(int nDate=start; nDate<=end; nDate++)
		{
			double fReturn = fLocalReturn( pastPathArray[0][nDate] / pastPathArray[0][nDate-1], nDate-1) ;
			perfArray[array_index++] = fReturn;
		}
		double fPPayoff = fPeriodPayoff(perfArray, nPeriod);

		m_fPastPayoff += fPPayoff * m_period[nPeriod].mDiscountFactor; 
	}
}



void GeneralCliquetPricer::payout(CMatrix& values, CMatrix& pathArray, CVector& discountFactors, MlEqMonteCarlo& mc)
{			
								// date index (indexed s.t. the return is spot[nDate] / spot[nDate - 1])
	double	fPayoff = m_fPastPayoff;
	CVector perfArray;

	for(int nPeriod=m_nStartPeriod; nPeriod<m_nPeriods; nPeriod++)
	{
		long start = m_period[nPeriod].mStart;
		long end   = m_period[nPeriod].mEnd;
		perfArray.resize(end-start+1,0.0);
		int array_index = 0;

		for(int nDate=start; nDate<=end; nDate++)
		{
			double fReturn = fLocalReturn( pathArray[0][nDate] / pathArray[0][nDate-1], nDate-1) ;
			perfArray[array_index++] = fReturn;
		}
		double fPPayoff = fPeriodPayoff(perfArray, nPeriod);

		fPayoff += fPPayoff * m_period[nPeriod].mDiscountFactor;
	}
	
	values[0][0] = fPayoff;
}


void GeneralCliquetPricer::setUp(CMatrix& value, MlEqMonteCarlo& mc)
{
	value.resize(1, 2);
	m_results.resize(1, 2);
	product::setUp(value,mc);
}


double GeneralCliquetPricer::fLocalReturn(double perf, int dateindex)
{
	perf = (this->*m_localPayoff)(perf, dateindex) ;
	return m_afWeights[dateindex] * std::min(m_afLocalCaps[dateindex], std::max(m_afLocalFloors[dateindex], perf) );
}

double	GeneralCliquetPricer::fPeriodPayoff(const CVector& perfArray, int nPeriod)
{
	double periodPayoff = (this->*m_periodPayoff)( perfArray, nPeriod );
	periodPayoff = m_period[nPeriod].mFixedCoupon + m_period[nPeriod].mGearing * periodPayoff ;
	
	return std::min(m_period[nPeriod].mCap, std::max(m_period[nPeriod].mFloor, periodPayoff) ) ;
}



double	GeneralCliquetPricer::fSum(const CVector& perfArray, int nPeriod)
{
	double sum = 0.0;
	for(int i=0; i<perfArray.getsize(); i++)
		sum += perfArray[i];

	return sum;
}

double	GeneralCliquetPricer::fProduct(const CVector& perfArray, int nPeriod)
{
	double product = 1.0;
	for(int i=0; i<perfArray.getsize(); i++)
		product *= perfArray[i];

	return product;
}


double	GeneralCliquetPricer::fMinimum(const CVector& perfArray, int nPeriod)
{
	double mmin = std::numeric_limits<double>::max();
	for(int i=0; i<perfArray.getsize(); i++)
		mmin = std::min(perfArray[i], mmin);

	return mmin;
}


double	GeneralCliquetPricer::fMaximum(const CVector& perfArray, int nPeriod)
{
	double mmax = std::numeric_limits<double>::min();
	for(int i=0; i<perfArray.getsize(); i++)
		mmax = std::max(perfArray[i], mmax);

	return mmax;
}

void GeneralCliquetPricer::fPointToSpecificPayoff()
{
	m_periodPayoff = &GeneralCliquetPricer::fSum ;
} 


