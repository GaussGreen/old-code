
#include "stdafx.h"
#include "MlEqObjects.h"
#include "CallableCliquetSwapPricer.h"
#include "MonteCarlo.h"
#include "MlEqSwap.h"
#include <set>
#include "MlEqPde.h"



void CallableCliquetSwapPricer::initialize(	CVector	&localCaps,
												CVector	&localFloors,
												CVector	&strikes,
												CVector	&callput,
												CVector	&weights,
												CVector &GlobalCap,
												CVector &GlobalFloor,
												CVector &FixedCoupon,
												double	callableLevel,
												GVector< long >& fixingDates,
												CVector	&periodId,
												GVector< long >& couponPayDates,
										//		MlEqZeroCurveHandle hZeroCurve,
										//		GVector< long >& fundingPayDates,
										//		double FundingValue,
										//		MlEqSwapHandle hswap,
												CMatrix& pastFixings)

{
	m_payoffDates	=	fixingDates;
	mCallableLevel	=	callableLevel ;

	int ndates = fixingDates.getsize();
	if(localCaps.getsize() != ndates)	throw "Local cap array has the wrong size";
	if(localFloors.getsize() != ndates)	throw "Local floor array has the wrong size";
	if(strikes.getsize() != ndates)		throw "Strikes array has the wrong size";
	if(callput.getsize() != ndates)		throw "Call/Put array has the wrong size";
	if(weights.getsize() != ndates)		throw "Weights array has the wrong size";
	if(periodId.getsize() != ndates)	throw "Period Identifier array has the wrong size";

	m_localCaps		=	localCaps;
	m_localFloors	=	localFloors;
	m_strikes		=	strikes;
	m_callput		=	callput;
	m_weights		=	weights;

	mNumberOfPeriods = couponPayDates.getsize();
	if(GlobalCap.getsize() != mNumberOfPeriods)		throw "Global Cap array has the wrong size	";
	if(GlobalFloor.getsize() != mNumberOfPeriods)	throw "Global Floor array has the wrong size ";
	if(FixedCoupon.getsize() != mNumberOfPeriods)	throw "Fixed Coupon array has the wrong size ";

	m_period.resize(mNumberOfPeriods);

	long nToday = MlEqDate::GetCurrentDate();		// ToDo - this is REALLY BAD. Can't have this once completed.
//	MlEqZeroCurveHandle hZeroCurve = hswap->GetA_ZeroCurve();	// get paying leg...

	for(int period = 0; period < mNumberOfPeriods; period++)
	{
		m_period[period].mCap			 = GlobalCap[period];
		m_period[period].mFloor			 = GlobalFloor[period];
		m_period[period].mFixedCoupon	 = FixedCoupon[period];
		m_period[period].mPayDate		 = couponPayDates[period];
		m_period[period].mDiscountFactor
			= 1.0;//hZeroCurve->GetDiscountFactor(nToday, couponPayDates[period]);	
		m_period[period].mFunding = 0.0;	// filled in later
	}

	std::set< int > idCopy;
	std::set< int >::iterator it = idCopy.begin() ;
	for(int i=0; i<ndates; i++)
		idCopy.insert( int(periodId[i]) );

	int idSize = idCopy.size() ;
	if(idSize != mNumberOfPeriods)	throw"Number of identifiers and periods mismatch";


	int checksize = 0;
	int idcheck = periodId[0];
	m_period[0].mStart = 1;

	for(int i=1; i<ndates; i++)
	{
		if( periodId[i] != idcheck )	// change period
		{
			if( periodId[i] == idcheck+1 )	
			{
				m_period[++checksize].mStart = i;

				if(checksize >= mNumberOfPeriods)
					throw"The fixing and coupon schedules sizes mismatch";

				idcheck++;
			}
			else
				throw" Two following identifiers values are inconsistant";
		}
		m_period[checksize].mEnd = i+1;
	}

	if(checksize != mNumberOfPeriods-1)	
		throw"The fixing and coupon schedules sizes mismatch, or an error in the identifiers occured";



	// find current period
	mCurrentPeriod = -1;	// forward start case

	for(int period = 0; period < mNumberOfPeriods; period++)
	{
		long start	=  fixingDates[ m_period[period].mStart-1 ] ;
		long end	=  fixingDates[ m_period[period].mEnd  -1 ] ;

		if( nToday >= start && nToday < end )
		{			
			mCurrentPeriod = period;

			int j = m_period[mCurrentPeriod].mStart - 1 ;	// readjust first index...
			while( nToday >= start )
				start = fixingDates[++j] ;

			m_period[mCurrentPeriod].mStart = j;
			break;			
		}	
	}

	int tmpCurrentPeriod = mCurrentPeriod ;
	mCurrentPeriod = MlEqMaths::Max( mCurrentPeriod, 0) ; // to avoid crashing in the forward start case...

	// swap part...

	for(int period = mCurrentPeriod; period < mNumberOfPeriods; period++)
	{
		long start = fixingDates[ m_period[period].mStart-1 ] ;
		if (period == mCurrentPeriod)
			start = MlEqMaths::Min( nToday, start ) ;
		long end = fixingDates[ m_period[period].mEnd-1 ] ;

//		m_period[period].mFunding = hswap->GetForwardValue(start, end);		
	}


	mCurrentPeriod = tmpCurrentPeriod ;

	init_payout(pastFixings, nToday);

	mCurrentPeriod = MlEqMaths::Max( mCurrentPeriod, 0) ;

}




void CallableCliquetSwapPricer::init_payout(CMatrix& pathArray, long nToday)
{

	isExpired = false ;	// still pay in the future if pay date > today...

	init_couponSum = 0.0;
	init_paidCoupons = 0.0 ;
	mCurrentCliquet = 0.0;

	int start, end;

	for (int iperiod = 0; iperiod <= mCurrentPeriod; iperiod++)
	{
		if(iperiod < mCurrentPeriod)
		{
			start	= m_period[iperiod].mStart;
			end		= m_period[iperiod].mEnd;
		}
		else	// current period
		{
			start	= MlEqMaths::Max( m_period[mCurrentPeriod-1].mEnd, 1) ;
			end		= m_period[mCurrentPeriod].mStart;
		}
	
		double funding	= m_period[iperiod].mFunding ;
		double cliquet = 0.0;

		for(int i=start; i<end; i++)
		{
			double tmpcliquet = pathArray[0][i] / pathArray[0][i-1] ;
			tmpcliquet = m_callput[i]*( tmpcliquet - m_strikes[i] );				
			tmpcliquet	= MlEqMaths::Min(MlEqMaths::Max(tmpcliquet,m_localFloors[i]),m_localCaps[i]);
	
			cliquet += tmpcliquet * m_weights[i] ;
		}

		if(iperiod < mCurrentPeriod)
			cliquet = m_period[iperiod].mFixedCoupon 
				+ MlEqMaths::Min(MlEqMaths::Max(cliquet,m_period[iperiod].mFloor),m_period[iperiod].mCap);
		
	
		if(iperiod == mCurrentPeriod) 
			mCurrentCliquet = cliquet;
		else
		{
			init_couponSum += cliquet ;		

			if(	init_couponSum > mCallableLevel )
			{				
				cliquet -= init_couponSum - mCallableLevel;
				init_couponSum = mCallableLevel;
				isExpired = true;
			}

			if( m_period[iperiod].mPayDate >= nToday)
				init_paidCoupons += cliquet * m_period[iperiod].mDiscountFactor + funding;
		}

		if( isExpired ) return;		
	}
}






void CallableCliquetSwapPricer::payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors, MlEqMonteCarlo& mc)
{

	values[0][0] = 0.0;

	if(isExpired)
	{
		values[0][0] =  init_paidCoupons ;
		fillValues(values,pathArray,mc);
		return;
	}

	double couponSum = init_couponSum;
	double paidCoupons = init_paidCoupons;	// for the swap...
	double df = 1.0;
	double funding = 0.0;
	double cliquet = mCurrentCliquet ;
		
	for (int iperiod = mCurrentPeriod; iperiod < mNumberOfPeriods; iperiod++)
	{
		int start	= m_period[iperiod].mStart;
		int end		= m_period[iperiod].mEnd;
		df			= m_period[iperiod].mDiscountFactor ;
		funding		= m_period[iperiod].mFunding ;

		for(int i=start; i<end; i++)
		{
			double tmpcliquet = pathArray[0][i] / pathArray[0][i-1] ;
			tmpcliquet = m_callput[i]*( tmpcliquet - m_strikes[i] );				
			tmpcliquet	= MlEqMaths::Min(MlEqMaths::Max(tmpcliquet,m_localFloors[i]),m_localCaps[i]);
	
			cliquet += tmpcliquet * m_weights[i] ;
		}

		cliquet = m_period[iperiod].mFixedCoupon 
				+ MlEqMaths::Min(MlEqMaths::Max(cliquet,m_period[iperiod].mFloor),m_period[iperiod].mCap);

		couponSum += cliquet ;

		if(	couponSum > mCallableLevel )
		{
			cliquet -= couponSum - mCallableLevel ;
			paidCoupons += cliquet * df + funding ; 
			values[0][0] =  paidCoupons ; 
			fillValues(values,pathArray,mc);
			return;
		}

		paidCoupons += cliquet * df + funding;
		cliquet = 0.0;
	}
	

	// the product hasn't been called at maturity
	paidCoupons += (mCallableLevel - couponSum) * df + funding ; 
	values[0][0] =  paidCoupons - df ; // redeem notional

		
	fillValues(values,pathArray,mc);
	return;
}






void CallableCliquetSwapPricer::setUp(CMatrix& value,MlEqMonteCarlo& mc)
{

	value.resize(1,2);
	m_results.resize(1,1);

	product::setUp(value,mc);
}












void AutoCallableSwapPricer::initialize(	const 	CVector	&rainbow_weights,	
											const 	CVector	&threshold_barriers,	
											const 	CVector	&coupon_up,	
											const 	CVector	&coupon_low,	
											const 	CVector	&call_barriers,	
											const	CVector	&rebates,
											const	GVector<long>& barrierDates,
											double	strike,			
											double	callput,
											double	gearing,
											double	fixed_coupon,				
											const	GVector< long >& couponPayDates,
											const	GVector< long >& couponResetDates,
											const	CVector& couponPayments,
											MlEqZeroCurveHandle hZeroCurve,
											MlEqConstDateHandle hDate,
											CMatrix& pastFixings)

{
	m_payoffDates	=	barrierDates;

	m_nBarrierDates = barrierDates.getsize();
	if(rebates.getsize() != m_nBarrierDates-1)		throw "Rebate array has the wrong size";
	if(call_barriers.getsize() != m_nBarrierDates-1)		throw "Barrier array has the wrong size";

	m_rebates		=	rebates;
	m_barriers		=	call_barriers;

	m_thresholds	=	threshold_barriers ;
	m_coupon_up		=	coupon_up ;
	m_coupon_low	=	coupon_low ;

	m_strike		=	strike;
	m_callput		=	callput;
	m_gearing		=	gearing;
	m_fixed_coupon	=	fixed_coupon;

	m_rainbow_weights = rainbow_weights ;

	if( !m_rainbow_weights.getsize() )
	{
		m_rainbow_weights.resize(1);
		m_rainbow_weights[0] = 1.0;
	}

	double sumWeights = 0.0;
	int nAsset	=	m_rainbow_weights.getsize() ;

	for(int iasset=0; iasset<nAsset; iasset++){		// rescale the whole thing...
		sumWeights += m_rainbow_weights[iasset] ;
	}
	CVector tmp = m_rainbow_weights ;
	for(int iasset=0; iasset<nAsset; iasset++){
		m_rainbow_weights[iasset] = tmp[nAsset-iasset-1] / sumWeights ;
	}


	long nToday = hDate->GetDate();

	// compute libor stuff

	m_libor_over_period.resize(m_nBarrierDates-1);

	int nPayDates = couponPayDates.getsize();
	if( couponResetDates.getsize() != nPayDates+1 ){
		throw "Reset and Pay schedule sizes mismatch";
	}
	if( couponPayments.getsize() != nPayDates ){
		throw "Coupons and Pay schedule sizes mismatch";
	}

	int ireset = 0;
	long prevBarrier = std::max( nToday , barrierDates[0] ) ;

	for( int ibarrier=1; ibarrier<m_nBarrierDates; ++ibarrier)
	{
		long nextBarrier = barrierDates[ibarrier];
		double libor_payment = 0.0;

		while(ireset < nPayDates && couponResetDates[ireset+1] <= nextBarrier )
		{
			if( couponResetDates[ireset+1] > prevBarrier ){
				libor_payment += couponPayments[ireset] * hZeroCurve->GetDiscountFactor( couponPayDates[ireset] );
			}		
			ireset++;
		}
		if( ireset < nPayDates && couponResetDates[ireset] < nextBarrier && couponResetDates[ireset+1] > nextBarrier ) 
		{			// accrued coupon in case this is called, added to the call rebate
		
			double tprec = hDate->GetYearFraction( couponResetDates[ireset] );
			double yFrac = hDate->GetYearFraction( couponResetDates[ireset+1] ) - tprec;
			double x = (hDate->GetYearFraction( nextBarrier) - tprec) / yFrac;
				
			m_rebates[ibarrier-1] += couponPayments[ireset] * x;
		}

		m_libor_over_period[ibarrier-1] = libor_payment;

		prevBarrier = nextBarrier ;
	}


	// check if the product expired...

	m_nStart = MlEqMaths::Max(pastFixings.cols(), 1);
	m_isExpired = false;
	std::vector<double> perfs( pastFixings.rows() );

	for(int ipast=1; ipast<m_nStart; ++ipast)
	{
		for(int ia=0; ia<nAsset; ++ia)
		{
			double perfa = pastFixings[ia][ipast] / pastFixings[ia][0];
			perfs[ia] = perfa;
		}

		double perf = 0.0;
		std::sort( perfs.begin(), perfs.end() ) ;

		for(int ia=0; ia<nAsset; ++ia)	{
			perf += perfs[ia] * m_rainbow_weights[ia] ;
		}

		if( perf > m_barriers[ipast-1] )
		{
			m_isExpired = true;
			break;
		}
	}
}



/*
void AutoCallableSwapPricer::payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MC& mc)
{
	values[0][0] = 0.0;

	if(m_isExpired)
	{
		values[0][0] =  m_rebates[m_nStart-1] ;
		fillValues(values,pathArray,mc);
		return;
	}

	double payoff = 0.0;;
		
	for (int ibarrier = m_nStart; ibarrier < m_nBarrierDates; ibarrier++)
	{
		payoff += m_libor_over_period[ibarrier-1]; 
		// Alex add a mapping and compute this properly in case stochastic rates are needed...

		double perf = pathArray[0][ibarrier] / pathArray[0][0];
		if( perf > m_barriers[ibarrier-1] )
		{
			payoff += m_rebates[ibarrier-1] * discountFactors[ibarrier];
			values[0][0] =  payoff ;
			fillValues(values,pathArray,mc);
			return;
		}
	}
	
	double final_payoff = pathArray[0][m_nBarrierDates-1] / pathArray[0][0] ;
	final_payoff = std::max(m_callput * (final_payoff - m_strike) , 0.0 );
	final_payoff = m_fixed_coupon + m_gearing * final_payoff ;

	values[0][0] =  payoff + final_payoff * discountFactors[m_nBarrierDates-1] ; 		
	fillValues(values,pathArray,mc);
	return;
}
*/


void AutoCallableSwapPricer::payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors, MlEqMonteCarlo& mc)
{
	values[0][0] = 0.0;

	if(m_isExpired)
	{
		values[0][0] =  m_rebates[m_nStart-1] ;
		fillValues(values,pathArray,mc);
		return;
	}

	double payoff = 0.0;;

	
	int nAsset = pathArray.rows();
	std::vector<double> perfs(nAsset) ; 
	double perf;
	

	for (int ibarrier = m_nStart; ibarrier < m_nBarrierDates; ibarrier++)
	{
		payoff += m_libor_over_period[ibarrier-1]; 
		// Alex add a mapping and compute this properly in case stochastic rates are needed...

		for(int ia=0; ia<nAsset; ++ia)
		{
			double perfa = pathArray[ia][ibarrier] / pathArray[ia][0];
			perfs[ia] = perfa;
		}

		perf = 0.0;
		std::sort( perfs.begin(), perfs.end() ) ;

		for(int ia=0; ia<nAsset; ++ia)	{
			perf += perfs[ia] * m_rainbow_weights[ia] ;
		}


		double coupon_t = ( perf < m_thresholds[ibarrier-1] )? m_coupon_low[ibarrier-1] : m_coupon_up[ibarrier-1] ;
		payoff += coupon_t * discountFactors[ibarrier] ;

		if( perf > m_barriers[ibarrier-1] )
		{
			payoff += m_rebates[ibarrier-1] * discountFactors[ibarrier];
			values[0][0] =  payoff ;
			fillValues(values,pathArray,mc);
			return;
		}
	}
	
	double final_payoff = std::max(m_callput * (perf - m_strike) , 0.0 );
	final_payoff = m_fixed_coupon + m_gearing * final_payoff ;

	values[0][0] =  payoff + final_payoff * discountFactors[m_nBarrierDates-1] ; 		
	fillValues(values,pathArray,mc);
	return;
}



void AutoCallableSwapPricer::setUp(CMatrix& value,MlEqMonteCarlo& mc)
{
	value.resize(1,2);
	m_results.resize(1,1);

	product::setUp(value,mc);
}














/*
void CallableCliquetSwapPricer::init_payout(CMatrix& pathArray, long nToday)
{

	isExpired = false ;	// still pay in the future if pay date > today...

	init_couponSum = 0.0;
	init_paidCoupons = 0.0 ;
	mCurrentCliquet = 0.0;

	int start, end;

	for (int iperiod = 0; iperiod <= mCurrentPeriod; iperiod++)
	{
		if(iperiod < mCurrentPeriod)
		{
			start	= m_period[iperiod].mStart;
			end		= m_period[iperiod].mEnd;
		}
		else	// current period
		{
			start	= MlEqMaths::Max( m_period[mCurrentPeriod-1].mEnd, 1) ;
			end		= m_period[mCurrentPeriod].mStart;
		}
	
		double funding	= m_period[iperiod].mFunding ;
		double cliquet = 0.0;

		for(int i=start; i<end; i++)
		{
			double tmpcliquet = pathArray[0][i] / pathArray[0][i-1] ;
			tmpcliquet = m_callput[i]*( tmpcliquet - m_strikes[i] );				
			tmpcliquet	= MlEqMaths::Min(MlEqMaths::Max(tmpcliquet,m_localFloors[i]),m_localCaps[i]);
	
			cliquet += tmpcliquet * m_weights[i] ;
		}

		if(iperiod < mCurrentPeriod)
			cliquet = m_period[iperiod].mFixedCoupon 
				+ MlEqMaths::Min(MlEqMaths::Max(cliquet,m_period[iperiod].mFloor),m_period[iperiod].mCap);
		
	
		if(iperiod == mCurrentPeriod) 
			mCurrentCliquet = cliquet;
		else
		{
			init_couponSum += cliquet ;		

			if(	init_couponSum > mCallableLevel )
			{				
				cliquet -= init_couponSum - mCallableLevel;
				init_couponSum = mCallableLevel;
				isExpired = true;
			}

			if( m_period[iperiod].mPayDate >= nToday)
				init_paidCoupons += cliquet * m_period[iperiod].mDiscountFactor + funding;
		}

		if( isExpired ) return;		
	}
}

*/


