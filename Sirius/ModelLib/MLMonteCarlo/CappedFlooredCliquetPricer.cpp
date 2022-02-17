//	cappedflooredcliquet.cpp : Implementation of CCappedFlooredCliquet
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqObjects.h"
#include "CappedFlooredCliquetPricer.h"
#include "MonteCarlo.h"
#include "MlEqPde.h"


/****************************************************************
**	Class  : CappedFlooredCliquetPricer 
**	Routine: initialize
**	Returns: nothing
**	Action : 
**           
****************************************************************/


void CappedFlooredCliquetPricer::initialize(	
						CVector	&localCaps,//[idate]
						CVector	&localFloors,//[idate]
						CVector	&strikes,//[idate]
						CVector	&callput,//[idate]	
						CVector	&weights,//[idate]
						double GlobalCap,
						double GlobalFloor,
						double GlobalRedemption,
						double GlobalGearing,
						double Notional,
						GVector< long >& payoffDates,
						CVector basketWeights,
						bool rebalancingBasket
						)
{

	m_localCaps			= localCaps;
	m_localFloors		= localFloors;
	m_strikes			= strikes;
	m_callput			= callput;	
	m_weights			= weights;

	m_GlobalCap			= GlobalCap;
	m_GlobalFloor		= GlobalFloor;
	m_GlobalRedemption	= GlobalRedemption;
	m_GlobalGearing		= GlobalGearing;
	m_Notional			= Notional;

	m_payoffDates		= payoffDates;
	m_basketWeights		= basketWeights;

/*	m_sumOfBasketWeights = 0.0;
	for ( int i = 0 ; i < m_basketWeights.getsize(); i++ ){
		m_sumOfBasketWeights += basketWeights[i];
	}
*/
	m_LockinPayoff = false ;

	m_rebalancingBasket = rebalancingBasket ;

}


void CappedFlooredCliquetPricer::initialize(	
						CVector	&localCaps,//[idate]
						CVector	&localFloors,//[idate]
						CVector	&strikes,//[idate]
						CVector	&callput,//[idate]	
						CVector	&weights,//[idate]
						double GlobalCap,
						double GlobalFloor,
						double GlobalRedemption,
						double GlobalGearing,
						double Notional,
						GVector< long >& payoffDates,
						CVector basketWeights,
						const CMatrix& LockinSchedule,
						bool rebalancingBasket
						)
{

	m_localCaps			= localCaps;
	m_localFloors		= localFloors;
	m_strikes			= strikes;
	m_callput			= callput;	
	m_weights			= weights;

	m_GlobalCap			= GlobalCap;
	m_GlobalFloor		= GlobalFloor;
	m_GlobalRedemption	= GlobalRedemption;
	m_GlobalGearing		= GlobalGearing;
	m_Notional			= Notional;

	m_payoffDates		= payoffDates;
	m_basketWeights		= basketWeights;

/*	m_sumOfBasketWeights = 0.0;
	for ( int i = 0 ; i < m_basketWeights.getsize(); i++ ){
		m_sumOfBasketWeights += basketWeights[i];
	}
*/
	m_LockinSchedule = LockinSchedule ;
	m_LockinPayoff = true;
	m_LockinCapped = (m_LockinSchedule.cols() == 2) ;

	m_rebalancingBasket = rebalancingBasket ;
}


/****************************************************************
**	Class  : CappedFlooredCliquetPricer 
**	Routine: evaluate
**	Returns: nothing
**	Action : calculates payout
**           
****************************************************************/



void CappedFlooredCliquetPricer::evaluate(CMatrix& Results,CForwardSkewMC& mc)
{	


	CVector					values,results;			/* index	SUM		lc/lf	GC/GF
														0		 N		 N		 N
														1		 Y		 N		 N
														2		 N		 Y		 N
														3		 Y		 Y		 N
														4		 N		 N		 Y
														5		 Y		 N		 Y
														6		 N		 Y		 Y
														7		 Y		 Y		 Y
													*/

	int ndates = m_weights.getsize();

	if ( m_localCaps.getsize()!= ndates || m_localFloors.getsize() != ndates || m_strikes.getsize()!=ndates || m_callput.getsize() != ndates || ndates > mc.m_nDates)
	{
		throw("error in setting up localCapFloor inputData");
	}
	
	
	
	values.resize(8);
	results.resize(8);
	int ipath, idate, i;
	for ( i = 0; i < m_weights.getsize(); i++ )
	{
		m_weights[i] *= m_GlobalGearing;
	}
	
	
	
	double xRet, x;
	
	
	
	for ( ipath = 0 ; ipath < mc.m_nPaths; ipath++ )
	{	
		
		values[1] = 0.0;
		values[3] = 0.0;
		values[5] = 0.0;
		values[7] = 0.0;
		
		values[0] = 1.0;
		values[2] = 1.0;
		values[4] = 1.0;
		values[6] = 1.0;
		
		for ( idate = 1 ; idate <= ndates; idate++ )
		{	
			xRet = calculateReturn(mc,idate,ipath);

			xRet = m_callput[idate-1]*(xRet-m_strikes[idate-1]);
			
			
			values[0] *= xRet;
			values[1] += xRet;
				
			x	= MlEqMaths::Min(MlEqMaths::Max(xRet,m_localFloors[idate-1]),m_localCaps[idate-1]);
			
			values[2] *= x;
			values[3] += x;
			
			for ( i = 0; i < 4; i++ )
			{
				values[i] *= m_weights[idate-1];
			}
		}	
		
		for ( i = 0; i < 4; i++ )
		{
			values[i] += m_GlobalRedemption;
		}
		
		values[4]	= MlEqMaths::Min(MlEqMaths::Max(values[0], m_GlobalFloor), m_GlobalCap);
		values[5]	= MlEqMaths::Min(MlEqMaths::Max(values[1], m_GlobalFloor), m_GlobalCap);
		values[6]	= MlEqMaths::Min(MlEqMaths::Max(values[2], m_GlobalFloor), m_GlobalCap);
		values[7]	= MlEqMaths::Min(MlEqMaths::Max(values[3], m_GlobalFloor), m_GlobalCap);
			
		for ( i = 0; i < 8; i++ )
		{
			results[i] += values[i];
			results[i] *= mc.GetDiscount(ipath,ndates)/mc.m_nPaths;

		}	
	}	
	

	Results.resize(1,results.getsize());
	Results[0] = results;

	return;
	
}	

/****************************************************************
**	Class  : CappedFlooredCliquetPricer 
**	Routine: payout
**	Returns: nothing
**	Action : calculates payout
**           
****************************************************************/


double CappedFlooredCliquetPricer::calculateReturn(CMatrix& pathArray,int idate)
{
	double xRet;

	if ( idate < 1 ){
		throw("return cannot be calculated");
	}

	if ( pathArray.rows() == 1 )
	{
		int iasset = 0 ;
		xRet = pathArray[iasset][idate] / pathArray[iasset][idate-1];
		return xRet;
	}


	if( !m_rebalancingBasket )
	{
		double xBasket = 0.0;
		double xBasketPrev = .0;
		for ( int iasset = 0 ; iasset < pathArray.rows(); iasset++ )
		{
			double wght = m_basketWeights[iasset];
	
			xBasket		+= wght * pathArray[iasset][idate] ;
			xBasketPrev	+= wght * pathArray[iasset][idate-1]  ;
		}

		xRet = xBasket/xBasketPrev;
		return xRet;
	}


	xRet = 0.0;
	for ( int iasset = 0 ; iasset < pathArray.rows(); iasset++ ){
		xRet += m_basketWeights[iasset]*pathArray[iasset][idate] / pathArray[iasset][idate-1]; 
	}

	return xRet;

}

/****************************************************************
**	Class  : CappedFlooredCliquetPricer 
**	Routine: payout
**	Returns: nothing
**	Action : calculates payout
**           
****************************************************************/

double CappedFlooredCliquetPricer::calculateReturn(CForwardSkewMC& mc,int idate,int ipath)
{
	double xRet;

	if ( idate < 1 ){
		throw("return cannot be calculated");
	}

	if ( mc.m_nAssets == 1 )
	{
		xRet = mc.GetPathValue(ipath,idate) / mc.GetPathValue(ipath,idate-1);
		return xRet;
	}

	if( !m_rebalancingBasket )
	{
		double xBasket = 0.0;
		double xBasketPrev = .0;
		for ( int iasset = 0 ; iasset < mc.m_nAssets; iasset++ )
		{
			double wght = m_basketWeights[iasset];

			xBasket		+= wght * mc.GetPathValue(ipath,idate,iasset)  ;
			xBasketPrev	+= wght * mc.GetPathValue(ipath,idate-1,iasset)  ;
		}

		xRet = xBasket/xBasketPrev;
		return xRet;
	}

	xRet = 0.0;
	for ( int iasset = 0 ; iasset < mc.m_nAssets; iasset++ ){
		xRet += m_basketWeights[iasset]*mc.GetPathValue(ipath,idate,iasset) / mc.GetPathValue(ipath,idate-1,iasset);
	}

	return xRet;

}


/****************************************************************
**	Class  : CappedFlooredCliquetPricer 
**	Routine: payout
**	Returns: nothing
**	Action : calculates payout
**           
****************************************************************/

void CappedFlooredCliquetPricer::payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc)
{
	if(m_LockinPayoff )		payoutLockin(values, pathArray, discountFactors, mc);
	else					payoutClassic(values, pathArray, discountFactors, mc);
}



//void CappedFlooredCliquetPricer::payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MC& mc)
void CappedFlooredCliquetPricer::payoutClassic(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc)
{

		/* index	SUM		lc/lf	GC/GF
			0		 N		 N		 N
			1		 Y		 N		 N
			2		 N		 Y		 N
			3		 Y		 Y		 N
			4		 N		 N		 Y
			5		 Y		 N		 Y
			6		 N		 Y		 Y
			7		 Y		 Y		 Y
										*/

	int n = m_weights.getsize();

	if ( m_localCaps.getsize()!= n || m_localFloors.getsize() != n || m_strikes.getsize()!=n || m_callput.getsize() != n )
	{
		throw("error in setting up localCapFloor inputData");
	}

	int ndates = pathArray.cols();
		


	int idate, i;

	
	
	double xRet, x;
		
	values[1][0] = 0.0;
	values[3][0] = 0.0;
	values[5][0] = 0.0;
	values[7][0] = 0.0;
		
	values[0][0] = 1.0;
	values[2][0] = 1.0;
	values[4][0] = 1.0;
	values[6][0] = 1.0;
		
	int iasset = 0;
	for ( idate = 1 ; idate < ndates; idate++ )
	{	
		xRet = calculateReturn(pathArray,idate);
		xRet = m_callput[idate-1]*(xRet-m_strikes[idate-1]);
			
			
		values[0][0] *= xRet;
		values[1][0] += xRet;
				
		x	= MlEqMaths::Min(MlEqMaths::Max(xRet,m_localFloors[idate-1]),m_localCaps[idate-1]);
			
		values[2][0] *= x*m_weights[idate-1];;
		values[3][0] += x*m_weights[idate-1];;
			
//		for ( i = 0; i < 4; i++ )
//		{
//			values[i][0] *= m_weights[idate-1];
//		}
	}	
		
	for ( i = 0; i < 4; i++ )
	{
		values[i][0] = m_GlobalGearing*values[i][0] + m_GlobalRedemption;
	}
		
	values[4][0]	= MlEqMaths::Min(MlEqMaths::Max(values[0][0], m_GlobalFloor), m_GlobalCap);
	values[5][0]	= MlEqMaths::Min(MlEqMaths::Max(values[1][0], m_GlobalFloor), m_GlobalCap);
	values[6][0]	= MlEqMaths::Min(MlEqMaths::Max(values[2][0], m_GlobalFloor), m_GlobalCap);
	values[7][0]	= MlEqMaths::Min(MlEqMaths::Max(values[3][0], m_GlobalFloor), m_GlobalCap);
			
	for ( i = 0; i < values.rows(); i++ )
	{
		values[i][0] *= m_Notional*discountFactors[ndates-1];
	}
		
	fillValues(values,pathArray,mc);


	return;

}

void CappedFlooredCliquetPricer::payoutLockin(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc)
{
	int ndates = pathArray.cols();
	int nlevels = m_LockinSchedule.rows(); 
	double globalFloor = m_GlobalFloor ;
	double payoff = 0.0;

	for (int idate = 1 ; idate < ndates; idate++ )
	{	
		double cliquet = calculateReturn(pathArray,idate) - m_strikes[idate-1];	
		cliquet = m_callput[idate-1] * cliquet ;
		cliquet	= std::min( std::max(cliquet,m_localFloors[idate-1]), m_localCaps[idate-1] );
			
		payoff += cliquet * m_weights[idate-1];

		globalFloor = std::max( globalFloor, payoff );
	}

	if( m_LockinCapped )
	{
		double locklevel = globalFloor ;
		globalFloor = m_GlobalFloor;
		int nextLevel = 0;
		while( nextLevel < nlevels && locklevel >= m_LockinSchedule[nextLevel][0] ){
			globalFloor = m_LockinSchedule[nextLevel++][1] ;
		}
	}
	else if( globalFloor < m_LockinSchedule[0][0]){
		globalFloor = m_GlobalFloor ;
	}
	
	globalFloor = std::max( globalFloor, m_GlobalFloor );		

	payoff = m_GlobalGearing * payoff + m_GlobalRedemption;		
	payoff = std::min( std::max(payoff, globalFloor), m_GlobalCap );		

	values[7][0] = m_Notional * discountFactors[ndates-1] * payoff;	// cliquet storage convention...
	fillValues(values,pathArray,mc);

	return;
}




/****************************************************************
**	Class  : CappedFlooredCliquetPricer 
**	Routine: setUp
**	Returns: nothing
**	Action : calculates payout
**           
****************************************************************/


void CappedFlooredCliquetPricer::setUp(CMatrix& value,MlEqMonteCarlo& mc)
{

	value.resize(8,2);
	m_results.resize(8,2);

	product::setUp(value,mc);
}





/****************************************************************
**	Class  : Himalaya 
**	Routine: payout
**	Returns: nothing
**	Action : calculates payout
**           
****************************************************************/


void Himalaya::payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc)
{
//  pathArray[iasset][idate]
	int ndates = pathArray.cols();		
	int idate,iasset;
	int currentWinner=0;
				
	CVector vals(m_nAssets),wals(m_nAssets);
	CVector payoff(m_nAssets);
	double payout = 0.0;

	for ( iasset = 0 ; iasset < m_nAssets; iasset++ ){
		m_exclude[iasset] = false;
	}

	double n = 0 ;
	for ( idate = 0 ; idate < ndates; idate++ )
	{	
		n++;
		for ( iasset = 0 ; iasset < m_nAssets; iasset++ ){
			vals[iasset] += pathArray[iasset][idate];
		}

		if ( callHimalaya(payout,idate,vals,discountFactors,n) ){
			values[0][0] = m_Notional*payout;
			return;
		}

		if ( m_isHimalayaDate[idate] )
		{

			for ( iasset = 0 ; iasset < m_nAssets; iasset++ ){
					wals[iasset] = vals[iasset]/(n*m_spots[iasset]);
			}

			payout = selectHimalayaWinner(idate,wals,currentWinner);
			payoff[currentWinner] = payout;

			if ( !m_asianFromStart )
			{
				for ( iasset = 0 ; iasset < m_nAssets; iasset++ ){
					vals[iasset] = 0.0;
				}
				n = 0;
			}

		}
	}

	if ( !m_isHimalaya ){
		for ( int i = 0 ; i < m_nAssets; i++ ){
			payoff[i] = vals[i]/(n*m_spots[iasset]);
		}
	}

	callRainbow(payout,payoff);

	values[0][0] = MlEqMaths::Max(m_callPut*(payout-m_strike),0.0);
 	values[0][0] *= m_Notional*discountFactors[ndates-1];
	;
	fillValues(values,pathArray,mc);

}

/****************************************************************
**	Class  : Himalaya 
**	Routine: callRainbow
**	Returns: double
**	Action : selects Himalaya winner
**           
****************************************************************/
/*
void Himalaya::callRainbow(double & res,CVector& returns)

	payoff = MlEqMaths::Max(m_callPut*(payoff-m_strike),0.0);

	if ( idate == ndates ){
		values[0][0] = m_Notional*payoff*discountFactors[ndates-1];	
	}else{
		//  payoff was called in this case
			values[0][0] = m_Notional*payoff;	
	}

	fillValues(values,pathArray,mc);

	return;
}
*/

/****************************************************************
**	Class  : Himalaya 
**	Routine: selectHimalayaWinner
**	Returns: double
**	Action : selects Himalaya winner
**           
****************************************************************/


double Himalaya::selectHimalayaWinner(int idate,CVector& values,int& currentWinner)
{

	double currMax = -1;

	for ( int iasset = 0 ; iasset < m_nAssets; iasset++ )
	{
		if ( m_exclude[iasset] ){
			continue;
		}

		if ( values[iasset] > currMax )
		{
			currentWinner = iasset;
			currMax = values[iasset];
		}
	}

	if ( m_exclude[currentWinner] == true ){
		throw(" himalay winner was selected already");
	}
	else{
		m_exclude[currentWinner] = true; 
	}

	return currMax;
}



/****************************************************************
**	Class  : Himalaya 
**	Routine: callHimalaya
**	Returns: bool
**	Action : check for callable feature
**           
****************************************************************/


bool Himalaya::callHimalaya(double& payoff,int idate,CVector& values,CVector& discountFactors,double divisor)
{

	int ndates = m_payoffDates.getsize();

	if ( !m_isCallable[idate] ){
		return false;
	}

	int iasset;
	for ( iasset = 0 ; iasset < m_nAssets; iasset++ )
	{
		if ( values[iasset]/divisor < m_callLevels[idate] ){
			break;
		}
	}

	if ( iasset < m_nAssets ){
		return false;
	}
	else
	{
		payoff = m_rebate[idate-1];
		if ( m_delayRebate ){
			payoff *= discountFactors[ndates-1];
		}else{
			payoff *= discountFactors[idate];
		}

		return true;
	}

}




/****************************************************************
**	Class  : Himalaya 
**	Routine: callHimalaya
**	Returns: void
**	Action : setUp
**           
****************************************************************/


void Himalaya::setUp(CMatrix& value,MlEqMonteCarlo& mc)
{
	value.resize(1,2);
	m_results.resize(1,2);

	if ( mc.m_nAssets != m_nAssets ){
		throw("number of assets in MonteCarlo is not consitent with Himalaya set up");
	}

	product::setUp(value,mc);

}

/****************************************************************
**	Class  : Himalaya 
**	Routine: initialize
**	Returns: void
**	Action : initializes
**           
****************************************************************/

void Himalaya::initialize(int nAssets, CVector& spot,GVector<long> payoffDates, GVector<long> isHimalayaDate, CVector rainbowWeights, double notional, double strike, double callPut, CVector callLevels, CVector rebate, bool delayRebate, bool isAsianFromStart)
{

	int ndates = payoffDates.getsize();
	GVector<bool> isCallable(ndates);

	if ( isHimalayaDate.getsize() == 0 ){
		isHimalayaDate.resize(ndates);
	}


	if ( callLevels.getsize() == 0 )
	{
		callLevels.resize(ndates);
		for ( int i = 0 ; i < ndates; i++ ){
			isCallable[i] = false;
		}
	}
	else
	{
		if ( callLevels.getsize() != ndates ){
			throw("arraysize of callLevel array is incorrect!");
		}

		for ( int i = 0 ; i < ndates; i++ )
		{
			if ( fabs(callLevels[i]) > 1e-3 ){
				isCallable[i] = true;
			}
			else{
				isCallable[i] = false;
			}
		}
	}

	m_exclude.resize(nAssets);
	m_spots = spot;

	if ( isCallable.getsize() != ndates ){
		throw("incorrect arraysize encountered in isCallable array");
	}

	m_isCallable = isCallable;
	
	if ( callLevels.getsize() != ndates ){
		throw("incorrect arraysize encountered in callLevels array");
	}
	m_callLevels = callLevels;

	if ( isHimalayaDate.getsize() != ndates ){
		throw(" isHimalayadate array has incorrect size");
	}

	m_nAssets = nAssets;
	m_isHimalayaDate = isHimalayaDate;


	int icount = 0 ;
	for ( int i = 0 ; i < ndates; i++ ){
		if ( m_isHimalayaDate[i] ){
			icount++;
		}
	}

	if ( icount == 0 ){
		m_isHimalaya = false;
		// this is a regular rainbow option
	}
	else{
		if ( icount != m_nAssets ){
			throw("number of Himalaya points do not coincide with number of assets");
		}
		m_isHimalaya = true;
	}

	m_asianFromStart	= isAsianFromStart;
	m_strike			= strike;
	m_Notional			= notional;
	m_callPut			= callPut;
	m_rebate			= rebate;
	m_delayRebate		= delayRebate;
	m_rainbowWeights	= rainbowWeights;
	m_payoffDates		= payoffDates;

	if ( m_rainbowWeights.getsize() != 0 ){
		if ( m_rainbowWeights.getsize() != m_nAssets ){
			throw("incorrect number of rainbow weights have been entered");
		}
	}
}


/****************************************************************
**	Class  : Himalaya 
**	Routine: callRainbow
**	Returns: void
**	Action : performs rainbow operation
**           
****************************************************************/

void Himalaya::callRainbow(double & res,CVector& returns)
{
	if ( m_rainbowWeights.getsize() == 0 )
	{
		res = 0.0;
		for ( int i = 0 ; i < returns.getsize(); i++ ){
			res += returns[i];
		}

		res /= returns.getsize();

	}
	else
	{

		MlEqMaths::Hpsort(returns);

		res = 0.0;
		for ( int i = 0 ; i < m_rainbowWeights.getsize(); i++ ){
			res += returns[i]*m_rainbowWeights[i];
		}

	}
}
	
	
	
/****************************************************************
**	Class  : CAsianMax 
**	Routine: callRainbow
**	Returns: void
**	Action : performs rainbow operation
**           
****************************************************************/
	
	
void CAsianMax::payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc)
{	
	
	int nassets = pathArray.rows();
	int ndates  = pathArray.cols();
	
	CVector MaxAsianReturns(nassets);

	double currentAverage;
	CVector	m_spotFixings;

	for ( int iasset = 0 ; iasset < nassets; iasset++ )
	{
		currentAverage = 0.0;
		for ( int idate = 0; idate < ndates; idate++ )
		{
			currentAverage = 1.0/(double)(idate+1)*(idate*currentAverage + pathArray[iasset][idate]/m_spotFixings[iasset]);
			MaxAsianReturns[iasset] = MlEqMaths::Max(MaxAsianReturns[iasset],currentAverage);					
		}
	}	

	double globalMax = MaxAsianReturns[0];
	for ( int iasset = 1; iasset < nassets; iasset++ ){
		globalMax = MlEqMaths::Max(MaxAsianReturns[1],globalMax);
	}

	for ( int istrike = 0 ; istrike < m_strike.getsize(); istrike++ ){
		values[istrike][0] = MlEqMaths::Max(m_callPut*(globalMax-m_strike[istrike]),0.0);
	}
		
}



/****************************************************************
**	Class  : CAsianMax 
**	Routine: callRainbow
**	Returns: void
**	Action : performs rainbow operation
**           
****************************************************************/
	
	
void CAsianMax::initialize(double callPut,CVector& strikes, GVector<long> Dates,CVector& spotFixings)
{	
	m_payoffDates	=	Dates;
	m_callPut		=	callPut;
	m_spotFixings	=	spotFixings;
	m_strike		=	strikes;
}



/****************************************************************
**	Class  : CappedFlooredCliquetPricer 
**	Routine: setUp
**	Returns: nothing
**	Action : calculates payout
**           
****************************************************************/


void CAsianMax::setUp(CMatrix& value,MlEqMonteCarlo& mc)
{
	value.resize(m_strike.getsize(),2);
	m_results.resize(m_strike.getsize(),2);
	product::setUp(value,mc);
}



