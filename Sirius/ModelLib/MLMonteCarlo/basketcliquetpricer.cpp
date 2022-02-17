
#include "stdafx.h"
#include "MlEqObjects.h"
#include "BasketCliquetPricer.h"
#include "MonteCarlo.h"
#include "MlEqPde.h"


/****************************************************************
**	Class  : CappedFlooredCliquetPricer 
**	Routine: initialize
**	Returns: nothing
**	Action : 
**           
****************************************************************/


void BasketCliquetPricer::initialize(	
						CVector	&localCaps,//[idate]
						CVector	&localFloors,//[idate]
						CVector	&strikes,//[idate]
						CVector	&callput,//[idate]	
						double GlobalCap,
						double GlobalFloor,
						double GlobalRedemption,
						double GlobalGearing,
						double Notional,
						GVector< long >& payoffDates,
						CVector rainbowWeights
						)
{

	int ndates = payoffDates.getsize();

	m_localCaps			= localCaps;
	m_localFloors		= localFloors;
	m_strikes			= strikes;
	m_callput			= callput;	

	m_GlobalCap			= GlobalCap;
	m_GlobalFloor		= GlobalFloor;
	m_GlobalRedemption	= GlobalRedemption;
	m_GlobalGearing		= GlobalGearing;
	m_Notional			= Notional;

	m_payoffDates		= payoffDates;


	if( m_localCaps.getsize() != ndates-1)	throw "wrong local caps size";
	if( m_localFloors.getsize() != ndates-1)	throw "wrong local floors size";
	if( m_strikes.getsize() != ndates-1)	throw "wrong strikes size";
	if( m_callput.getsize() != ndates-1)	throw "wrong call put size";



	m_rainbowWeights	= rainbowWeights;

	double sumOfWeights = 0.0;
	for ( int i = 0 ; i < m_rainbowWeights.getsize(); i++ ){
		sumOfWeights += m_rainbowWeights[i];
	}
	if( sumOfWeights < 1e-6 )
		throw "You must enter some weights (from worst to best)";

	for ( int i = 0 ; i < m_rainbowWeights.getsize(); i++ ){	// renormalise weights
		m_rainbowWeights[i] /= sumOfWeights ;
	}

}


void BasketCliquetPricer::payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc)
{

	int ndates = pathArray.cols();
	int nasset = pathArray.rows();

	double xRet, x;

	std::vector<double> cliquetPerf;
		
	for( int iasset = 0; iasset<nasset; iasset++ )
	{
		x = 0.;
		for (int idate = 1 ; idate < ndates; idate++ )
		{	
			xRet = pathArray[iasset][idate] / pathArray[iasset][idate-1]; 
			xRet = m_callput[idate-1]*(xRet-m_strikes[idate-1]);				
		
			x += std::min( std::max(xRet,m_localFloors[idate-1]), m_localCaps[idate-1] );
		}
		cliquetPerf.push_back( x );
	}

	std::sort( cliquetPerf.begin(), cliquetPerf.end() );

	double payout = 0.;
	for( int iworst = 0; iworst<nasset; iworst++ )	{
		payout += cliquetPerf[iworst] * m_rainbowWeights[iworst] ;
	}
		

	payout = m_GlobalGearing * payout + m_GlobalRedemption;		
	payout = std::min( std::max(payout, m_GlobalFloor), m_GlobalCap );

	values[0][0] = payout * m_Notional * discountFactors[ndates-1];		
	fillValues(values,pathArray,mc);

}



void BasketCliquetPricer::setUp(CMatrix& value,MlEqMonteCarlo& mc)
{

	value.resize(1,2);
	m_results.resize(1,1);

	product::setUp(value,mc);
}





//////////////////////////////////////
//									//
//			Asian Pricer			//
//									//
//////////////////////////////////////




void AsianPricer::initialize(	double Strike,
								double Cp,
								const GVector< long >& payoffDates,
								const CVector& fixings )
{
	m_payoffDates		= payoffDates;
	m_strike			= Strike;
	m_cp				= Cp;	
	m_nAvgDates			= m_payoffDates.getsize();

	m_pastSum = 0.0;
	m_nStart = 0;

	int npast = MlEqMaths::Min( int(fixings.getsize()), m_nAvgDates);
	for(int i=0; i<npast; ++i)
	{
		if( fixings[i] > 0.0 )
		{
			m_pastSum += fixings[i];
			m_nStart++;
		}
		else break;
	}

}

/*
void AsianPricer::payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MC& mc)
{
	double sum = m_pastSum;
	for (int idate = m_nStart ; idate < m_nAvgDates; ++idate ){	
		sum += pathArray[0][idate] ; 
	}	
	sum /= double(m_nAvgDates) ;

	double payout = std::max(m_cp * ( sum - m_strike), 0.0 );

	values[0][0] = payout * discountFactors[m_nAvgDates-1];		
	fillValues(values,pathArray,mc);
}



void AsianPricer::setUp(CMatrix& value,MC& mc)
{
	value.resize(1,2);
	m_results.resize(1,1);

	product::setUp(value,mc);
}
*/


void AsianPricer::payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc)
{
	double df = discountFactors[m_nAvgDates-1];	
	int ngreeks = pathArray.rows();

	for(int k=0; k<ngreeks; ++k)
	{
		double payoff = payout( pathArray[k] );
		values[k][0] = payoff * df;
	}

	fillValues(values,pathArray,mc);
}



double AsianPricer::payout(const CVector& pathArray)
{
	double sum = m_pastSum;

	for (int idate = m_nStart ; idate < m_nAvgDates; ++idate )
	{	
		sum += pathArray[idate] ; 
	}	
	sum /= double(m_nAvgDates) ;


	double payout = std::max(m_cp * ( sum - m_strike), 0.0 );
	return payout ;
}


void AsianPricer::setUp(CMatrix& value,MlEqMonteCarlo& mc)
{
	value.resize(4,2);
	m_results.resize(4,1);

	product::setUp(value,mc);
}


