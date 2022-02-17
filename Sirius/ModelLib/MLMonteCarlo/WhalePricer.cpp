

#include "stdafx.h"
#include "MlEqObjects.h"
#include "WhalePricer.h"
#include "MonteCarlo.h"
#include "MlEqPde.h"



void WhalePricer::initialize(	CVector	&weights,
								double strike,
								double cp,
								long nToday,
								const GVector<long>&	payoffDates,
								const CMatrix& fixings )
{
	m_payoffDates		= payoffDates;
	m_strike			= strike ;
	m_cp				= cp ;

	m_nDates = m_payoffDates.getsize(); 

	if(weights.getsize() != m_nDates-1)		throw "Weights array has the wrong size";
	m_weights		=	weights;


	m_runningAverage = 0.0;
	m_sumWeights = 0.0;
	m_start = 1;

	for(int i=1; i<m_nDates; i++)
	{
		if( payoffDates[i] < nToday )
		{
			m_start++;
			m_runningAverage += fixings[0][i] * m_weights[i-1];
		}
		m_sumWeights += m_weights[i-1];
	}

}


void WhalePricer::payout(CMatrix& values, CMatrix& pathArray, CVector& discountFactors, MlEqMonteCarlo& mc)
{			
	double	average = m_runningAverage;

	for(int ndate=m_start; ndate<m_nDates; ndate++)	{
		average += pathArray[0][ndate] * m_weights[ndate-1]  ;
	}
	average /= m_sumWeights * pathArray[0][0];
	average = 1./average ;

	double payoff = std::max( m_cp*(average - m_strike), 0.0);
	
	values[0][0] = payoff;
}


void WhalePricer::setUp(CMatrix& value, MlEqMonteCarlo& mc)
{
	value.resize(1, 2);
	m_results.resize(1, 2);
	product::setUp(value,mc);
}

