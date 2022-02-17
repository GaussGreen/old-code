

#include "stdafx.h"
#include "MlEqObjects.h"
#include "MonteCarlo.h"

#include "IguanaProduct.h"

void IguanaPricer::initialize(						
						double Strike,
						double Barrier,
						GVector< long >& payoffDates,
						CVector& spotFixings
						)
{
	m_Strike	=	Strike ;
	m_Barrier	=	Barrier ;
	m_payoffDates	= payoffDates;

	m_scheduleSize = payoffDates.getsize();
	m_start = 1;

	long today = MlEqDate(MlEqDate::GetCurrentDate()).GetDate();
	m_pastMin = 1e9;
	m_pastMax = -1.0;

	for(int i=1; i<m_scheduleSize; ++i)
	{
		if(payoffDates[i] < today)
		{
			m_start++;
			if(m_start > payoffDates.getsize())
				throw "some fixing values are missing, you dick wad";

			double fixedSpot = spotFixings[i] / spotFixings[0] ;
			m_pastMin = MlEqMaths::Min( m_pastMin, fixedSpot );
			m_pastMax = MlEqMaths::Max( m_pastMax, fixedSpot );
		}
		else break;
	}
}




void IguanaPricer::payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc)
{

	double mmax = m_pastMax;
	double mmin = m_pastMin;
	double refSpot = pathArray[0][0] ;
		
	values[0][0] = 0.0;

	for (int idate = m_start ; idate < m_scheduleSize; idate++ )
	{	
		mmax = MlEqMaths::Max(mmax, pathArray[0][idate]/refSpot) ;
		mmin = MlEqMaths::Min(mmin, pathArray[0][idate]/refSpot) ;
	}
	
	if(mmin < m_Barrier && mmax > m_Strike)	
		values[0][0] = discountFactors[m_scheduleSize-1]*(mmax - m_Strike) ;
	
	return;
}




void IguanaPricer::setUp(CMatrix& value,MlEqMonteCarlo& mc)
{

	value.resize(1,1);
	m_results.resize(1,1);

	product::setUp(value,mc);
}



