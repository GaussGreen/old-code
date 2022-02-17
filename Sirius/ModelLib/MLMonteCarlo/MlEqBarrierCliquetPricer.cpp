//	MlEqBarrierCliquetPricer.cpp : Implementation of MlEqBarrierCliquetPricer
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqObjects.h"
#include "MlEqBarrierCliquetPricer.h"
#include "MonteCarlo.h"
#include "utility.h"
#include "MlEqPde.h"


//#undef max
//#undef min

void MlEqBarrierCliquetPricer::setUp(CMatrix& value, MlEqMonteCarlo& mc)
{
	value.resize(1, 2);
	m_results.resize(1, 2);
	product::setUp(value,mc);

	int nBarrierDates = m_Barrier.rows();

	m_BarrierPrices.resize(m_nNumberOfBarriersPerTimeSlice,nBarrierDates);

	if ( m_modelFlag == 0 ){
		m_stressForwardBarrierVolSlope = m_VolatilityStrikePlace;
	}

	if ( m_modelFlag != 1 ){
		return;
	}

//	m_drift.resize(nBarrierDates);


	const vector<MlEqAssetHandle> & assets = mc.getAssets();

	if ( assets.size() > 1 || assets.size() == 0 ){
		throw("do not know how to deal with this case of zero assets or basket in MlEqBarrierCliquetPricer set up");
	}


	const std::vector<long>& mcDates = mc.GetMCDates();

	m_mapPayoffDateToMcDate.resize(m_payoffDates.getsize());

	int imcdate;
	for ( int idate = 0 ; idate < m_payoffDates.getsize(); idate++ )
	{
		if ( m_payoffDates[idate] < mcDates[0] ){
				continue;
		}

		for ( imcdate = 0 ; imcdate < mcDates.size(); imcdate++ )
		{


			if ( mcDates[imcdate] == m_payoffDates[idate] ){
				m_mapPayoffDateToMcDate[idate] = imcdate;
				break;
			}
		}

		if ( imcdate == mcDates.size() ){
			throw("payoff date not found in barrier cliquet setup");
		}
	}


	int npayoffDates = m_payoffDates.getsize();

	int nToday = mcDates[0];
	bool hitToday = false;
	for ( int idate = 0 ; idate < npayoffDates; idate++ )
	{
		if ( m_payoffDates[idate] == nToday  ){
			hitToday = true;
		}
	}
	int offset;
	if ( hitToday ){
		offset = 0;
	}else{
		offset = 1;
	}

	int mcdate = 0;
	if ( !hitToday )
	{
		int idateBarrier = 0;
		for ( int idate = 0 ; idate < npayoffDates-1; idate++ )
		{
			if ( !m_hasbarrier[idate] ){
				continue;
			}

			if ( m_BarrierStartDates[idateBarrier] < nToday && m_BarrierEndDates[idateBarrier] > nToday ){
				m_mapPayoffDateToMcDate[idate] = 0;
			}

			idateBarrier++;
		}
	}



}





void MlEqBarrierCliquetPricer::payout(CMatrix& values, CMatrix& pathArray, CVector& discountFactors, MlEqMonteCarlo& mc)
{
	CMatrix HitProb;
	double price;

	int ndates = pathArray.cols();

	double payoff=0;	
	int idateBarrier = 0;


	for(int idate=0;idate < ndates-1; idate++) 
	{

		if(!m_hasbarrier[idate] || m_BarrierPayDates[idateBarrier] < m_nToday){
			continue;
		}

		for(int ibarrier=0;ibarrier<m_nNumberOfBarriersPerTimeSlice;ibarrier++) 
		{

			if ( m_BarrierEndDates[idate] <= m_nToday )
			{
				if ( m_BarrierHasHit[idateBarrier][ibarrier] )
				{
					payoff += m_Rebate[idateBarrier][ibarrier] * m_discount[idateBarrier] ;
				}
			}
			else 
			{
				int mcidate = m_mapPayoffDateToMcDate[idate];

				price = mc.AmericanDigital(m_BarrierPrices[ibarrier][idateBarrier],mc.GetCurrentPath(),mcidate,
										   m_Barrier[idateBarrier][ibarrier],mc.GetPathValue(mc.GetCurrentPath(),mcidate+1)/mc.GetPathValue(mc.GetCurrentPath(),mcidate),
										   0,m_modelFlag,m_stressForwardBarrierVol,
								           m_stressForwardBarrierVolSlope,m_npoints,m_blend);

				if(!m_bKnockIn){
					price = 1.- price;
				}

				payoff += price * /*m_Rebate[idateBarrier][ibarrier] */ m_discount[idateBarrier] ; 

			}

		}
					
		idateBarrier++;
	}

	values[0][0] = payoff;

 	fillValues(values,pathArray,mc);

}



void MlEqBarrierCliquetPricer::initialize( MlEqZeroCurveHandle hZCurve,

											CMatrix&					Barriers,
											CMatrix&					Rebate,
											GMatrix< bool > &			BarrierHasHit,
											GVector< long>				BarrierStartDates,
											GVector< long>				BarrierEndDates,
											GVector< long>				BarrierPayDates,

											int							pricingMethod,
											double						stressForwardBarrierVol,
											double						stressForwardBarrierVolSlope,
											double						VolatilityStrikePlace,
											double						blend,
											int							npoints,
											MlEqConstDateHandle			hDate,
											CVector&					bfixings,
											double						spot,
											bool						bKnockIn

											) 
{

	
	m_nToday						= hDate->GetDate();

	m_Barrier						= Barriers;
	m_Rebate						= Rebate;
	m_VolatilityStrikePlace			= VolatilityStrikePlace;
	m_blend							= blend;
	m_npoints						= npoints;
	m_nNumberOfBarriersPerTimeSlice	= Barriers.cols();
	m_modelFlag						= pricingMethod;
	m_stressForwardBarrierVol		= stressForwardBarrierVol;
	m_stressForwardBarrierVolSlope	= stressForwardBarrierVolSlope;
	m_BarrierHasHit					= BarrierHasHit;
	m_BarrierPayDates				= BarrierPayDates;
	m_BarrierStartDates				= BarrierStartDates;
	m_BarrierEndDates				= BarrierEndDates;

	m_bKnockIn						= bKnockIn;
	

	for ( int i = 0 ; i < m_BarrierHasHit.rows(); i++ )
	{
		for ( int ibarrier = 0 ; ibarrier < m_BarrierHasHit.cols(); ibarrier++ )
		{
			if ( m_BarrierHasHit[i][ibarrier] && BarrierStartDates[ibarrier] > m_nToday ){
				throw("barrier has hit cannot be true for a barrier that has not even started");
			}
		}
	}

	for ( int idate = 0 ; idate < BarrierStartDates.getsize(); idate++ )
	{
		if ( BarrierEndDates[idate] <= BarrierStartDates[idate] ){
			throw("barrier end date must be greater than barrier start date");
		}

		if ( idate )
		{
			if ( BarrierStartDates[idate] <=  BarrierStartDates[idate-1] ){
				throw("barrier startdates must be in ascending order");
			}

			if ( BarrierEndDates[idate] <=  BarrierEndDates[idate-1] ){
				throw("barrier enddates must be in ascending order");
			}
		}
	}


	int npayoffDates = 1;
	for (int idate = 0; idate < BarrierStartDates.getsize(); idate++ )
	{

		npayoffDates++;

		if ( idate ){
			if ( BarrierStartDates[idate] > BarrierEndDates[idate-1] ){
				npayoffDates++;
			}
		}
	}


	m_payoffDates.resize(npayoffDates);

	m_hasbarrier.resize(npayoffDates-1);
	
	int n = 0;
	int k = 0;
	m_payoffDates[k]	= BarrierStartDates[0];

	k++;

	for (int idate = 0; idate < BarrierStartDates.getsize(); idate++ )
	{

		n++;

		m_payoffDates[k]	= BarrierEndDates[idate];
		m_hasbarrier[k-1]		= true;

		k++;


		if ( idate < BarrierStartDates.getsize()-1)
		{
			if ( BarrierStartDates[idate+1] > BarrierEndDates[idate] )
			{
				m_hasbarrier[k-1]			= false;
				m_payoffDates[k]		= BarrierStartDates[idate+1];
				k++;


			}
		}


	}

	bool hitToday = false;
	for ( int idate = 0 ; idate < npayoffDates; idate++ )
	{
		if ( m_payoffDates[idate] == m_nToday  ){
			hitToday = true;
		}
	}
	int offset;
	if ( hitToday ){
		offset = 0;
	}else{
		offset = 1;
	}

	int mcdate = 0;
/*	m_mapPayoffDateToMcDate.resize(npayoffDates);
	for ( int idate = 0 ; idate < npayoffDates; idate++ )
	{
		if ( m_payoffDates[idate] < nToday ){
			m_mapPayoffDateToMcDate[idate] = -1000;
			continue;
		}

		if ( m_payoffDates[idate] == nToday  ){
			 m_mapPayoffDateToMcDate[idate] = mcdate;

			 mcdate++;
		}
		
		if ( m_payoffDates[idate] > nToday  ){
			 m_mapPayoffDateToMcDate[idate] = mcdate+offset;
			 mcdate++;
		}
	}



	if ( !hitToday )
	{
		int idateBarrier = 0;
		for ( int idate = 0 ; idate < npayoffDates-1; idate++ )
		{
			if ( !m_hasbarrier[idate] ){
				continue;
			}

			if ( BarrierStartDates[idateBarrier] < nToday && BarrierEndDates[idateBarrier] > nToday ){
				m_mapPayoffDateToMcDate[idate] = 0;
			}

			idateBarrier++;
		}
	}

*/

//	barriers can be entered either as dollar values or in percent;
//  future starting barriers however must be entered as percentage value
//  not that as far as past barriers are concerned we don't care because valuation depends only on 
//  m_BarrierHasHit flag
//  hence only case that matters is of type is where barrier started in the past and extends into the future

	int idateBarrier = 0 ;
	for (int idate = 0; idate < npayoffDates-1; idate++ )
	{
		if ( m_hasbarrier[idate] )
		{
			if ( BarrierStartDates[idateBarrier] <= m_nToday && BarrierEndDates[idateBarrier]> m_nToday )
			{
				for ( int ibarrier = 0 ; ibarrier < m_nNumberOfBarriersPerTimeSlice; ibarrier++ ){
						m_Barrier[idate][ibarrier] *= bfixings[idateBarrier]/spot;
				}
			}

			idateBarrier++;
		}
	}


	m_discount.resize(idateBarrier);
	for (int idate = 0; idate < idateBarrier; idate++ )	{	
		m_discount[idate] = hZCurve->GetDiscountFactor(m_nToday, m_BarrierPayDates[idate]);
		
	}



}