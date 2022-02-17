//	MlEqBarrierCliquetPricer.h : Barrier Cliquet product pricer
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __MLEQBARRIERCLIQUETPRICER_H_
#define __MLEQBARRIERCLIQUETPRICER_H_

class MlEqBarrierCliquetPricer : public product
{
protected:							
	

	CMatrix								m_Rebate;//[idateBarrier][ibarrier]
	CMatrix								m_Barrier;//[idateBarrier][ibarrier]
	double								m_VolatilityStrikePlace;
	int									m_MCnDates;
	double								m_blend;
	double								m_pastcontributions;
	GVector<int>						m_hasbarrier;
	int									m_nDates;
	int									m_nNumberOfBarriersPerTimeSlice;
	vector <double>						m_DiscountFactors;

	int									m_npoints;
//	CVector								m_drift;
	int									m_modelFlag;
	double								m_stressForwardBarrierVol;
	double								m_stressForwardBarrierVolSlope;
	GMatrix < CVector > 				m_BarrierPrices;//[ibarrier][idateBarrier][ivolstate]
	GMatrix	< bool >					m_BarrierHasHit;//[ndate][nbarrier]
	GVector	< long >					m_BarrierPayDates;//[ibarrierDate]
	GVector	< long >					m_BarrierStartDates;//[ibarrierDate]
	GVector	< long >					m_BarrierEndDates;//[ibarrierDate]

	GVector	< long >					m_mapPayoffDateToMcDate;

	bool								m_bKnockIn;
//	bool								m_bPayAtEnd;
	CVector								m_discount;
	long								m_nToday;


	void initialize(
					CMatrix&					init_Barrier,
					CMatrix&					init_Rebate,
					vector<bool>&				init_hasbarrier,
					double						init_VolatilityStrikePlace,
					int							init_nDates,
					int							init_MCnDates,
					int							init_nRuns,
					vector<double>&				init_discountfactors,
					double						init_blend,
					double						init_pastcontributions
					);

	void initialize(	MlEqZeroCurveHandle			hZCurve,

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
						CVector		&				bfixings,
						double						spot,
						bool						bKnockIn

						) ;


	void payout(CMatrix& values, CMatrix& pathArray, CVector& discountFactors, MlEqMonteCarlo& mc);
	void setUp(CMatrix& value, MlEqMonteCarlo& mc);
};

#endif


