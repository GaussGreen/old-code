//	NapoleonPricer.h : Napoleon product pricer
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __NAPOLEONPRICER_H_
#define __NAPOLEONPRICER_H_

#include "MlEqInterpolator.h"

class MlEqNapoleonPricer : public product
{
protected:							
	GVector<bool>						m_abStartsPeriod;
	GVector<long>						m_anCallPuts;
	CVector								m_afStrikes;
	CVector								m_afLocalFloors;
	CVector								m_afLocalCaps;	
	CVector								m_afWeights;		
	CVector								m_afCoupons;
	CVector								m_afGearings;	
	CVector								m_afPeriodFloors;
	CVector								m_afPeriodCaps;
	int									m_nPeriods;						// number of 'cliquette' steps in the Napoleon
	CVector								m_afPayDiscountFactors;			// discount factors corresponding to each payment event	

	void initialize(MlEqAssetHandle				hUnderlying,
					const GVector<long>&		anPayoffDates,
					const GVector<bool>&		abStartsPeriod,		// [date]
					const GVector<long>&		anCallPuts,			// [date]
					const CVector&				afStrikes,			// [date]
					const CVector&				afLocalFloors,		// [date]
					const CVector&				afLocalCaps,		// [date]
					const CVector&				afWeights,			// [date]
					const CVector&				afCoupons,			// [period]
					const CVector&				afGearings,			// [period]
					const CVector&				afPeriodFloors,		// [period]
					const CVector&				afPeriodCaps,		// [period]
					GVector<long>				anPayDatesArr,		// [period, optional]
					bool						bDelayPaymentsToEnd);
	
	void payout(CMatrix& values, CMatrix& pathArray, CVector& discountFactors, MlEqMonteCarlo& mc);
	void setUp(CMatrix& value, MlEqMonteCarlo& mc);
};



class DeltaHedgedOptionPortfolio
{

	friend void integratePortfolio(double, double[],void *vp);
	friend void integrateExpectedPnlRungeKutta(double x, double z[]  ,void *vp);

protected:

	bool						m_gaussIntegrator;
	double						m_nstdev;
	double						m_eps;
	long						m_approximationDate;
	int							m_ntimeGauss;

	double						m_sq;
	bool						m_calcAddOn;

protected:

	MlEqAssetHandle				m_hAsset;

	GVector< MlEqDateHandle >	m_recordingDates;//[idateRecord]
	CVector						m_divsToMaturity;//[idateRecordRemaining]
	CVector						m_ratesToMaturity;//[idateRecordRemaining]
	int							m_previousSettlementIndex;

	MlEqInterpolatorHandle		m_bookingRates;
	MlEqInterpolatorHandle		m_bookingDivs;

	double						m_bookingRate;
	double						m_bookingDiv;

	GVector<int>				m_mapToSettlement;//[ipnldate]
	GVector<int>				m_mapToRecordingDate;//[isettlement]

	long						m_nMaturity;


	int							m_ngauss;
	CVector						m_gaussPoint;
	CVector						m_gaussWeight;
	int							m_nDates;

	double		calcPortfolio(int idate,MlEqAssetHandle hAsset,double spot);
	double		calcPnl(int idate,MlEqAssetHandle hAsset,double spot,double refSpot);
	double		calcExpectedPnlGauss(bool calcAddOn,int idate,MlEqAssetHandle hAsset,double refSpot);
	double		calcExpectedPnlGauss(CMatrix& values,bool calcAddOn,int idateStart,int idateEnd,MlEqAssetHandle hAsset,double nstdev);
	double		calcExpectedPnlRungeKutta(bool calcAddOn,int idateStart,int idateEnd,MlEqAssetHandle hAsset,double nstdev,double eps);

	double		calcExpectedPnl(CMatrix& values,bool calcAddOn,bool gaussIntegrator,int idateStart,int idateEnd,MlEqAssetHandle hAsset,double nstdev,double eps);
	double		calcPnlB(bool calcAddOn,int idate,MlEqAssetHandle hAsset,double spot,double refSpot);
	double		calcPnlC(bool calcAddOn,int idate,MlEqAssetHandle hAsset,double refSpot);

	double calcPortfolioB(double& deltaTerm,bool calcPortfolioB,int idate,MlEqAssetHandle hAsset,double spot,double refSpot);


public:


	long						m_startDate;
	GVector<long>				m_settlementDates;//[idate]  last one being the maturity

//  the following is needed for the P/L recording freuqency

	std::string					m_szTenor;
	BusinessDayConventionEnum	m_bdc;
	std::string					m_szCalendar;
	DayCountConventionEnum		m_dcc;

	CMatrix						m_portfolio;//[iweight][icall/iput]
	CVector						m_bookingVols;//[iweight];
	long						m_useControlVariate;
	double						m_includeDeltaTerm;




	double price(CMatrix& values,MlEqAssetHandle hAsset);

	



	void initialise(	
						MlEqAssetHandle hAsset,
						long startDate,
						GVector<long>& settlementDates,
						std::string	&szTenor,
						BusinessDayConventionEnum	bdc,
						std::string					szCalendar,
						CMatrix	&					portfolio,
						CVector	&					bookingVols,
						double						bookingRate,
						double						bookingDiv,
						DayCountConventionEnum		dcc,
						long						includeDeltaTerm,
						bool						gaussIntegrator,
						int							ngauss,
						double						nstdev,
						double						eps,
						long						approximationDate,
						int							ntimeGauss
						);


};


#endif