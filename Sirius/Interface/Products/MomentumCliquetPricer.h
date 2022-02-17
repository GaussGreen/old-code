//	Momentum.h : General Cliquet product
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __MOMENTUMCLIQUET_PRICER_H_
#define __MOMENTUMCLIQUET_PRICER_H_

#pragma once

#undef max


class GeneralCliquetPricer : public product
{

protected:
	typedef double (GeneralCliquetPricer::*periodPayoff)(const CVector& perfArray, int nPeriod);
	periodPayoff						m_periodPayoff;
	virtual double								fPeriodPayoff(const CVector& perfArray, int nPeriod);

	double fSum(    const CVector& perfArray, int nPeriod);
	double fProduct(const CVector& perfArray, int nPeriod);
	double fMinimum(const CVector& perfArray, int nPeriod);
	double fMaximum(const CVector& perfArray, int nPeriod);

	virtual void fPointToSpecificPayoff();	// for more exotic payoffs

protected:
	typedef double (GeneralCliquetPricer::*localPayoff)(double, int);
	localPayoff							m_localPayoff;
	double								fLocalReturn(double, int);

	double fForward(double perf, int date){ return perf - m_afStrikes[date]; }
	double fCall(double perf, int date){ return std::max(perf - m_afStrikes[date],0.0); } 
	double fPut(double perf, int date) { return std::max(m_afStrikes[date] - perf,0.0); } 
	double fStraddle(double perf, int date){ return fabs(perf - m_afStrikes[date]); } 

public:

	struct CouponPeriod	// to store values
	{
		long   mPayDate ; 
		double mDiscountFactor;
		double mFloor;
		double mCap;
		double mFixedCoupon;
		double mGearing;
		int mStart;
		int mEnd;
	};

protected:							
	GVector<bool>						m_abStartsPeriod;
	CVector								m_afStrikes;
	CVector								m_afLocalFloors;
	CVector								m_afLocalCaps;	
	CVector								m_afWeights;

	GVector< CouponPeriod >				m_period;
	int									m_nPeriods;	

	double								m_fPastPayoff;
	int									m_nStartPeriod;
	virtual void						init_payout(const CMatrix& pastFixings);

protected:

	virtual void initialize(MlEqAssetHandle				hUnderlying,
							const GVector<long>&		anPayoffDates,		// [date]
							const GVector<bool>&		abStartsPeriod,		// [date]
							const CVector&				afStrikes,			// [date]
							const CVector&				afLocalFloors,		// [date]
							const CVector&				afLocalCaps,		// [date]
							const CVector&				afWeights,			// [date]
							const CVector&				afCoupons,			// [period]
							const CVector&				afGearings,			// [period]
							const CVector&				afPeriodFloors,		// [period]
							const CVector&				afPeriodCaps,		// [period]
							GVector<long>				anPayDatesArr,		// [period]
							MinMaxEnum					minmaxtype,
							PayoffTypeEnum					potype,
							const CMatrix&				pastFixings);		
	
	virtual void payout(CMatrix& values, CMatrix& pathArray, CVector& discountFactors, MlEqMonteCarlo& mc);
	virtual void setUp(CMatrix& value, MlEqMonteCarlo& mc);
};



#endif