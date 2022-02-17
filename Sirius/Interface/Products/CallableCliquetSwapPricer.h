
#pragma once

class MlEqMonteCarlo;
class CMatrix;

class CallableCliquetSwapPricer : public product
{
public:

	struct CouponPeriod	// only to store values
	{
		long   mPayDate ; 
		double mDiscountFactor;
		double mFloor;
		double mCap;
		double mFixedCoupon;
		double mFunding;
		int mStart;
		int mEnd;
	};

protected:

	CVector	m_localCaps;
	CVector	m_localFloors;
	CVector	m_strikes;
	CVector	m_callput;
	CVector	m_weights;

	double  mCallableLevel;

	int mNumberOfPeriods;
	std::vector< CouponPeriod > m_period;

	void init_payout(CMatrix& pathArray, long today);
	double init_couponSum ;
	double init_paidCoupons ;
	double mCurrentCliquet;
	int mCurrentPeriod ;
	bool isExpired ;

public:
	CallableCliquetSwapPricer(){};
	~CallableCliquetSwapPricer(){};

	void initialize(	CVector	&localCaps,
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
			/*			MlEqZeroCurveHandle hZeroCurve,
						GVector< long >& fundingPayDates,
						double FundingVAlue,
			*///			MlEqSwapHandle hswap,
						CMatrix& pastFixings
						);

	void payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors, MlEqMonteCarlo& mc);
	void setUp(CMatrix& value, MlEqMonteCarlo& mc);
};

class AutoCallableSwapPricer : public product
{
protected:

	CVector	m_thresholds;
	CVector	m_coupon_low;
	CVector	m_coupon_up;

	CVector	m_rainbow_weights ;

protected:

	CVector	m_barriers;
	CVector	m_rebates;
	CVector	m_libor_over_period;

	double	m_strike;
	double	m_callput;
	double	m_gearing;
	double	m_fixed_coupon;

	bool m_isExpired ;
	int m_nStart;
	int m_nBarrierDates;

public:
	AutoCallableSwapPricer(){};
	~AutoCallableSwapPricer(){};

	void initialize(	const 	CVector	&rainbow_weights,	
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
						CMatrix& pastFixings
						);

	void payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc);
	void setUp(CMatrix& value,MlEqMonteCarlo& mc);
};


