
#pragma once

#ifndef	WHALEPRICERH
#define WHALEPRICERH

class MlEqMonteCarlo;
class CMatrix;

class WhalePricer : public product
{
public:
	void initialize(CVector	&weights,
					double strike,
					double cp,
					long nToday,
					const GVector<long>&	payoffDates,
					const CMatrix& fixings	);

	void payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc);
	void setUp(CMatrix& value,MlEqMonteCarlo& mc);

protected:
	CVector	m_weights;
	double m_strike;
	double m_cp;

	double m_runningAverage;
	double m_sumWeights;
	int    m_start;
	int	   m_nDates;
};

#endif