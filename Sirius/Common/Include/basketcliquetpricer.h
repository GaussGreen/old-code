
#pragma once

class MlEqMonteCarlo;
class CMatrix;

class BasketCliquetPricer : public product
{
public:
	CVector	m_localCaps;//[idate]
	CVector	m_localFloors;//[idate]
	CVector	m_strikes;//[idate]
	CVector	m_callput;//[idate]	

	double m_GlobalCap;
	double m_GlobalFloor;
	double m_GlobalRedemption;
	double m_GlobalGearing;
	double m_Notional;

	CVector m_rainbowWeights;		// Rainbow Weights...

	void initialize(	CVector	&localCaps,//[idate]
						CVector	&localFloors,//[idate]
						CVector	&strikes,//[idate]
						CVector	&callput,//[idate]	
						double GlobalCap,
						double GlobalFloor,
						double GlobalRedemption,
						double GlobalGearing,
						double Notional,
						GVector< long >& payoffDates,
						CVector basketWeights
						);


	void payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc);
	void setUp(CMatrix& value,MlEqMonteCarlo& mc);
};



class AsianPricer : public product
{
public:
	void initialize( double Strike, double Cp, const GVector< long >& payoffDates, const CVector& fixings );

	void payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc);
	void setUp(CMatrix& value,MlEqMonteCarlo& mc);

protected:
	int		m_nStart;
	double	m_pastSum;
	int		m_nAvgDates;

	double	m_strike;
	double	m_cp;	

protected:
	double	payout(const CVector& pathArray);
};


