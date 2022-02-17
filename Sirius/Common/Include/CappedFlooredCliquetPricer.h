
#ifndef	CAPPEDFLOOREDCLIQUETPRICERH
#define CAPPEDFLOOREDCLIQUETPRICERH

class MlEqMonteCarlo;
class CMatrix;

class CappedFlooredCliquetPricer : public product
{
protected:

	double calculateReturn(CForwardSkewMC& mc,int idate,int ipath);
	double calculateReturn(CMatrix& pathArray,int idate);

public:
	CVector	m_localCaps;//[idate]
	CVector	m_localFloors;//[idate]
	CVector	m_strikes;//[idate]
	CVector	m_callput;//[idate]	
	CVector	m_weights;//[idate]

	double m_GlobalCap;
	double m_GlobalFloor;
	double m_GlobalRedemption;
	double m_GlobalGearing;
	double m_Notional;

	CVector m_basketWeights;

	double m_sumOfBasketWeights;

	void initialize(	CVector	&localCaps,//[idate]
						CVector	&localFloors,//[idate]
						CVector	&strikes,//[idate]
						CVector	&callput,//[idate]	
						CVector	&weights,//[idate]
						double GlobalCap,
						double GlobalFloor,
						double GlobalRedemption,
						double GlobalGearing,
						double Notional,
						GVector< long >& payoffDates,
						CVector basketWeights,
						bool rebalancingBasket //= true
						);

	void initialize(	CVector	&localCaps,//[idate]
						CVector	&localFloors,//[idate]
						CVector	&strikes,//[idate]
						CVector	&callput,//[idate]	
						CVector	&weights,//[idate]
						double GlobalCap,
						double GlobalFloor,
						double GlobalRedemption,
						double GlobalGearing,
						double Notional,
						GVector< long >& payoffDates,
						CVector basketWeights,
						const CMatrix& LockinSchedule,
						bool rebalancingBasket //= true
						);

	void payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc);
	void evaluate(CMatrix& results,CForwardSkewMC& mc);
	void setUp(CMatrix& value,MlEqMonteCarlo& mc);

protected:
	void payoutClassic(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc);	// old feature..
	void payoutLockin(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc);	// new feature..

protected:		// new payoff specific
	CMatrix		m_LockinSchedule ;
	bool		m_LockinCapped ;
	bool		m_LockinPayoff ;
	bool		m_rebalancingBasket ;
};


class Himalaya : public product
{

	GVector<bool>	m_exclude;
	double			selectHimalayaWinner(int idate,CVector& values,int& currentWinner);
	bool			callHimalaya(double& payoff,int idate,CVector& values,CVector& discountFactors,double divisor);
	void			callRainbow(double & res,CVector& returns);

	bool			m_asianFromStart;
	bool			m_isHimalaya;

	CVector			m_spots;

public:

	CVector			m_rainbowWeights;
	GVector<long>	m_isHimalayaDate;// [m_nDates]
	double			m_strike;
	int				m_nAssets;
	double			m_Notional;
	double			m_callPut;
	CVector			m_rebate;
	bool			m_delayRebate;
	CVector			m_callLevels;//[m_nDates]
	GVector<bool>	m_isCallable;
		
	void payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc);
	void setUp(CMatrix& value,MlEqMonteCarlo& mc);

	void initialize(int nAssets, CVector& spotLevels,GVector<long> payoffDates, GVector<long> isHimalayaDate, CVector rainbowWeights, double notional, double strike, double callPut, CVector m_callLevels, CVector rebate, bool delayRebate, bool isAsianFromStart);
};



class CAsianMax : public product
{

	public:

	CVector		m_spotFixings;
	CVector m_strike;//[idate]
	double m_callPut;

	void initialize(double callPut,CVector& strikes, GVector<long> Dates,CVector& spotFixings);

	void payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc);
//	void evaluate(CMatrix& results,CForwardSkewMlEqMonteCarlo& mc);
	void setUp(CMatrix& value,MlEqMonteCarlo& mc);
};


#endif
