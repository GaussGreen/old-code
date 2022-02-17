
#pragma once

#include "StaticReplication.h"


class CVarianceSwap		:	public CStaticReplication, public RCObject
{
public:
	CVarianceSwap(	MlEqAssetHandle			hUdly,
					DATE					dateStart,
					DATE					dateMaturity,
					MlEqDateScheduleHandle	hFixingSchedule,
					std::string				szCalendar,
					PerformanceTypeEnum		perfType,
					bool					useCurrentSpot,
					StrikesTypeEnum			boundType,
					double					cutUp,
					double					cutDown,
					long					nDaysAYear = 252);


	~CVarianceSwap(){};

	void		fComputePrice(CMatrix& results);

protected:
	long					m_dateStart;		// this is the first observation date
	MlEqDateScheduleHandle	m_hFixingSchedule;
	PerformanceTypeEnum		m_perfType;
	StrikesTypeEnum			m_boundType;
	bool					m_useCurrentSpot;
	std::string				m_szCalendar ;

protected:
	double		fComputePastVariance();
	void		updateMaturity( DATE maturity, double forward );

	double		m_strikeUp;		
	double		m_strikeDown;
	double		m_nDaysStartMat;
	double		m_nDaysTodayMat;
	double		m_nDaysTodayStart;
	double		m_nDaysAYear ;
	double		m_adjStartMat;
	double		m_adjTodayMat;
	double		m_adjTodayStart;
	double		m_forwardStart;

protected:
	double f( double x );
	double df( double x );
	double ddf( double x );	

// the following is useless...
/*
protected:	
	void		shiftSpot( double factor );
	double		fComputeDiscreteVariance(long& nbBusDays, long dayShift);

	void		fComputeDiscretePrice(CMatrix& results, long dayShift);
	void		fComputePriceAndGreeks(CMatrix& result);
	double		fComputePrice(double&, double&);

	void		CalibrateVanilla(const CVector& strikes, const CVector& cp, CVector& weights);
	double		CalculateVanillaPortfolio(const CVector& strikes, const CVector& cp, const CVector& weights);
*/

};


/*
class CVarianceSwapDiscreteSlice		:	public CStaticReplication, public RCObject
{
public:
	CVarianceSwapDiscreteSlice(	MlEqAssetHandle			hUdly,
								DATE					dateMaturity):
	CStaticReplication( hUdly, dateMaturity, 50 ){};

	~CVarianceSwapDiscreteSlice(){};

	double		fComputeCrossTerm( DATE dateStart, DATE dateMaturity, double forwardStart, double forwardMat );
	double		fComputeSquareTerm( DATE dateMaturity, double forwardMat );

protected:		
	bool		m_isCrossTerm;
	bool		m_integrateVolTerm;	

	double f( double x );
	double df( double x );
	double ddf( double x );	
};
*/