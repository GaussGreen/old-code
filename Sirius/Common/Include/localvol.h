
#ifndef __LOCALVOL_H__
#define __LOCALVOL_H__

#include "MlEqInterpolator.h"

class CLocalVolHelper	:	public RCObject // this class will be able to perform interpolation
{
public:
	void initialize( DupireLocalVolHandle lvh );
	double getLocalVol( int slice, double y );
	double getBumpedLocalVol( int slice, double y );

protected:
	long		m_nToday;
	int			m_nSlices;
	CVector		m_forwards;
	double		m_spot;

	RCPtr<MlEqMonotonicSplineInterpolator>	m_interpolator;
	RCPtr<MlEqMonotonicSplineInterpolator>	m_binterpolator;
};



class CBumpVolHelper	:	public RCObject // this class will be able to perform of bump scenarios to guarantee smoothness and no arbitrage...
{
public:
	void initialize( MlEqAssetHandle pAsset, long dateMaturity, const CVector& mktForwards  );
	void	getBumpDerivatives( CVector& bderivs, double reducedStrike, double yearFraction );

protected:
	MlEq2DInterpolatorHandle			m_scenarioBump;
//	CVector		m_refForwards;	// as opposed to market forwards...
};



class FittedLocalVol;
class FittedLocalVol2;

class DupireLocalVol : public RCObject
{

	friend class CLocalVolHelper;	// a helper is a friend...
	friend class FittedLocalVol;
	friend class FittedLocalVol2;

protected:	// local vol calculation members
	int	 m_gridTypeFlag; // 0 is Spot, 1 is Spot / Forward, 2 is log( s/f )
	virtual void computeLocalVol( int t, int k, double rStrike );

protected:	// grid calculation members

	GVector< long >		m_Dates;
	GVector<CVector>	m_localVols;//[idate][ixspace]


	GVector<CVector>	m_bLocalVols;
	CVector				m_times;
	CVector				m_loanspread;
	CVector				m_quantoDrift;
	CVector				m_rate;
	CVector				m_forwards;
	int	m_nSpace;	
	int m_nTime;

	CVector				m_fixedStrikes;
	CVector				m_reducedStrikes;
	GVector<CVector>	m_yStrikes;
	int  m_nInterpolationSpace;

	double	m_volbump ;
	bool	m_vbumped;

protected:
	MlEqAssetHandle		m_pAsset;
	double				m_spot;
	long				m_nToday;

	RCPtr<CLocalVolHelper> m_pLocalVolGrid;	// does nothing if not initialized

	RCPtr<CBumpVolHelper> m_pBumpHelper ;


protected:
	virtual void createTimeGrid( long matDate, int nt );
	virtual void createTimeGrid( const GVector<long>& dates );
	virtual void createTimeGrid( const GVector<long>& dates, int nt );
	void		 createSpaceGrid( double lSpot, double hSpot, int nx );
	void		 createLocalDrift();
	void		 initQuantoDrift();
	virtual void initInterpolationGrid( double stddev = 4., int nxx = 30 );
	void		 createInterpolatedLocalVolGrid();

	double yStrike( int slice, double X );
	double reducedStrike( int slice, double X );

	CVector		m_mktForwards;

public:

	double				getlocalVol(int idate,int ixspace);
	void				setlocalVol(double value,int idate,int ixspace);


	const GVector<CVector>&			getLocalVolGrid()const{return m_localVols;};
	const GVector<CVector>&			getBumpedLocalVolGrid()const{return m_bLocalVols;};
	const CVector&					getRates()const {return m_rate;};
	const CVector&					getLoanSpreads()const {return m_loanspread;};
	const CVector&					getQuantoDrift(){ return m_quantoDrift; };
	const CVector&					getTimes()const {return m_times;};
	const GVector<long>&			getDates()const {return m_Dates;};
	const CVector&					getSpotGrid()const {return m_fixedStrikes;};
	double							getSpot()const {return m_spot;};

	void createLocalVolGrid();
	void bumpLocalDrift(double spotbump);
	void parallelVegaBump();

	void getLocalVolOnGrid( int slice, double low, double up, int nx, CVector& grid, int gridType = 0);

public:	

	DupireLocalVol();
	DupireLocalVol( const DupireLocalVol& rhs );
	~DupireLocalVol(){}

	virtual void initialize(MlEqAssetHandle pAsset,	double lSpot,double hSpot,GVector<long> simDates,int nx, double vb = 0.01, int gridTypeFlag = 2 );// for MC
	virtual void initialize(MlEqAssetHandle pAsset, double lSpot,double hSpot,long mat,int nt,int nx, double vb = 0.01, int gridTypeFlag = 0);// for PDE
	virtual void initialize(MlEqAssetHandle pAsset,	double lSpot,double hSpot,GVector<long> simDates,int nt, int nx, double vb = 0.01, int gridTypeFlag = 0);

	void bumpLocalVol( int t, int k, double bump );
};


inline double DupireLocalVol::yStrike( int slice, double X )
{
	switch( m_gridTypeFlag )
	{
	case 0:	return	log( X / m_forwards[slice] );
	case 1: return	log( X );
	case 2:	return		 X ;
	default:
		throw "wrong grid type defined";
	}
}


inline double DupireLocalVol::reducedStrike( int slice, double X )
{
	switch( m_gridTypeFlag )
	{
	case 0:	return	X / m_forwards[slice] ;
	case 1: return	X ;
	case 2:	return	exp(X) ;
	default:
		throw "wrong grid type defined";
	}
}


/*
class ALLocalVol : public RCObject
{

	GVector<MlEqAnalyticCurveWithTanhWingEdgeHandle> m_localVolSlice;// [icalibDates]
	MlEqDateHandle m_startDate;


public:


	double				getlocalVol(int idate,int ixspace);
	void				setlocalVol(double value,int idate,int ixspace);

	GVector<long> m_calibDates;// [icalibDates]

	void initialize(GVector<MlEqAnalyticCurveWithTanhWingEdgeHandle>& lv);

  	void reinitialize(int it,MlEqAnalyticCurveWithTanhWingEdgeHandle& lv);

	void reinitialize(int it,CVector& lvols);


};
*/



#endif


