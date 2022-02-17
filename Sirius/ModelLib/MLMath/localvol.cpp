
#include "stdafx.h"
#include "mleqobjects.h"
#include "localvol.h"



/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: DupireLocalVol
**	Returns: 
**	Action : 
**           
****************************************************************/


DupireLocalVol::DupireLocalVol()
{
}



/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: DupireLocalVol
**	Returns: 
**	Action : 
**           
****************************************************************/

DupireLocalVol::DupireLocalVol( const DupireLocalVol& rhs )
{
	m_gridTypeFlag	=	rhs.m_gridTypeFlag; 
	m_Dates			=	rhs.m_Dates;;
	m_localVols		=	rhs.m_localVols;
	m_bLocalVols	=	rhs.m_bLocalVols;
	m_times			=	rhs.m_times;
	m_loanspread	=	rhs.m_loanspread;
	m_rate			=	rhs.m_rate;
	m_forwards		=	rhs.m_forwards;
	m_nSpace		=	rhs.m_nSpace;	
	m_nTime;

	m_fixedStrikes	=	rhs.m_fixedStrikes;
	m_reducedStrikes
					=	rhs.m_reducedStrikes;
	m_yStrikes		=	rhs.m_yStrikes;
	m_nInterpolationSpace
					=	rhs.m_nInterpolationSpace;

	m_volbump		=	rhs.m_volbump ;
	m_vbumped		=	rhs.m_vbumped;

	m_pAsset		=	rhs.m_pAsset;
	m_spot			=	rhs.m_spot;
	m_nToday		=	rhs.m_nToday;

	m_pLocalVolGrid	=	rhs.m_pLocalVolGrid;
	
	m_quantoDrift	=	rhs.m_quantoDrift;
}


/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: initialize
**	Returns: 
**	Action : 
**           
****************************************************************/

void DupireLocalVol::initialize(MlEqAssetHandle pAsset,
								double lowSpot,
								double highSpot,
								GVector<long> simDates,
								int nt,
								int nx, 
								double vb,
								int gridTypeFlag )
{
	m_pAsset	= pAsset;
	m_nSpace	= nx;
	m_volbump	= vb;
	m_nToday	= m_pAsset->GetDateHandle()->GetDate();
	m_spot		= pAsset->GetSpot(m_nToday) ;

	m_gridTypeFlag = gridTypeFlag ;

	createTimeGrid( simDates, nt );
	m_nTime		= m_times.getsize();

	createSpaceGrid( lowSpot, highSpot, nx );
	createLocalDrift();	

	initQuantoDrift();

	m_localVols.resize(m_nTime);
	m_bLocalVols.resize(m_nTime);

	if( !m_pLocalVolGrid ){
		initInterpolationGrid();
	}

	createInterpolatedLocalVolGrid();

	m_vbumped = false;
}


/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: initialize
**	Returns: 
**	Action : 
**           
****************************************************************/

void DupireLocalVol::initialize(MlEqAssetHandle pAsset,
								double lowSpot,
								double highSpot,
								GVector<long> simDates,
								int nx, 
								double vb,
								int gridTypeFlag )
{
	m_pAsset	= pAsset;
	m_nSpace	= nx;
	m_volbump	= vb;
	m_nToday	= m_pAsset->GetDateHandle()->GetDate();
	m_spot		= pAsset->GetSpot(m_nToday) ;

	m_gridTypeFlag = gridTypeFlag ;

	m_nTime		= simDates.getsize();
	createTimeGrid( simDates );

	createSpaceGrid( lowSpot, highSpot, nx );
	createLocalDrift();	
	initQuantoDrift();

	m_localVols.resize(m_nTime);
	m_bLocalVols.resize(m_nTime);

	if( !m_pLocalVolGrid ){
		initInterpolationGrid();
	}

	createInterpolatedLocalVolGrid();

	m_vbumped = false;
}


/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: initialize
**	Returns: 
**	Action : 
**           
****************************************************************/

void DupireLocalVol::initialize(MlEqAssetHandle pAsset,
								double lowSpot, 
								double highSpot,
								long maturityDate,
								int nt,
								int nx, 
								double vb,
								int gridTypeFlag )
{
	m_pAsset	= pAsset;
	m_nSpace	= nx;
	m_volbump	= vb;
	m_nToday	= m_pAsset->GetDateHandle()->GetDate();
	m_spot		= pAsset->GetSpot(m_nToday) ;

	m_gridTypeFlag = gridTypeFlag ;

	m_nTime		= nt;
	createTimeGrid( maturityDate, m_nTime );

	createSpaceGrid( lowSpot, highSpot, nx );
	createLocalDrift();

	initQuantoDrift();

	m_localVols.resize(m_nTime);
	m_bLocalVols.resize(m_nTime);

	if( !m_pLocalVolGrid ){
		initInterpolationGrid();
	}

	createInterpolatedLocalVolGrid();
		
	m_vbumped = false;
}


/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: createTimeGrid
**	Returns: 
**	Action : 
**           
****************************************************************/

void DupireLocalVol::createTimeGrid( const GVector<long>& dates, int nt )
{
	MlEqConstDateHandle	dateh = m_pAsset->GetDateHandle();

	int addnt = dates.getsize();
	long maturityDate = dates[addnt-1];

	std::vector<long> sdates;
	sdates.push_back(m_nToday);

	double T = double( maturityDate - m_nToday ) ;
	
	long date = m_nToday;
	long nextDate = dates[1], prevDate = dates[0];

	for ( int i = 1 ; i < addnt; i++ )
	{
		nextDate = dates[i];

		if( nextDate > prevDate && nextDate > m_nToday)
		{
			double Tlap = double(nextDate - prevDate);
			int nlapt = int(nt *  Tlap / T);
			double dt = Tlap / double(nlapt);
			date = long( date + dt );

			while( date > prevDate && date < nextDate )
			{			
				sdates.push_back( date );
				date = long( date + dt );
			}
			sdates.push_back( nextDate );
		}
		date = prevDate = nextDate;
	}

	
	int totalnt = sdates.size();
	m_times.resize(totalnt);
	m_Dates.resize(totalnt);

	for(int i=0; i<totalnt; ++i)
	{
		date = sdates[i];
		m_Dates[i] = date;
		m_times[i] = dateh->GetYearFraction( date );
	}
}


/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: createTimeGrid
**	Returns: 
**	Action : 
**           
****************************************************************/


void DupireLocalVol::createTimeGrid( long maturityDate, int nt )	// we should have a time grid class...
{
	MlEqConstDateHandle	dateh = m_pAsset->GetDateHandle();

	m_times.resize(nt);
	m_Dates.resize(nt);

	if ( nt >= maturityDate-m_nToday )
		throw("case not implemented yet");

	m_Dates[0] = m_nToday;
	m_times[0] = 0.;
	double dt = double( maturityDate - m_nToday ) / (nt-1.) ;
	for ( int i = 1 ; i < nt; i++ )
	{
		m_Dates[i] = long( m_nToday + i*dt );
		m_times[i] = dateh->GetYearFraction(m_Dates[i]);
	}
	m_Dates[nt-1] = maturityDate;
	m_times[nt-1] = dateh->GetYearFraction( maturityDate );
}


/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: createTimeGrid
**	Returns: 
**	Action : 
**           
****************************************************************/


void DupireLocalVol::createTimeGrid( const GVector<long>& dates )
{
	m_nTime = dates.getsize();
	m_Dates = dates;
	m_times.resize(m_nTime);

	MlEqConstDateHandle	dateh = m_pAsset->GetDateHandle();

	for (int i=0 ; i<m_nTime; i++){
		m_times[i] = dateh->GetYearFraction(m_Dates[i]);
	}
}

/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: createSpaceGrid
**	Returns: 
**	Action : 
**           
****************************************************************/


void DupireLocalVol::createSpaceGrid( double lowSpot, double highSpot, int nx )
{
	m_nSpace = nx;
	m_fixedStrikes.resize(m_nSpace);

	if( m_gridTypeFlag == 0 )
	{
		double ux = log(highSpot/m_spot);
		double lx = log(lowSpot/m_spot);		
		double dlogspot = (ux - lx) / (m_nSpace-1.);

		for ( int i = 0 ; i < m_nSpace; i++ ){	
			m_fixedStrikes[i] = m_spot * exp( lx + i*dlogspot );
		}
	}
	else if( m_gridTypeFlag == 1 )
	{
		double ux = log(highSpot);
		double lx = log(lowSpot);			
		double dlogspot = (ux - lx) / (m_nSpace-1.);

		for ( int i = 0 ; i < m_nSpace; i++ ){	
			m_fixedStrikes[i] = exp( lx + i*dlogspot );
		}
	}
	else if( m_gridTypeFlag == 2 )
	{	
		double ux = highSpot;
		double lx = lowSpot;			
		double dlogspot = (ux - lx) / (m_nSpace-1.);

		for ( int i = 0 ; i < m_nSpace; i++ ){	
			m_fixedStrikes[i] = lx + i*dlogspot ;
		}
	}
}


/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: createLocalDrift
**	Returns: 
**	Action : 
**           
****************************************************************/

void DupireLocalVol::createLocalDrift()
{
	const MlEqZeroCurveHandle	curve = m_pAsset->GetPayZeroCurve(true);

	m_loanspread.resize(m_nTime);
	m_rate.resize(m_nTime);
	m_forwards.resize(m_nTime);
		 
	m_forwards[0] = m_spot;

	double dt, mu, r, disc, discprev = 1., fwd, fwdprev = m_spot, t, tprev = 0.;

	for ( int idate = 0 ; idate < m_nTime-1; idate++ )
	{
		fwd		   = m_pAsset->GetNaturalForward(m_nToday,m_Dates[idate+1], false);
		disc	   = curve->GetDiscountFactor(m_Dates[idate+1]);
		t		   = m_times[idate+1];

		dt  = t - tprev;	
		mu	= fwd / fwdprev;
		r	= disc / discprev;	 	

		m_loanspread[idate] = log(mu)/dt;
		m_rate[idate]		= -log(r)/dt;

//		if ( m_rate[idate] < 0.0 ){
//			throw("negative discount rate encountered");	// this may actually happen in Japan !!!
//		}

		m_forwards[idate+1] = fwd ;
	
		fwdprev = fwd;	
		discprev = disc;	
		tprev = t;
	}

	m_loanspread[m_nTime-1]	= m_loanspread[m_nTime-2];
	m_rate[m_nTime-1]		= m_rate[m_nTime-2];

	initQuantoDrift();
}


/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: initInterpolationGrid
**	Returns: 
**	Action : 
**           
****************************************************************/

void DupireLocalVol::initInterpolationGrid( double stdDev, int nx )
{
	m_nInterpolationSpace	= nx;

// build space grid
	long maturityDate = m_Dates[m_nTime-1];

	double stdy = stdDev * m_pAsset->GetVolatility( MlEqStrike(m_spot), m_nToday, maturityDate);
	CVector nstrike(m_nInterpolationSpace);

	double dstd = 2.5/ double(m_nInterpolationSpace-1);

	for (int j=0 ; j<m_nInterpolationSpace; j++){
		nstrike[j] = stdy * ( dstd*j - 1.5);
	}

	m_yStrikes.resize(m_nTime);
	for ( int i = 0 ; i < m_nTime; i++ )
	{
		m_yStrikes[i].resize(m_nInterpolationSpace);
		double sqT = sqrt( std::max(m_times[i], 0.03) );

		for ( int j = 0 ; j < m_nInterpolationSpace; j++ )
		{
			m_yStrikes[i][j] = nstrike[j] * sqT ;
		}
	}


	std::vector<double> volDates = m_pAsset->GetVolatilityStructure()->getDates();
	int nVolDates = volDates.size();
	m_mktForwards.resize(nVolDates);

	for(int i=0; i<nVolDates; ++i)	{
		m_mktForwards[i] = m_pAsset->GetNaturalForward(m_nToday, volDates[i], false);
	}

	// init bump interpolator,,,

//	m_pBumpHelper = new CBumpVolHelper();
//	m_pBumpHelper->initialize( m_pAsset, maturityDate, m_mktForwards );
	

// build lv grid
	m_localVols.resize(m_nTime);
	m_bLocalVols.resize(m_nTime);

	for(int t=0; t<m_nTime; t++)
	{
		m_localVols[t].resize(m_nInterpolationSpace);
		m_bLocalVols[t].resize(m_nInterpolationSpace);

		for(int k=0; k<m_nInterpolationSpace; k++)
		{ 
			double rstrike = exp( m_yStrikes[t][k] );
			computeLocalVol(t,k, rstrike);
		}
	}

	m_pLocalVolGrid = new CLocalVolHelper;
	m_pLocalVolGrid->initialize( this );
}


/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: createInterpolatedLocalVolGrid
**	Returns: 
**	Action : 
**           
****************************************************************/


void DupireLocalVol::createInterpolatedLocalVolGrid()
{
	for(int t=0; t<m_nTime; t++)
	{
		m_localVols[t].resize(m_nSpace);
		m_bLocalVols[t].resize(m_nSpace);

		for(int k=0; k<m_nSpace; k++)
		{
			double y = yStrike( t, m_fixedStrikes[k] );

			m_localVols[t][k]  =	m_pLocalVolGrid->getLocalVol( t, y );	
			m_bLocalVols[t][k] =	m_pLocalVolGrid->getBumpedLocalVol( t, y );
		}
	}
}


/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: computeLocalVol
**	Returns: 
**	Action : 
**           
****************************************************************/



void DupireLocalVol::computeLocalVol( int t, int k, double reducedStrike )
{
	if( t<0 || t>m_nTime-1)	throw "Wrong time index";

	MlEqVolatilityStructureHandle hVolStruct = m_pAsset->GetVolatilityStructure();
	
	CVector res;
	double T = m_times[t];
	long date = m_Dates[t];
	double fwd = m_forwards[t];

	if( T < 1e-4 )
	{
		double atmvol = hVolStruct->getVol( MlEqStrike(fwd), DATE(date) ) ;

		m_localVols[t][k]  = atmvol;
		m_bLocalVols[t][k] = atmvol + m_volbump;
		return;
	}	
	
	
	double dsdk, d2sd2k, dsdt, vol, sqT  = sqrt(T);

	double rk =	reducedStrike; 


	MlEqForwardBasedStrike fstk(rk, fwd, date);
	hVolStruct->getVolDeriv( res, fstk, date, m_mktForwards) ;
	dsdk	= res[0];
	d2sd2k	= res[1];
	vol		= res[2];
	dsdt	= res[3];
/*
	CVector bres(4);
	m_pBumpHelper->getBumpDerivatives( bres, rk, T );
	dsdk	+= bres[0];
	d2sd2k	+= bres[1];
	vol		+= bres[2];
	dsdt	+= bres[3];
*/

	double SsqT = vol*sqT ;
	double d	= -log(rk)/SsqT + 0.5*SsqT ;

	double num = vol * ( vol + 2.* T * dsdt ) ;

	double denom = 1.+ rk * d * dsdk * sqT ;
	denom *= denom ;

	double temp = rk*rk *vol*T * ( d2sd2k - d*sqT * dsdk*dsdk );

	denom += temp;

	double local_vol =  num/denom ;
	local_vol = sqrt( std::max(local_vol, 1e-4) );

	m_localVols[t][k] = local_vol;

// compute bumped local vol

	vol		+= m_volbump;

/*	
	m_pBumpHelper->getBumpDerivatives( bres, rk, T );
	dsdk	+= bres[0];
	d2sd2k	+= bres[1];
	vol		+= bres[2];
	dsdt	+= bres[3];
*/


	SsqT = vol*sqT ;
	d	= -log(rk)/SsqT + 0.5*SsqT ;
	num = vol * ( vol + 2.* T * dsdt ) ;
	denom = 1.+ rk * d * dsdk * sqT ;
	denom *= denom ;
	temp = rk*rk *vol*T * ( d2sd2k - d*sqT * dsdk*dsdk );
	denom += temp;

	local_vol =  num/denom ;
	local_vol = sqrt( std::max(local_vol, 1e-4) );

	m_bLocalVols[t][k] = local_vol;
}




/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: parallelVegaBump
**	Returns: 
**	Action : 
**           
****************************************************************/


void DupireLocalVol::initQuantoDrift()
{
	m_quantoDrift.resize(m_nTime,0.0);

	std::string	payCurr = m_pAsset->GetPayZeroCurve(true)->GetName() ;
	std::string	natCurr = m_pAsset->GetNaturalZeroCurve(true)->GetName() ;
	std::string	compoFX = m_pAsset->GetFxCompositeAsset()->GetName() ;


	if( compoFX != std::string( natCurr + ".USD" ) ){
		throw "Composite case not supported in Local Vol";
	}

	if( natCurr == payCurr ){
		return;	
	}// no quanto adjustment


	if( payCurr == "USD" || natCurr == "USD" )
	{
		int SpotOrForward = 1;
		BidAskMidEnum middle = Middle;
		MlEqAssetHandle quanto_fx ;
		double corSX = 0.;

		if( payCurr == "USD" )
		{		
			quanto_fx = m_pAsset->GetFxNaturalAsset() ;

			std::string assetName = m_pAsset->GetName() ;
			std::string fxName = natCurr + ".USD" ;
			corSX = m_pAsset->GetCorrelationMatrix()->GetCorrelation( assetName, fxName );
		}
		else
		{
			quanto_fx = m_pAsset->GetFxPayAsset() ;
		
			std::string assetName = m_pAsset->GetName() ;
			std::string fxName = payCurr + ".USD" ;
			corSX = - m_pAsset->GetCorrelationMatrix()->GetCorrelation( assetName, fxName );
		}

		long start = m_nToday;
		for(int t=1; t<m_nTime; ++t)
		{
			long date = m_Dates[t];
			double volX = quanto_fx->getNaturalATMVol(start,date,SpotOrForward,middle);
			m_quantoDrift[t-1] = corSX * volX ;

			start = date;
		}
	}
	else	// we could have had only this case...
	{
		std::string assetName = m_pAsset->GetName() ;
		std::string fx1Name = natCurr + ".USD" ;
		std::string fx2Name = payCurr + ".USD" ;
		double corSX1 = m_pAsset->GetCorrelationMatrix()->GetCorrelation( assetName, fx1Name );
		double corSX2 = m_pAsset->GetCorrelationMatrix()->GetCorrelation( assetName, fx2Name );

		int SpotOrForward = 1;
		BidAskMidEnum middle = Middle;
		const MlEqAssetHandle quanto_fx1 = m_pAsset->GetFxNaturalAsset() ;
		const MlEqAssetHandle quanto_fx2 = m_pAsset->GetFxPayAsset() ;

		long start = m_nToday;
		for(int t=1; t<m_nTime; ++t)
		{
			long date = m_Dates[t];
			double volX1 = quanto_fx1->getNaturalATMVol(start,date,SpotOrForward,middle);
			double volX2 = quanto_fx2->getNaturalATMVol(start,date,SpotOrForward,middle);
		
			m_quantoDrift[t-1] = corSX1 * volX1 - corSX2 * volX2 ;

			start = date;
		}
		
	}

	m_quantoDrift[m_nTime-1] = m_quantoDrift[m_nTime-2];
}



/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: parallelVegaBump
**	Returns: 
**	Action : 
**           
****************************************************************/


void DupireLocalVol::parallelVegaBump()	// swap local vol and bumped local vol
{
	GVector<CVector> tmp = m_bLocalVols;
	m_bLocalVols = m_localVols;
	m_localVols = tmp;

	m_vbumped = !m_vbumped ;
}


/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: bumpLocalDrift
**	Returns: 
**	Action : 
**           
****************************************************************/


void DupireLocalVol::bumpLocalDrift( double driftbump )
{
	for ( int idate = 0 ; idate < m_nTime; idate++ ){
		m_loanspread[idate] += driftbump;
	}
}


/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: createLocalVolGrid
**	Returns: 
**	Action : 
**           
****************************************************************/

void DupireLocalVol::createLocalVolGrid()
{
	if(!!m_pLocalVolGrid)
	{
		createInterpolatedLocalVolGrid();	
	}
	else	// do the usual stuff...
	{
		m_localVols.resize(m_nTime);
		m_bLocalVols.resize(m_nTime);

		for(int t=0; t<m_nTime; t++)
		{
			m_localVols[t].resize(m_nSpace);
			m_bLocalVols[t].resize(m_nSpace);

			for(int k=0; k<m_nSpace; k++)
			{
				double rstrike = reducedStrike( t, m_fixedStrikes[k] );
				computeLocalVol(t,k, rstrike);
			}
		}
	}

}


/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: getLocalVolOnGrid
**	Returns: 
**	Action : 
**           
****************************************************************/

void DupireLocalVol::getLocalVolOnGrid( int slice, double low, double up, int nx, CVector& grid, int gridType )
{
	grid.resize(nx);

	CVector tmpFixed(nx,0.);
	if( gridType == 1100)
	{
		double dx = log(up/low) / double(nx-1.) ;
		for(int i=0; i<nx; ++i)
			tmpFixed[i] = low * exp( i*dx );
	}
	else
	{
		double dx = (up-low) / double(nx-1.) ;
		for(int i=0; i<nx; ++i)
			tmpFixed[i] = m_spot * exp(low + i*dx) ;
	}


	for(int k=0; k<nx; k++)
	{
		double y = yStrike( slice, tmpFixed[k] );
		grid[k]  =	m_pLocalVolGrid->getLocalVol( slice, y );	
	}
}


/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: bumpLocalVol
**	Returns: 
**	Action : 
**           
****************************************************************/
	

void DupireLocalVol::bumpLocalVol( int t, int k, double bump )
{
	m_localVols[t][k] += bump;
	double lv = getlocalVol(t,k);
	setlocalVol(lv+bump,t,k);
}

/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: getlocalVol
**	Returns: 
**	Action : 
**           
****************************************************************/

double	DupireLocalVol::getlocalVol(int idate,int ixspace)
{
		double vol = m_localVols[idate][ixspace];
		return vol;
}

/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: setlocalVol
**	Returns: 
**	Action : 
**           
****************************************************************/

void DupireLocalVol::setlocalVol(double value,int idate,int ixspace)
{
	m_localVols[idate][ixspace] = value;

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

	void initialize(GVector<long>& calibDates,GVector<MlEqAnalyticCurveWithTanhWingEdgeHandle>& lv);


  	void reinitialize(int it,MlEqAnalyticCurveWithTanhWingEdgeHandle& lv);
	void reinitialize(int it,CVector& lvols);


};
*/









/****************************************************************
**	Class  : DupireLocalVol 
**	Routine: initialize
**	Returns: 
**	Action : 
**           
****************************************************************/


void CLocalVolHelper::initialize( DupireLocalVolHandle lvh )
{
	m_nSlices	= lvh->m_nTime;
	m_nToday	= lvh->m_nToday ;
	m_forwards	= lvh->m_forwards;
	m_spot		= lvh->m_spot;

	int			addTanhWings = 1;	// as in mc model info
	double		cL = 0.01;	
	double		cR = 3.;			
	int			yPower = 2;	
	
	m_interpolator = new MlEqMonotonicSplineInterpolator( 	lvh->m_yStrikes,
															lvh->m_localVols,
															addTanhWings, cL, cR, yPower );

	m_binterpolator = new MlEqMonotonicSplineInterpolator( lvh->m_yStrikes,
															lvh->m_bLocalVols,
															addTanhWings, cL, cR, yPower );
}

/****************************************************************
**	Class  : CLocalVolHelper 
**	Routine: getLocalVol
**	Returns: 
**	Action : 
**           
****************************************************************/


double CLocalVolHelper::getLocalVol( int slice, double y )
{
	if( slice >= m_nSlices || slice < 0)	throw "Wrong date for local vol calculation...";

	double local_vol = m_interpolator->getValue( y, slice );
	return local_vol;
}

/****************************************************************
**	Class  : CLocalVolHelper 
**	Routine: getBumpedLocalVol
**	Returns: 
**	Action : 
**           
****************************************************************/

double CLocalVolHelper::getBumpedLocalVol( int slice, double y )
{
	if( slice >= m_nSlices || slice < 0)	throw "Wrong date for local vol calculation...";

	double local_vol = m_binterpolator->getValue( y, slice );
	return local_vol;
}




void CBumpVolHelper::initialize( MlEqAssetHandle pAsset, long dateMaturity, const CVector& mktForwards )
{

	CVector m_times;
	CVector m_rStrikes;
	GVector<CVector> m_bumpGrid;

	CVector res ;

// first set a few dates / strikes

	double dt = 10. ; // 30 days...

	MlEqConstDateHandle hDate = pAsset->GetDateHandle() ;
	long nToday = hDate->GetDate() ;

	long date = nToday ;
	std::vector<long> dates;
	dates.push_back( date);

	while( date < dateMaturity )
	{
		date = long( date + dt );
		dates.push_back( date );
	}
	int nt = dates.size();
	m_times.resize( nt );


	double dk = 0.01;	// 10% forward
	double kmin = 0.1, kmax = 2.0 ;

	int nk = int( (kmax - kmin) / dk ) + 1;

	m_rStrikes.resize(nk);

	for(int k=0; k<nk; ++k){
		m_rStrikes[k] = kmin + k * dk ;
	}


	m_bumpGrid.resize( nt );

	MlEqVolatilityStructureHandle hVol = pAsset->GetVolatilityStructure() ;

	for(int t=0; t<nt; ++t)
	{
		date = dates[t] ;
		double fwd = pAsset->GetForward( date, false );
		m_times[t] = hDate->GetYearFraction(date) ;

		m_bumpGrid[t].resize( nk );

		for(int k=0; k<nk; ++k)
		{
			MlEqForwardBasedStrike fwdStrike( m_rStrikes[k] , fwd, date, hDate );

			double bvol = hVol->getVol(fwdStrike, DATE(date) );

			hVol->getVolDeriv( res, fwdStrike, DATE(date), mktForwards );
			double vol = res[2];

			double tst = bvol/vol - 1;

			m_bumpGrid[t][k] = bvol - vol ;
		}
	}

	// now initialize interpolator....

	GVector<CVector> xVals(1);
	xVals[0] = m_rStrikes ;

	int			addTanhWings = 1;	// as in mc model info
	double		cL = 0.01;	
	double		cR = 3.;			
	int			yPower = 1;	

	MlEqInterpolatorHandle ShiftMatrix = new MlEqCubicSplineInterpolator(xVals, m_bumpGrid, addTanhWings, cL, cR, yPower );

	GVector<CVector> TIMES(1);
	TIMES[0] = m_times ;
//	MlEqInterpolatorHandle MaturityInterpolator = new MlEqCubicSplineInterpolator(TIMES, TIMES, addTanhWings, cL, cR, yPower );
	MlEqInterpolatorHandle MaturityInterpolator = new MlEqInterpolator(TIMES, TIMES);


	m_scenarioBump = new MlEq2DInterpolator(MaturityInterpolator, ShiftMatrix);		
}

void	CBumpVolHelper::getBumpDerivatives( CVector& bderivs, double reducedStrike, double yearFraction )
{
	const double eps = 1e-6;

	double dvol		= m_scenarioBump->getValue( yearFraction, reducedStrike );
	double dvol_	= m_scenarioBump->getValue( yearFraction, reducedStrike + eps );
	double _dvol	= m_scenarioBump->getValue( yearFraction, reducedStrike - eps );
	double dvol_t	= m_scenarioBump->getValue( yearFraction + eps, reducedStrike );

	bderivs.resize(4);

	bderivs[0]	=	(dvol_-_dvol) / (2.*eps) ;
	bderivs[1]	=	(dvol_+_dvol - 2.*dvol) / (eps*eps);
	bderivs[2]	=	dvol ;

	bderivs[3]	=	(dvol_t - dvol) / eps ;
}