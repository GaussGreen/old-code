

#include "stdafx.h"
#include "MlEqObjects.h"

#include "MonteCarlo.h"
#include <math.h>

#include <assert.h>
#include "SciSobol.h"
#include "hpsort.h"
#include "threemoments.h"
//#include "mleqshortmodels.h"
#include "localvol.h"
#include "heston.h"
#include "MlEqPde.h"



void getCorrelationMatrix( const std::vector<MlEqAssetHandle>& assets, CMatrix& correl )
{
	int nasset = assets.size();
	correl.resize(nasset, nasset);

	for(int iasset=0; iasset<nasset; ++iasset)
	{
		MlEqCorrelationMatrixHandle pCorrel = assets[iasset]->GetCorrelationMatrix();
		correl[iasset][iasset] = 1.;

		for(int jasset=iasset+1; jasset<nasset; ++jasset){
			double fcorrel = pCorrel->GetCorrelation(assets[iasset]->GetName(), assets[jasset]->GetName());	
			correl[iasset][jasset] = correl[jasset][iasset] = fcorrel ;
		}
	}
}





/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: initialize
**	Returns : nothing
**	Comment : initializes forward skew MonteCarlo
****************************************************************/


void CForwardSkewMC::multiplyForwardVols(int idate,double factor)
{

	for ( int i = 0 ; i < m_vols.cols(); i++ )
	{
		m_vols[idate][i] *= factor;
	}

	m_VolAtm[idate] *= factor;

}


/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: initialize
**	Returns : nothing
**	Comment : initializes forward skew MonteCarlo
****************************************************************/


void CForwardSkewMC::Initialize(double spot,const vector< long> & simdates,const CVector& fwds,const CVector& discounts,const CMatrix& vols,const CVector& beta,const CVector& volvol,const CVector& meanRevRate,const CVector& meanRevLevel,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes,CMatrix & controlVariatePrices, MlEq2DInterpolatorHandle pAccessVolInterp)
{

	m_iasset = 0;
	int m_nDates = volvol.getsize();

	if ( ( meanRevRate.getsize() != m_nDates ) || ( meanRevLevel.getsize() != m_nDates )
		|| ( beta.getsize() != m_nDates ) )
	{
		throw("inconsitent model market data entered");
	}


	CMatrix modelMktParams(m_nDates,4);

	for (int n = 0 ; n < m_nDates; n++ )
	{	
		modelMktParams[n][0] = beta[n]				;
		modelMktParams[n][2] = meanRevRate[n];
		modelMktParams[n][3] = meanRevLevel[n];
		modelMktParams[n][1] = volvol[n]		;
	}	

	Initialize(spot,simdates,fwds,discounts,vols,modelMktParams,modelParameters,calibTimeIndex,calibVols,pCalibStrikes,controlVariatePrices, pAccessVolInterp);
}

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: initialize
**	Returns : nothing
**	Comment : initializes forward skew MonteCarlo
****************************************************************/


void CForwardSkewMC::Initialize(double spot,const vector<long>& simdates,const CVector& fwds,const CVector& discounts,const CMatrix& vols,const CMatrix& modelMarketParameters,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes,CMatrix & controlVariatePrices, MlEq2DInterpolatorHandle pAccessVolInterp)
{

	CVector times;

	int i0=0;
	if ( m_hDate->GetYearFraction(simdates[0]) > 1e-3 ){
		times.resize(simdates.size());}
	else{
		i0=1;
		times.resize(simdates.size()-1);}

	int i;
	for ( i = 0 ; i < times.getsize(); i++ ){	
		times[i] = m_hDate->GetYearFraction(simdates[i+i0]);
	}

	int n,m;
	
//  retrieve basic market data


	if ( ( fwds.getsize() != discounts.getsize() ) || 
		 (fwds.getsize() != vols.rows()-1 ) || 
		 ( fwds.getsize() != modelMarketParameters.rows()) || 
		 ( fwds.getsize() != times.getsize() ))
		throw("intputs are inconsistent");

	m_nDates	= times.getsize();
	m_discount  = discounts;
	m_nAssets = 1;
	m_initialSpot = spot;




	if ( m_nDates == 0 )
	{
		m_nPaths = 1;
		if ( !m_randomsAreExternal ){
			m_path_array.resize(m_nPaths,1,m_nDates+1);
			m_pPath_array	= &m_path_array;
			m_iasset		= 0;
		}

		for ( n = 0; n < m_nPaths; n++ ){
			(*m_pPath_array)[n][m_iasset][0] = m_initialSpot;
		}

		m_discountPrepended.resize(1);
		m_discountPrepended[0] = 1.0;

		return;
	}

	m_fwds.resize(fwds.getsize(),1);
	for ( int i = 0 ; i < fwds.getsize(); i++ ){
		m_fwds[i][0] = fwds[i];
	}

	m_discountPrepended.resize(m_discount.getsize()+1);
	m_discountPrepended[0] = 1.0;
	for ( i = 0 ; i < m_discount.getsize(); i++ ){
		m_discountPrepended[i+1] = m_discount[i];
	}

	m_normalizedStrikes.resize(m_nDates,vols.cols());

	m_termFwds.resize(fwds.getsize()+1,1);
	m_termFwds[0][0] = m_initialSpot;
	for ( n = 0 ; n < fwds.getsize(); n++ ){
		m_termFwds[n+1][0] = fwds[n];
	}

	// m_fwds is the forwardforward

	m_fwds.resize(m_nDates,1);
	m_fwds[0][0] = fwds[0]/spot;
	for ( n = 1; n < fwds.getsize(); n++ ){
		m_fwds[n][0] = fwds[n]/fwds[n-1];
	}


	for ( n = 0 ; n < m_nDates; n++ ){
		for ( m = 0 ; m < vols.cols(); m++ ){
			m_normalizedStrikes[n][m] = vols[0][m];
		}
	}

	m_vols.resize(m_nDates,vols.cols());
	for ( n = 0 ; n < m_nDates; n++ )
	{
		for ( m = 0 ; m < vols.cols(); m++ ){
			m_vols[n][m] = vols[n+1][m];
		}
	}

	m_VolAtm.resize(m_nDates);
	for ( n = 0 ; n < m_nDates; n++ ){
		m_VolAtm[n] = MlEqMaths::linearInterp(m_normalizedStrikes[n],m_vols[n],0.0,0);
	}

	
	m_logFwdFactors.resize(m_fwds.rows());

	for ( n = 0 ; n < m_fwds.rows();n++){
		m_logFwdFactors[n] = log(m_fwds[n][0]);
	}

	m_beta.resize(m_nDates);
	m_meanReversionRate.resize(m_nDates);
	m_meanReversionLevel.resize(m_nDates);
	m_volvol.resize(m_nDates);
	
	for ( n = 0 ; n < m_nDates; n++ )
	{	
		m_beta[n]				= modelMarketParameters[n][0];
		m_meanReversionRate[n]	= modelMarketParameters[n][2];
		m_meanReversionLevel[n]	= modelMarketParameters[n][3];
		m_volvol[n]				= modelMarketParameters[n][1];
	}	

	m_dt.resize(m_nDates);
	m_sqrtdt.resize(m_nDates);
	
	for ( n = 0 ; n < m_nDates; n++ )
	{
		m_dt[n] = times[n];

		if ( n ){
			m_dt[n] -= times[n-1];
			assert(m_dt[n] > 0.0 );
		}
	
		m_sqrtdt[n] = sqrt(m_dt[n]);
	}

	m_t.resize(times.getsize()+1);
	m_t[0] = 0.0;
	for ( i = 0 ; i < times.getsize();i++ ){
		m_t[i+1] = times[i];
	}

	
// start setting up forwardSkewCurves
	
	m_ForwardVolcurve.resize(m_nDates);
	m_bareForwardVolcurve.resize(m_nDates);
	
	double cL			= modelParameters[0];
	double cR			= modelParameters[1];
	double addTanhWing	= modelParameters[2];
	double yPower		= modelParameters[3];
	m_seed				= modelParameters[4];
	m_nPaths			= modelParameters[5];
	int calibflag		= modelParameters[6];
	m_numberVolStates	= modelParameters[7];
	m_localVolFlag		= modelParameters[8];
	m_saveVolofVolInfo  = modelParameters[9];
	m_numberGridpoints  = modelParameters[10];
	m_saveVolGrid		= modelParameters[11];
	m_randomNumberFlag	= modelParameters[12];
	m_controlVariate	= modelParameters[13];
	m_globalCalib		= modelParameters[14];
	m_useHermites		= modelParameters[15];

	m_controlVariateCalcFlag = m_controlVariate;
	m_calibflag				 = calibflag;

	if ( m_saveVolofVolInfo ) 
	{
		m_averageBeta.resize(m_nDates);
		m_termvolvol.resize(m_nDates);
		m_spotvolvol.resize(m_nDates);
		m_prevVol.resize(m_nPaths);
		m_saveDiscreteVols.resize(m_nDates,m_numberVolStates);
	}

	for ( n = 0 ; n < m_nDates; n++ )
	{
		
		m_ForwardVolcurve[n].initialize( 
					m_normalizedStrikes[n],
					m_vols[n],
					addTanhWing,
					cL,				
					cR,		
					yPower);
	}


	m_bareForwardVolcurve = m_ForwardVolcurve;

	if ( !m_randomsAreExternal ){
		m_path_array.resize(m_nPaths,1,m_nDates+1);
		m_pPath_array	= &m_path_array;
		m_iasset		= 0;
	}

	for ( n = 0; n < m_nPaths; n++ ){
		(*m_pPath_array)[n][m_iasset][0] = m_initialSpot;
	}


	if ( m_saveVolGrid == 1 ){
		m_volMap.resize(m_nDates,m_nPaths);
		m_discreteVolMatrix.resize(m_nDates,m_numberVolStates);
	}
	else if (  m_saveVolGrid == 2 ){
		m_continuousVolMatrix.resize(m_nDates,m_nPaths);
	}

	m_volmultiplyer = 1.0;

	m_calibTimeIndex.resize(calibTimeIndex.rows());
	for ( i = 0 ; i < calibTimeIndex.rows(); i++ ){
		m_calibTimeIndex[i] = calibTimeIndex[i][0];
	}
	m_calibVols		 = calibVols;
	m_pCalibStrikes	 = pCalibStrikes;

	int fitsize = 2;
	CVector accuracy(fitsize);

	accuracy[0] = 0.0005;
	accuracy[1] = 0.0003;

	m_accuracy  = accuracy;

	int k;
	int ncols = calibTimeIndex.cols(); 
	if ( calibTimeIndex.rows() )
	{
		m_calibrationInfo.resize(calibTimeIndex.rows(),4);

		for (int i = 0 ; i < calibTimeIndex.rows(); i++ )
		{
			if ( ncols>1 )
			{
				for ( k = 0 ; k < ncols-1; k++ ){
					m_calibrationInfo[i][k] = calibTimeIndex[i][k+1];
				}
			}
			else
			{
				// define defaults
				m_calibrationInfo[i][0] = 1;// this is calibration to level
				m_calibrationInfo[i][1] = 1;// this is calibration to beta
				m_calibrationInfo[i][2] = 0;// this is calibration to volvol
				m_calibrationInfo[i][3] = 0;// this is calibration to meanreversion

			}
		}
	}


	m_controlVariate = m_controlVariate;

	if ( m_controlVariate )
	{
		m_controlVariateCurrentSpot.resize(m_nPaths);
		for ( n = 0; n < m_nPaths; n++ )
		{
			m_controlVariateCurrentSpot[n] = m_initialSpot;
		}

		m_controlVariateVols.resize(m_VolAtm.getsize());
		m_controlVariateVols[0] = m_VolAtm[0];
		for ( n = 1 ; n < m_controlVariateVols.getsize(); n++ )
		{
m_controlVariateVols[n] = sqrt((pow(m_controlVariateVols[n-1],2.0)*m_t[n]+pow(m_VolAtm[n],2.0)*m_dt[n])/m_t[n+1] );
		}

	}

	m_controlVariatePrices = controlVariatePrices;


	if ( m_controlVariatePrices.rows() )
	{
		m_controlVariateVols.resize(m_VolAtm.getsize());
		m_controlVariateVols[0] = m_VolAtm[0];
		for ( n = 1 ; n < m_controlVariateVols.getsize(); n++ ){
m_controlVariateVols[n] = sqrt((pow(m_controlVariateVols[n-1],2.0)*m_t[n]+pow(m_VolAtm[n]*1.0,2.0)*m_dt[n])/m_t[n+1] );
		}

	}
	
	m_rev0			= -1;	
	m_rev1			= -1;
	m_xslope0		= -1;
	m_xslope1		= -1;
	m_mrUpdateOnly	= -1;
	
	long idum = m_seed;
	
	if ( !!pAccessVolInterp )
	{	
		m_excessVolVolInterp.resize(m_nDates);
		
		for ( int i = 0 ; i < m_nDates; i++ )
		{	
			
			MlEqConstInterpolatorHandle pInterp = pAccessVolInterp->getXInterpolator(i);
			const CVector& xData = pInterp->getXData(i);
			CVector yData(xData.getsize());
			for ( int j = 0 ; j < xData.getsize(); j++ ){				
				yData[j] = pAccessVolInterp->getValue(xData[j],m_dt[i]);
			}

			m_excessVolVolInterp[i]	=	new MlEqInterpolator(xData,yData);
		}	
			
	}		

	m_numberOfFactors = 2;
	initializeNumberOfFactors();
	m_centerFwds = false;
}

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: initialize
**	Returns : nothing
**	Comment : initializes forward skew MonteCarlo
****************************************************************/

void CForwardSkewMC::initializeNumberOfFactors()
{
	m_numberOfFactors = 2;;
}

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: initialize
**	Returns : nothing
**	Comment : initializes forward skew MonteCarlo
****************************************************************/

void CForwardSkewMC::Initialize(MlEqAsset& asset,product& deriv,MlEqStochBetaVolatilityStructure& betaVolStruct,CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes)
{
	const std::vector<long> additionalSimDates;
	Initialize(asset,deriv,additionalSimDates,betaVolStruct,fixings,modelParameters, calibTimeIndex,calibVols,pCalibStrikes);
}

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: initialize
**	Returns : nothing
**	Comment : initializes forward skew MonteCarlo
****************************************************************/

void CForwardSkewMC::Initialize(MlEqAsset& asset,product& deriv,const std::vector<long>& additionalSimDates,MlEqStochBetaVolatilityStructure& betaVolStruct,CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes)
{	
	initializeMcHelper(deriv, fixings, additionalSimDates);
	const std::vector<long>& mcDates	= m_mMCHelper.GetMCDates();;
	Initialize(asset,mcDates,betaVolStruct,modelParameters, calibTimeIndex,calibVols,pCalibStrikes);
}


/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: initialize
**	Returns : nothing
**	Comment : initializes forward skew MonteCarlo
****************************************************************/

void CForwardSkewMC::Initialize(MlEqAsset& asset,product& deriv,const std::vector<long>& additionalSimDates, CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes)
{
	MlEqVolatilityStructureHandle			pVol =   asset.GetVolatilityStructure();
	MlEqStochBetaVolatilityStructureHandle	pBetaVol(dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*pVol));

	if (!pBetaVol){
		throw(" asset handle must contain stochastic beta volatility structure for monte carlo calculation");
	}

	Initialize(asset,deriv,additionalSimDates, *pBetaVol,fixings,modelParameters, calibTimeIndex,calibVols,pCalibStrikes);
}


/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: initialize
**	Returns : nothing
**	Comment : initializes forward skew MonteCarlo
****************************************************************/

void CForwardSkewMC::Initialize(MlEqAsset& asset,product& deriv, CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes)
{
	const std::vector<long> additionalSimDates;
	Initialize(asset,deriv,additionalSimDates,fixings,modelParameters,calibTimeIndex,calibVols,pCalibStrikes);
}

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: initialize
**	Returns : nothing
**	Comment : initializes forward skew MonteCarlo
****************************************************************/

void CForwardSkewMC::Initialize(MlEqAsset& asset,const std::vector<long>& mcDates,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes)
{
	MlEqVolatilityStructureHandle			pVol =   asset.GetVolatilityStructure();
	MlEqStochBetaVolatilityStructureHandle	pBetaVol(dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*pVol));

	if (!pBetaVol){
		throw(" asset handle must contain stochastic beta volatility structure for monte carlo calculation");
	}

	long nToday = mcDates[0];//sos
	CMatrix fixings;
	m_mMCHelper.initialize(mcDates,fixings,nToday);

	Initialize(asset,mcDates,*pBetaVol,modelParameters,calibTimeIndex,calibVols,pCalibStrikes);
}

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: CForwardSkewMC
**	Returns : nothing
**	Comment : constructor
****************************************************************/

CForwardSkewMC::CForwardSkewMC(MlEqConstDateHandle hDate,MlEqAsset& asset,product& deriv,MlEqStochBetaVolatilityStructure& betaVolStruct,CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes) : MlEqMonteCarlo(hDate)
{
	Initialize(asset,deriv,betaVolStruct,fixings,modelParameters, calibTimeIndex,calibVols,pCalibStrikes);
}

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: CForwardSkewMC
**	Returns : nothing
**	Comment : constructor
****************************************************************/

CForwardSkewMC::CForwardSkewMC(MlEqConstDateHandle hDate,MlEqAsset& asset,product& deriv,CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes) : MlEqMonteCarlo(hDate)
{
	Initialize(asset,deriv,fixings,modelParameters,calibTimeIndex,calibVols,pCalibStrikes);
}

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: CForwardSkewMC
**	Returns : nothing
**	Comment : constructor
****************************************************************/

CForwardSkewMC::CForwardSkewMC(MlEqConstDateHandle hDate,MlEqAsset& asset,const std::vector<long>& mcDates,MlEqStochBetaVolatilityStructure& betaVolStruct,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes) : MlEqMonteCarlo(hDate)
{
	Initialize(asset,mcDates,betaVolStruct,modelParameters, calibTimeIndex,calibVols,pCalibStrikes);
}

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: CForwardSkewMC
**	Returns : nothing
**	Comment : constructor
****************************************************************/

CForwardSkewMC::CForwardSkewMC(MlEqConstDateHandle hDate,MlEqAsset& asset,const std::vector<long>& mcDates,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes) : MlEqMonteCarlo(hDate)
{
	Initialize(asset,mcDates,modelParameters,calibTimeIndex,calibVols,pCalibStrikes);
}


/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: initialize
**	Returns : nothing
**	Comment : initializes forward skew MonteCarlo
****************************************************************/

void CForwardSkewMC::Initialize(MlEqAsset& asset,const std::vector<long>& mcDates,MlEqStochBetaVolatilityStructure& betaVolStruct,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes)
{
	
	m_nAssets = 1;
	m_assets.resize(m_nAssets);
	m_assets[0] = &asset;

	
	CVector ctimes;
	CVector discounts;
	
	MlEqZeroCurveHandle	payZeroCurve = asset.GetPayZeroCurve(true);
	
	int i0=0;
	if ( m_hDate->GetYearFraction(mcDates[0]) > 1e-3 ){
		ctimes.resize(mcDates.size());
		discounts.resize(mcDates.size());
	}
	else{
		i0=1;
		ctimes.resize(mcDates.size()-1);
		discounts.resize(mcDates.size()-1);
	}
	
	int i;
	for ( i = 0 ; i < ctimes.getsize(); i++ )
	{
		ctimes[i]		=	m_hDate->GetYearFraction(mcDates[i+i0]);
		discounts[i]	=	payZeroCurve->GetDiscountFactor(m_hDate->GetDate(),mcDates[i+i0]);
	}
	
	int ndates = ctimes.getsize();
	m_nDates = ndates;
	CVector fwds(ndates);
	
	double spot = asset.GetSpot(m_hDate->GetDate());
	for ( i = 0 ; i < ctimes.getsize(); i++ ){
		fwds[i] = asset.GetQuantoForward(m_hDate->GetDate(),mcDates[i+i0], false);
	}
	
	CVector normalizedStrikes(7);
	
	normalizedStrikes[0] = -2.0;
	normalizedStrikes[1] = -1.5;
	normalizedStrikes[2] = -1.0;
	normalizedStrikes[3] = -0.5;
	normalizedStrikes[4] =  0.0;
	normalizedStrikes[5] =  0.5;
	normalizedStrikes[6] =  1.0;

  
	CMatrix fwdVols(ndates+1,normalizedStrikes.getsize());

	for ( i = 0 ; i < normalizedStrikes.getsize(); i++ ){	
		fwdVols[0][i] = normalizedStrikes[i];
	}

	MlEqVolatilityStructureHandle   pVol =   asset.GetVolatilityStructure();

	double mat,fwdParam;
	CVector strikeVol(ndates);
	for ( int n = 0 ; n < ndates; n++ )
	{
		mat = m_hDate->GetYearFraction(mcDates[n+1]) ;
		fwdParam = fwds[n];
		mat -= m_hDate->GetYearFraction(mcDates[n]);

		double fwdShort,fwdNaturalShort;
		double fwdLong,fwdNaturalLong;

		if (  n == 0 )
		{ 
			fwdParam /= spot;
			fwdShort = spot;
			fwdNaturalShort = asset.GetNaturalSpot(mcDates[0]); 
		}
		else
		{
			fwdParam /= fwds[n-1];
			fwdShort   = fwds[n-1];
			fwdNaturalShort = asset.GetNaturalForward(mcDates[0], mcDates[n], false); 
		}

		fwdLong   = fwds[n];
		fwdNaturalLong = asset.GetNaturalForward(mcDates[0],mcDates[n+1], false);

		MlEqStrike shortStk(fwdShort);
		MlEqStrike longStk(fwdLong);
		MlEqStrike shortNaturalStk(fwdNaturalShort);
		MlEqStrike longNaturalStk(fwdNaturalLong);

		strikeVol[n] = pVol->getNaiveFutureVol(shortNaturalStk,longNaturalStk,mcDates[n],mcDates[n+1]); 

		for ( i = 0 ; i < normalizedStrikes.getsize(); i++ )
		{
			double vol = asset.GetCompositeVolatility(mcDates[n],(double)mcDates[n+1],strikeVol[n] );
			double fwdfwd = asset.GetNaturalForward(mcDates[n],mcDates[n+1], false);

			MlEqDateHandle startDate = new MlEqDate(mcDates[n],m_hDate->GetDayCountConvention());
			MlEqNormalizedStrike normStrike(normalizedStrikes[i],mcDates[n+1],vol,fwdfwd,true,startDate);

			fwdVols[n+1][i] = asset.GetCompositeVolatility(normStrike,mcDates[n],(double)mcDates[n+1], Middle);

		}
	}


	CMatrix modelMarketParameters(ndates,4);

	DATE startdate;
	for ( int n = 0 ; n < m_nDates; n++ )
	{	

		mat = m_hDate->GetYearFraction(mcDates[n+1]) ;
		fwdParam = fwds[n];
		mat -= m_hDate->GetYearFraction(mcDates[n]);

		if (  n == 0 ){ 
			fwdParam  /= spot;
		}
		else{
			fwdParam /= fwds[n-1];
		}
 
		startdate = mcDates[n];

		double fwdfwd = asset.GetNaturalForward(mcDates[n], mcDates[n+1], false);
		MlEqStrike atmStrike(fwdfwd);

		double shortFwd = asset.GetNaturalForward(mcDates[n], false);

		double vols = pVol->getFutureVol(atmStrike,mcDates[n],mcDates[n+1],shortFwd); 

		double vol = asset.GetCompositeVolatility(mcDates[n],(double)mcDates[n+1],vols, Middle);

		modelMarketParameters[n][0] = betaVolStruct.getBeta(startdate,mcDates[n+1])*vol;
		modelMarketParameters[n][2] = betaVolStruct.getMeanReversion(startdate,mcDates[n+1]);
		modelMarketParameters[n][3] = strikeVol[n];
		modelMarketParameters[n][1] = betaVolStruct.getVolVol(startdate,mcDates[n+1]);		

	}	


//  start to create simulation dates

	CMatrix  controlVariatePrices;
	// deal with the following later	

	Initialize(spot,mcDates,fwds,discounts,fwdVols,modelMarketParameters,modelParameters,calibTimeIndex,calibVols,pCalibStrikes,controlVariatePrices,betaVolStruct.m_excessVolOfVol);

	m_useHermites = modelParameters[15];
	if ( !m_useHermites || !betaVolStruct.m_hermiteElasticity){
		return;
	}


	m_hermiteElasticity.resize(m_nDates);
	for ( int i = 0 ; i < m_nDates; i++ ){
		m_hermiteElasticity[i] = betaVolStruct.m_hermiteElasticity->getValue(m_dt[i]);
	}

//  start setting up hermite Skew dynamics here

	setupHermites(asset,mcDates,strikeVol );

}

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: initialize
**	Returns : nothing
**	Comment : initializes forward skew MonteCarlo
****************************************************************/

void CForwardSkewMC::setupHermites(MlEqAsset& asset,const std::vector<long>& mcDates,CVector& strikeVol )
{
	double slope;
	CVector targetSlopes(m_nDates);
	for ( int n = 1 ; n < m_nDates; n++ )
	{
		int idate = n;
		double nsUp = 0.01;
		double nsDo = -0.01;

		double mat = m_hDate->GetYearFraction(mcDates[n+1]) ;
		mat -= m_hDate->GetYearFraction(mcDates[n]);

		double vol = m_VolAtm[idate];
//		double fwdfwd = asset.GetNaturalForward(mcDates[0],mcDates[n+1]-mcDates[n]+mcDates[0], false);


		double fwdfwd = asset.GetNaturalForward(mcDates[n],mcDates[n+1], false);


		MlEqDateHandle hStart = new MlEqDate(mcDates[n],m_hDate->GetDayCountConvention());
		MlEqNormalizedStrike normStrikeUp(nsUp,mcDates[n+1],vol,fwdfwd,false,hStart);


		
		/*MlEqNormalizedStrike(double strike, 
			long nMaturity,
			double vol,
			double fwd,
				bool isFloating = true,
				MlEqDateHandle startDate= NULL,double volscale=1.0,double elasticity=1.0,
						 double maturityShift=0.0);*/



		

		slope = asset.GetCompositeVolatility(normStrikeUp,mcDates[n],(double)mcDates[n+1], Middle);


		MlEqNormalizedStrike normStrikeDo(nsDo,mcDates[n+1],vol,fwdfwd,false,hStart);

		slope -= asset.GetCompositeVolatility(normStrikeDo,mcDates[n],(double)mcDates[n+1], Middle);
		slope /= (nsUp-nsDo)*vol;

		targetSlopes[idate] = slope;


/*		int idate = n;
		double nsUp = 0.01;
		double nsDo = -0.01;

		double mat = m_hDate->GetYearFraction(mcDates[n+1]) ;
		mat -= m_hDate->GetYearFraction(mcDates[n]);
		double vol = asset.GetCompositeVolatility(mcDates[n],(double)mcDates[n+1],strikeVol[n] );
		double fwdfwd = asset.GetNaturalForward(mcDates[n],mcDates[n+1], false);

		double shift = 0.01;

		MlEqStrike StrikeUp(fwdfwd*(1.0+shift));
		slope = asset.GetCompositeVolatility(StrikeUp,mcDates[n],(double)mcDates[n+1],MlEqVolatilityStructure::MidVol);

		MlEqStrike Strike(fwdfwd);
		slope -= asset.GetCompositeVolatility(Strike,mcDates[n],(double)mcDates[n+1],MlEqVolatilityStructure::MidVol);

		double ns = log(1.0+shift)/(m_VolAtm[idate]*sqrt(mat));
		slope /= ns;

		targetSlopes[idate] = slope;
*/


	}



	setupHermites(m_VolAtm,targetSlopes);

}


/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: initialize
**	Returns : nothing
**	Comment : initializes forward skew MonteCarlo
****************************************************************/

void CForwardSkewMC::setupHermites(CVector& volatm,CVector& targetSlopes)
{

	m_hermiteSkew.resize(m_nDates);// first element not used
	m_currentHermiteSkews.resize(m_numberVolStates);

	CVector hermitCoeff(4);

	hermitCoeff[1] = -0.2;


	double eps = 0.05;
	for ( int idate = 1; idate < m_nDates; idate++ )
	{

		hermitCoeff[0] = m_VolAtm[idate];
		m_hermiteSkew[idate].initialize(m_dt[idate],m_fwds[idate][0],hermitCoeff,false);
		m_hermiteSkew[idate].adjustSlope(targetSlopes[idate],m_fwds[idate][0],eps);

	}
}


/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: generatePath
**	Returns: double
**	Action : generates MonteCarlo paths
****************************************************************/

void CForwardSkewMC::GeneratePath(int calibflag,long idum,CVector& accuracy)
{
	if ( m_nDates == 0 ){
		return;
	}

	if ( m_globalCalib )
	{
		long idum = m_seed;

		m_calibResults.resize(m_calibTimeIndex.getsize(),6);

		GlobalCalibrate(m_calibVols,m_pCalibStrikes,
						idum,accuracy);

	}
	else if ( m_calibTimeIndex.getsize() )
	{
		long idum = m_seed;

		int idateCalibIndex = 0;
		int idateCalib = m_calibTimeIndex[idateCalibIndex];
		int calibGlobalFlag = 1;

		m_calibResults.resize(m_calibTimeIndex.getsize(),6);
		calibrateTimeSlice(idateCalibIndex,m_calibVols[idateCalibIndex],m_pCalibStrikes[idateCalib-1],idum,accuracy,calibGlobalFlag,m_calibrationInfo);

		int maxiter=12, succeedType=3;

		for ( int i = 1; i < m_calibTimeIndex.getsize(); i++ )
		{
			calibGlobalFlag = 0;
			calibrateTimeSlice(i,m_calibVols[i],m_pCalibStrikes[m_calibTimeIndex[i]-1],idum,accuracy,calibGlobalFlag,m_calibrationInfo);
		}
		for ( int i = m_calibTimeIndex[m_calibTimeIndex.getsize()-1]; i < m_nDates; i++ )
		{
			generatePath(i,calibflag,idum);
		}

	}
	else
	{
		GeneratePath(calibflag);
	}


	m_pcurrentVols = NULL;
	m_pnewcurrentVols = NULL;
	m_psavecurrentVols = NULL;
	m_currentVols.resize(0);
	m_oneStepRandoms.resize(0);
	m_controlVariateCurrentSpot.resize(0);
	m_allRandoms.resize(0,0);

	GMatrix < CVector >		m_allMultiRandoms;//[ipath][iasset][m_numberOfFactors*idate]
	for ( int i = 0; i <  m_allMultiRandoms.rows(); i++ ){
		for ( int j = 0; j < m_allMultiRandoms.rows(); j++ ){
			m_allMultiRandoms[i][j].resize(0);
		}
	}

	m_allMultiRandoms.resize(0,0);

}

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: GetForwardDiscount
**	Returns : double
**	Comment : returns discount at idate
****************************************************************/


double MlEqMonteCarlo::GetBrownianBridgeHitProb(int ipath,int idate,double U,double L)
{
	throw "Function not supported for the Monte Carlo type '" + GetName() + "'";
}

std::string MlEqMonteCarlo::GetName(void) const
{
	std::string sz = typeid(*this).name();	
	return sz.substr(6); // To remove the 'Class ' prefix.
}
	
double MlEqMonteCarlo::GetDiscount(int ipath,int idate)
{
	if ( idate == 0 ){
		return 1.0;
	}

	if ( idate > m_nDates ){
		throw("error out of bound in GetForwardDiscount");
	}

	return	m_discount[idate-1];
}

double	MlEqMonteCarlo::GetBridgeVolatility(int ipath,int idate,double finalSpot,int iasset)
{
	throw("function not implemented");
	return -1.0;
}

double MlEqMonteCarlo::GetImpliedLogContractVol(int ipath,int idate,int iasset,int ngauss,double lowerStdev,double upperStdev)
{

	throw("function not implemented");
	return -1.0;

}

double MlEqMonteCarlo::GetDiscount(int ipath,int idate,int idateAsOf)
{
	throw("this was a bug; I am curious where this function was called previously");

	double disc = GetDiscount(ipath,idate);
	disc /= GetDiscount(ipath,idateAsOf);

	return disc;
}

double	MlEqMonteCarlo::GetDiscount(int ipath,int idate,int jdate,int idateAsOf)
{
	double z ;

	if ( idate > jdate || idateAsOf > idate ){
		throw("indexing problem in MlEqMonteCarlo::GetDiscount");
	}

	z = GetDiscount(ipath,jdate)/GetDiscount(ipath,idate);
	return z;
}

double	MlEqMonteCarlo::GetDiscountToDate(int ipath,long To,int idateAsOf)
{
	double disc = getAssets()[0]->GetPayZeroCurve(true)->GetDiscountFactor(GetMCDates()[idateAsOf],To);
	return disc;
}

double MlEqMonteCarlo::GetStochasticDiscount(int ipath,int idate)
{
	return GetDiscount(ipath,idate);
}

const CVector& MlEqMonteCarlo::GetDiscounts(int ipath)
{
	return m_discountPrepended;
}




/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: findAsset
**	Returns: void
**	Action : 
****************************************************************/

int	MlEqMonteCarlo::findAsset(MlEqAssetHandle & asset)
{
	for ( int i = 0 ; i < m_assets.size(); i++ ){

		if ( asset->GetName() == m_assets[i]->GetName() ){	
//		if (asset == m_assets[i]) return i;
			return i;}
	}

	throw "Asset '" + asset->GetName() + "' not found.";
}


/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: GetDt
**	Returns : double
**	Comment : returns time differences
****************************************************************/

double	CForwardSkewMC::GetDt(int idate)
{
	return m_dt[idate];
}





/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: GetDiscreteVol
**	Returns : double
**	Comment : returns discrete vols
****************************************************************/
	
double CForwardSkewMC::GetSimulatedVols(int ipath,int idate)
{	
	double vol;
	if ( m_saveVolGrid == 0)
	{
		throw("ForwardVolHandle has not been created with saveVolMatrix flag on");
	}
	else if ( m_saveVolGrid == 1)
	{
		vol = m_discreteVolMatrix[idate][m_volMap[idate][ipath]];
	}
	else if ( m_saveVolGrid == 2)
	{
		vol = m_continuousVolMatrix[idate][ipath];
	}
	else
	{
		throw("ForwardVolHandle has not been created with incorrect savevolflag");
	}

	return vol;

}

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: GetDiscreteVol
**	Returns : double
**	Comment : returns discrete vols
****************************************************************/

double	CForwardSkewMC::GetSimulatedVolSquareDt(int ipath,int iStart,int iEnd,int normalize)
{

	double val = 0.0;
	double t = 0.0;
	for ( int i = iStart; i <= iEnd; i++ )
	{
		val += pow(GetSimulatedVols(ipath,i),2.0)*m_dt[i];
		t += m_dt[i];
	}

	if ( normalize )
	{
		val /= t;
	}
	return val;
}

	

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: GetBrownianBridgeHitProb
**	Returns : double
**	Comment : returns Brownain Bridge hitprob
****************************************************************/
	
double CForwardSkewMC::GetBridgeVolatility(int ipath,int idate,double finalSpot,int iasset)
{
	if ( idate == m_nDates )
	{
		throw("CForwardSkewMC::GetBridgeVolatility ** incorrect date index entered **");
	}
	
	if ( m_saveVolGrid== 0 )
	{
		throw("CForwardSkewMC::GetBridgeVolatility ** ForwardVolHandle has not been created with saveVolGrid flag on **");
	}

	double newVolAtm		=	GetSimulatedVols(ipath,idate);
	double initialSpot		=	GetPathValue(ipath,idate);

	double strikeVol		=	GetBridgeVolatility( newVolAtm, idate, finalSpot,initialSpot, iasset);

	return strikeVol;
}	

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: GetBrownianBridgeHitProb
**	Returns : double
**	Comment : returns Brownain Bridge hitprob
****************************************************************/
	
double CForwardSkewMC::GetImpliedLogContractVol(int ipath,int idate,int iasset,int ngauss,double lowerStdev,double upperStdev)
{
	if ( idate-1 == m_nDates ){
		throw("CForwardSkewMC::GetImpliedLogContractVol ** incorrect date index entered **");
	}
	
	if ( m_saveVolGrid== 0 ){
		throw("CForwardSkewMC::GetImpliedLogContractVol ** ForwardVolHandle has not been created with saveVolGrid flag on **");
	}

	if ( m_logContractVols.getsize() == 0 ){
		m_logContractVols.resize(m_nDates);
	}

	if ( m_logContractVols[idate].getsize() > 0 )
	{
		// assumes that value has aleady been precomputed
		double res = m_logContractVols[idate][m_volMap[idate][ipath]];
		return res;
	}


	m_logContractVols[idate].resize(GetNumberOfVolStates(ipath,idate));

	if ( m_gaussPoint.getsize() == 0 )
	{

		double lower = lowerStdev;
		double upper = upperStdev;

		m_gaussWeight.resize(ngauss);
		m_gaussPoint.resize(ngauss);

		MlEqMaths::dGauleg(lower,upper,m_gaussPoint,m_gaussWeight,ngauss,true);			
	}

	if ( m_bareForwardVolcurve.size() == 0 ){
		throw("method not implemented for hermite skews yet");
	}


	double initialSpot		=	GetPathValue(ipath,idate);
	double bareVolAtm		=	m_bareForwardVolcurve[idate].getValue(0.0);	

	double vol,st,mat,excessVolOfVol;

	mat = m_dt[idate];

	for ( int istate = 0 ; istate < m_logContractVols[idate].getsize(); istate++ )
	{
		double newVolAtm		=	m_discreteVolMatrix[idate][istate];

		double res = 0.0;
		for ( int i = 0 ; i < m_gaussPoint.getsize(); i++ )
		{
			st = exp(m_gaussPoint[i]*newVolAtm*m_sqrtdt[idate]);
			excessVolOfVol	=	getExcessVolOfVol(idate,m_gaussPoint[i]);
			vol				=	getNewVolSkew(m_bareForwardVolcurve[idate],m_gaussPoint[i],newVolAtm,bareVolAtm,excessVolOfVol);

			if ( m_gaussPoint[i] < 0 ){
				res += m_gaussWeight[i]/pow(st,1.0)*Bs(1,vol,mat,st,1,-1);
			}
			else{
				res += m_gaussWeight[i]/pow(st,1.0)*Bs(1,vol,mat,st,1,1);
			}
		}


		res *= newVolAtm*sqrt(mat);
		res *= 2.0/mat;
		res = sqrt(res);

		m_logContractVols[idate][istate] = res;
	}

	double res =  m_logContractVols[idate][m_volMap[idate][ipath]];

	return res;

}	




/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: GetBrownianBridgeHitProb
**	Returns : double
**	Comment : returns Brownain Bridge hitprob
****************************************************************/
	
double CForwardSkewMC::GetBridgeVolatility(double newVolAtm,int idate,double finalSpot,double initialSpot,int iasset)
{	
	if ( idate == m_nDates )
	{
		throw("CForwardSkewMC::GetBridgeVolatility ** incorrect date index entered **");
	}
	
	if ( m_saveVolGrid== 0 )
	{
		throw("CForwardSkewMC::GetBridgeVolatility ** ForwardVolHandle has not been created with saveVolGrid flag on **");
	}

	double bareVolAtm		=	m_bareForwardVolcurve[idate].getValue(0.0);	
	double normStrike		=	log(finalSpot/(initialSpot*m_fwds[idate][0]))/(newVolAtm*m_sqrtdt[idate]);
	double excessVolOfVol	=	getExcessVolOfVol(idate,normStrike);
	double strikeVol = getNewVolSkew(m_bareForwardVolcurve[idate],normStrike,newVolAtm,bareVolAtm,excessVolOfVol);

	return strikeVol;

}	



/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: GetBrownianBridgeHitProb
**	Returns : double
**	Comment : returns Brownain Bridge hitprob
****************************************************************/
	
double CForwardSkewMC::AmericanDigital(CVector& prices,int ipath,int idate,double DownBarrier,double finalSpot,
									 int iasset,int method,double stressForwardBarrierVol,
									 double stressForwardBarrierVolSlope,int npoints,double blend)
{	

	double strike,vol,fwdfwd,spot,barrier,dt,xprice;

	if ( method <= 1 )
	{
		if ( m_saveVolGrid == 0){
			throw("ForwardVolHandle has not been created with saveVolMatrix flag on");
		}

		if ( prices.getsize() > 0 )
		{
			// assumes that value has aleady been precomputed

			double res = prices[m_volMap[idate][ipath]];
			return res;
		}
	}

	if ( method == 0 )
	{

			// compute all vols here
			prices.resize(GetNumberOfVolStates(ipath,idate));

			double m_VolatilityStrikePlace = stressForwardBarrierVolSlope;

			strike		=	m_VolatilityStrikePlace*DownBarrier*GetPathValue(ipath,idate)
							+(1-m_VolatilityStrikePlace)*GetPathValue(ipath,idate);

			for ( int istate = 0 ; istate < prices.getsize(); istate++ )
			{
				vol			=	GetBridgeVolatility(ipath,idate,strike);
				fwdfwd		=	GetForwardForward(idate,0);
				dt			=	Getdt(idate);
				spot		=	GetPathValue(ipath,idate);
				barrier		=	DownBarrier;
								
				double_knock_out_rebate(xprice,dt,1,fwdfwd,1,vol,blend,spot*50,barrier,1);
				prices[istate] = xprice;

			}

			double res = prices[m_volMap[idate][ipath]];
			return res;


	}
	else if ( method == 1 )
	{

		// compute all vols here
		prices.resize(GetNumberOfVolStates(ipath,idate));

		for ( int istate = 0 ; istate < prices.getsize(); istate++ )
		{
			double eps					=	1e-2;
			spot						=	GetPathValue(ipath,idate);
			fwdfwd						=	GetForwardForward(idate,0);
			strike						=	DownBarrier*spot;
			dt							=	Getdt(idate);

			double currentVol			=	m_discreteVolMatrix[idate][istate];

		
			double SpotBarrierVols			= GetBridgeVolatility((double)currentVol,idate,strike,spot,0);
			double shiftedSpotBarrierVols	= GetBridgeVolatility((double)currentVol,idate,strike*(1.0+eps),spot,0);
			
			double SpotBarrierVolSlope	=   (shiftedSpotBarrierVols-SpotBarrierVols)/(strike*eps)*spot;


			double FwdBarrierVols		=	stressForwardBarrierVol*SpotBarrierVols;

			double FwdBarrierVolsSlope	=   ( GetBridgeVolatility((double)currentVol,idate,spot*(1.0+eps),spot,0)-
											  GetBridgeVolatility((double)currentVol,idate,spot,spot,0) )/(spot*eps)*spot;


//			double FwdBarrierVolsSlope  =  stressForwardBarrierVolSlope*SpotBarrierVolSlope;

			FwdBarrierVolsSlope  *=  stressForwardBarrierVolSlope;



			CMatrix HitProb;
			double xspot = 1.0;
			double drift = log(fwdfwd)/dt;
			BootstrapStoppingTimes(
								xprice,HitProb,DownBarrier,xspot,dt,
								SpotBarrierVols,drift,0,
								SpotBarrierVolSlope,
								FwdBarrierVols,FwdBarrierVolsSlope,
								true,npoints,false
							);

/*
			BootstrapStoppingTimes2(
								xprice,HitProb,DownBarrier,xspot,dt,
								SpotBarrierVols,drift,0,
								SpotBarrierVolSlope,
								FwdBarrierVols,FwdBarrierVolsSlope,
								SpotBarrierVols, // vol barrier --> strike
								1.,				 // strike
								1.,			 	 // cp
								0.05,			 // rebate
								true,npoints,false
							);

*/			prices[istate] = xprice;
		}

		double res = prices[m_volMap[idate][ipath]];
		return res;

	}
	else
	{
		spot		=	GetPathValue(ipath,idate);
		double res = 
		GetBrownianBridgeHitProb(ipath,idate,50.0*spot,finalSpot);
		return  res;
	}


}



/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: GetBrownianBridgeHitProb
**	Returns : double
**	Comment : returns Brownain Bridge hitprob
****************************************************************/
	
double CForwardSkewMC::GetBrownianBridgeHitProb(int ipath,int idate,double U,double L)
{	

	double initialSpot		=	GetPathValue(ipath,idate);
	double finalSpot		=	GetPathValue(ipath,idate+1);	

	double prob;
/*	double currentVolAtm  = GetSimulatedVols(ipath,idate);

	double x = log(finalSpot/initialSpot)/currentVolAtm;

	double eps = 1e-3;
	double xdens = getCummulativeProb(log(finalSpot*(1+eps)/initialSpot),currentVolAtm,idate);
	xdens -= getCummulativeProb(log(finalSpot/initialSpot),currentVolAtm,idate);
	xdens /= eps*finalSpot;

	double xx = getCummulativeProb(log(finalSpot/initialSpot),currentVolAtm,idate);
	xx = normal_inv(xx);
//	xx -= 0.5*strikeVol;

	double aa = getCummulativeProb(log(U/initialSpot),currentVolAtm,idate);
	aa = normal_inv(aa);

//	aa -= 0.5*strikeVol;
	double a = log(U/initialSpot)/currentVolAtm-0.5*currentVolAtm;

	double newLogSpot	=	getNewLogSpotIncrement(2.0*aa-xx,m_volMap[idate][ipath],idate,ipath);
	double mirrorSpot   =   initialSpot*exp(newLogSpot);

	double ydens = getCummulativeProb(log(mirrorSpot*(1+eps)/initialSpot),currentVolAtm,idate);
	ydens -= getCummulativeProb(log(mirrorSpot/initialSpot),currentVolAtm,idate);
	ydens /= eps*mirrorSpot;


	if ( finalSpot < U )
	{
		prob = 1.0;
	}
	else
	{
		prob = ydens/xdens;

		double test = exp(2.0*(x-a)*a);
	}

	return prob;
*/


//	double prob = bridge(m_dt[idate], U/initialSpot, L, 1,finalSpot/initialSpot, strikeVol);

//	return prob;

	double currentVolAtm  = GetSimulatedVols(ipath,idate);
	double strikeVol		=	GetBridgeVolatility(ipath,idate,0.9*finalSpot);//1.255);//1.3254);//1.235),//*1.33);


	double x = log(finalSpot/initialSpot)/strikeVol;
	int i=1;
	double a = log(U/initialSpot)/strikeVol;

//	double xx = getCummulativeProb(log(finalSpot/initialSpot),currentVolAtm,idate);
//	xx = normal_inv(xx);
//	xx -= 0.5*strikeVol;

//	double aa = getCummulativeProb(log(U/initialSpot),currentVolAtm,idate);
//	aa = normal_inv(aa);
//	aa -= 0.5*strikeVol;


//	if ( finalSpot > U )
	if ( finalSpot < U )
	{
		prob = 1.0;
	}
	else
	{
//		prob = exp(2.0*(xx-aa)*aa);//exp(2.0*(x-a)*a);
		prob = exp(2.0*(x-a)*a);
	}


	return prob;
}	


	
/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: 
**	Returns : nothing
**	Comment :  calculates forward vols between setting iStart and iEnd
****************************************************************/
	
int CForwardSkewMC::GetNumberOfVolStates(int ipath,int idate,int iasset)
{
	return  m_numberVolStates;
}

	
/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: 
**	Returns : nothing
**	Comment :  calculates forward vols between setting iStart and iEnd
****************************************************************/

void CForwardSkewMC::CalculateForwardVols(CVector& impliedVols,const CVector& Strikes,int iStart,int iEnd)
{

//	strikes assumes in terms of forward

	int nstrikes = Strikes.getsize();
	impliedVols.resize(nstrikes);
	
	CVector vals(nstrikes);
	int istrike;

	for ( int ipath = 0 ; ipath < m_nPaths; ipath++ )
	{
		for ( istrike = 0.0 ; istrike < nstrikes; istrike++ )
		{
	vals[istrike] += MlEqMaths::Max((*m_pPath_array)[ipath][m_iasset][iEnd]/(*m_pPath_array)[ipath][m_iasset][iStart]-Strikes[istrike],0.0);
		}
	}

	for ( istrike = 0.0 ; istrike < nstrikes; istrike++ )
	{
		vals[istrike] /= m_nPaths;
	}


// back out implied vol

	double forward	= 1.0;
	double maturity = 0.0;

	for ( int idate = iStart; idate < iEnd; idate++ )
	{
		maturity +=	m_dt[idate];
	}

	double accuracy		= 1e-5;
	double lower_x		= 0.0;
	double upper_x		= 1.5;

	for ( istrike = 0.0 ; istrike < nstrikes; istrike++ )
	{


		impliedVols[istrike] = MlEqBSImpliedVol(
								vals[istrike],		
								forward,			
								maturity,	
								Strikes[istrike],
								1.0,
								1,					// 1 : call, -1 : put, 0 : forward
								eROOTFIND_BRENT_GROWRANGE,		// rootfinder flag (optional) <DEFAULT> 8 </DEFAULT>
								accuracy,		// accuracy (optional) <DEFAULT> 1e-5 </DEFAULT>
								lower_x,			// lower bound to search (optional) <DEFAULT> 0.0 </DEFAULT>
								upper_x			// upper bound to search (optional) <DEFAULT> 1.5 </DEFAULT>
								);
	}

}


/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: destructor
**	Returns : nothing
**	Comment :  
****************************************************************/

 CForwardSkewMC::~CForwardSkewMC()
{
	m_path_array.resize(0,0,0);
	m_volMap.resize(0.0,0.0);
	m_discreteVolMatrix.resize(0.0,0.0);
	m_continuousVolMatrix.resize(0,0);

	m_averageBeta.resize(0);
	m_termvolvol.resize(0);
	m_spotvolvol.resize(0);
	m_prevVol.resize(0);
	m_saveDiscreteVols.resize(0,0);
}

/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: getNewVolSkew
**	Returns : double
**	Comment : calculates the vol skew for newVolAtm 
****************************************************************/

double CForwardSkewMC::getNewVolSkew(EDCubicSpline& baseVolcurve,double xStrike,double newVolAtm,double baseVolAtm,double excessVolVol)
{

	double vol;
	
	if (m_excessVolVolInterp.getsize()== 0 ){
		vol = baseVolcurve.getValue(xStrike)*newVolAtm/baseVolAtm;
	} else {
		vol = baseVolcurve.getValue(xStrike)*pow(newVolAtm/baseVolAtm,1.0+excessVolVol);
	}

	if ( vol <  0.0 && xStrike < -4.0 ){
		vol = 0.005;
	}

	return vol;
}


/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: getExcessVolOfVol
**	Returns: double
**	Action : 
****************************************************************/

double CForwardSkewMC::getExcessVolOfVol(int idate ,double normStrike)
{
	if (m_excessVolVolInterp.getsize()== 0 ){
		return 0.0;
	}

	double val = m_excessVolVolInterp[idate]->getValue(normStrike);
	return val;
}


/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: cummulativeProb
**	Returns: double
**	Action : returns cummulative probability for input strike logK
****************************************************************/
	
	
double CForwardSkewMC::getCummulativeProb(double logK,double currentAtmVol,int idate)
{	
	
	double y = currentAtmVol*m_sqrtdt[idate];
	double x = logK/y;
	double bareVol = m_bareVolAtm;	
	double excessVolOfVol  = getExcessVolOfVol(idate,x);
	double vol = getNewVolSkew(m_bareForwardVolcurve[idate],x,currentAtmVol,bareVol,excessVolOfVol);	
	double eps = 0.001;
	x = (logK+log(1.0+eps))/y;	
	excessVolOfVol  = getExcessVolOfVol(idate,x);

	double deriv = getNewVolSkew(m_bareForwardVolcurve[idate], x, currentAtmVol, bareVol,excessVolOfVol);	
	deriv = (deriv-vol)/(eps*exp(logK));
	
	double d1,d2,cumn;
	double sq = 1.0/sqrt(2.0*3.141592654);
	
    vol *= m_sqrtdt[idate];
	
	d1 = (-logK+0.5*vol*vol)/vol;
	d2 = d1 - vol;
	
	cumn = normal(-d2)+m_sqrtdt[idate]*sq*exp(-0.5*d1*d1)*deriv;	
	return cumn;

}	


////////////////

/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: initHermiteTable
**	Returns: double
**	Action : calibrates forward vols to the cliquet market
****************************************************************/

void CForwardSkewMC::scaleHermiteSkew(CVector& hermiteCoeff, double currentVol,double elasticity)
{
	double oldVol = hermiteCoeff[0];
	hermiteCoeff[0] = currentVol;

	for ( int i= 1 ; i < hermiteCoeff.getsize(); i++ ){
		hermiteCoeff[i] *= pow(currentVol/oldVol,elasticity);
	}

}


/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: initHermiteTable
**	Returns: double
**	Action : calibrates forward vols to the cliquet market
****************************************************************/


void  CForwardSkewMC::initHermiteTable(int numberGridpoints,int idate)
{	
		
	if ( idate == 0 ){
		return;
		// first slice handled as before
	}

	double currentAtmVol;
	for (int j = 0 ; j < m_numberVolStates; j++ )
	{
		currentAtmVol = m_discreteVol[j];
		CVector hermiteCoeff;
		hermiteCoeff = m_hermiteSkew[idate].m_hermitCoeff;
		scaleHermiteSkew(hermiteCoeff,currentAtmVol,m_hermiteElasticity[idate]);

		m_currentHermiteSkews[j].initialize(m_dt[idate],m_fwds[idate][0],hermiteCoeff);
	}

}	






/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: calibrate
**	Returns: double
**	Action : calibrates forward vols to the cliquet market
****************************************************************/


void  CForwardSkewMC::initTable(int numberGridpoints,int idate)
{	
		
	if ( m_hermiteSkew.getsize() > 0 && idate > 0 ){
		initHermiteTable(numberGridpoints,idate);
		return;
	}

//  start setting up guess grids:

	m_lowestNd2 = -5;
	m_highestNd2 = 3;

	m_spotGrid.resize(numberGridpoints,m_numberVolStates);
	m_omegaGridPoints.resize(m_numberVolStates,numberGridpoints);
	
	double currentAtmVol,Nd2;
	int i;

	CVector logK(numberGridpoints);


//  try to find optimal table grid size
		
	for (int j = 0 ; j < m_numberVolStates; j++ )
	{
		m_lowestNd2 = -5;
		m_highestNd2 = 3;

		currentAtmVol = m_discreteVol[j];

		double initialLow = m_lowestNd2;
		double yy,xx;
		int ndim = 500;

		Nd2 = m_lowestNd2;

		for ( i = 0 ; i < ndim ; i++ )	
		{	
			double LogK = (double)Nd2*m_sqrtdt[idate]*currentAtmVol;
			xx= getCummulativeProb(LogK,currentAtmVol,idate);
			yy = normal_inv(xx); 
			
			if ( fabs(yy-m_lowestNd2) < 0.25 || yy < initialLow )
			{
				break;
			}

			Nd2--;
		}

		if ( i == ndim )
		{
			throw( "error in setting up tablegrid" );
		}

		double initialHigh = m_highestNd2;
		Nd2 = m_highestNd2;

		
// check whether Nd2 is overadjusted



		for ( i = 0 ; i < ndim ; i++ )	
		{	
			double LogK = (double)Nd2*m_sqrtdt[idate]*currentAtmVol;
			xx= getCummulativeProb(LogK,currentAtmVol,idate);

			if ( i== 0 && fabs(xx-1.0) < 1e-5)
			{
				for (;;)
				{
					Nd2--;
					LogK	= (double)Nd2*m_sqrtdt[idate]*currentAtmVol;
					xx		= getCummulativeProb(LogK,currentAtmVol,idate);
					yy = normal_inv(xx); 

					if ( yy < m_highestNd2+0.5 )
						break;
				}	
			}

			yy = normal_inv(xx); 
			
			if ( fabs(yy-m_highestNd2) < 0.25 || yy > m_highestNd2 )//sosNd2 )
			{
				break;
			}
			Nd2++;
		}

		if ( i == ndim )
		{
			throw( "error in setting up tablegrid" );
		}

//		approximate optimal grid size found; start setting up spot grid

		for ( i = 0 ; i < numberGridpoints; i++ )
		{
				Nd2 = m_lowestNd2+(double)i/((double)numberGridpoints-1)*(m_highestNd2-m_lowestNd2);

				logK[i] = (double)Nd2*m_sqrtdt[idate]*currentAtmVol;

				double xval=  getCummulativeProb(logK[i],currentAtmVol,idate);
				double val = normal_inv(xval); 
			 
				m_omegaGridPoints[j][i] = val;
				m_spotGrid[i][j] = logK[i];
		}

		CVector temp(m_spotGrid.rows());
		for ( i = 0 ; i < temp.getsize(); i++ )
		{
			temp[i] = m_spotGrid[i][j];
		}

	}			
}	





/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: getExcessVolOfVol
**	Returns: excessVolofVol vector
**	Action : 
****************************************************************/
/*
const CVector CForwardSkewMC::getExcessVolOfVolArray(int idate)
{

	int nsize = m_excessVolofVol.rows();
	static CVector tmp;

	if ( nsize == 0 )
	{
		return tmp;
	}
	else if ( nsize == 1 )
	{
		return m_excessVolofVol[0];
	}
	else
	{
		return m_excessVolofVol[idate];	
	}
}

*/
/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: getLocalVolFactor
**	Returns: nothing
**	Action : calculates spot dependence of forward vols
****************************************************************/

		
void CForwardSkewMC::getLocalVolFactor(double& localVolFactor,double& beta,int idate,int ipath,double volold)		
{

	double maxLocalVolFactor = 3.0;

	if ( idate == 0 )
	{	
		localVolFactor = 1.0;
	}	
	else		 
	{	
		double z;
			
		if ( m_localVolFlag == 1)
		{
			beta = m_beta[idate];
			z = (*m_pPath_array)[ipath][m_iasset][idate]/m_initialSpot;

			localVolFactor = MlEqMaths::Min(pow(z,(double)beta),maxLocalVolFactor);

		}
		else if ( m_localVolFlag == 2)
		{				 
			beta =  m_beta[idate]/volold;
			z = (*m_pPath_array)[ipath][m_iasset][idate]/(*m_pPath_array)[ipath][m_iasset][idate-1];
			localVolFactor = MlEqMaths::Min( pow(z,(double)beta),maxLocalVolFactor);
				
		}	
		else if ( m_localVolFlag == 0)
		{					 
			beta =  m_beta[idate];
			z = (*m_pPath_array)[ipath][m_iasset][idate]/(*m_pPath_array)[ipath][m_iasset][idate-1];
			localVolFactor = MlEqMaths::Min(pow(z,(double)beta),maxLocalVolFactor);
				
		}							
	}	
}	



/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: generatePath
**	Returns: nothing
**	Action : executes one Monte Carlo timestep
****************************************************************/

void  CForwardSkewMC::calculateOptions(CVector& results,MlEqAssetHandle asset, int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor)
{	
	if ( m_assets.size() == 0 && m_nAssets == 1 ){
		// asset have not been set ( e.g. forward skews must have been entered explicitely
		int iasset = 0;
		calculateOptions(results,startIndex,endIndex,Strikes,controlVariate,impliedVolFlag,factor,iasset);
		return;
	}

	int iasset = findAsset(asset);
	if ( iasset != 0 ){
		throw("only one asset should be set in single factor monte carlo");
	}

	calculateOptions(results,startIndex,endIndex,Strikes,controlVariate,impliedVolFlag,factor,iasset);
}


/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: generatePath
**	Returns: nothing
**	Action : executes one Monte Carlo timestep
****************************************************************/

void CForwardSkewMC::calculateOptions(CVector& results,int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor,int iasset)
{	
	
	if ( iasset > 0 ){
		throw("iasset for a single factor monte carlo need to be set to zero");
	}

	int ipath,istrike;
	results.resize(Strikes.size());
	CVector values(Strikes.size());
	CVector ControlVarValues(Strikes.size());
	
	vector< MlEqStrike >  strikes; 
	strikes.resize(Strikes.size());
	
	for ( istrike = 0 ; istrike < strikes.size(); istrike++ )
	{
		MlEqStrike::convertStrikes(strikes[istrike],*Strikes[istrike]);
	}


	double payoff;
	for ( ipath = 0 ; ipath < m_nPaths; ipath++ )
	{
		for ( istrike = 0 ; istrike < strikes.size(); istrike++ )
		{
			payoff = MlEqMaths::Max(factor*GetPathValue(ipath,endIndex)/GetPathValue(ipath,startIndex)-strikes[istrike].m_strike,0.0);
			completePayoff(payoff,ipath);
			values[istrike] += payoff;
		}
	}

	averageMonteCarloResults(values);
	for ( int i = 0; i < values.getsize(); i++ ){
		values[i] /= factor;
	}

//  add controlVariates	
	if ( controlVariate )
	{

		if ( m_controlVariate == 0 )
		{
			throw("trying to use control variate which is not set");
		}


		for ( ipath = 0 ; ipath < m_nPaths; ipath++ )
		{
			for ( istrike = 0 ; istrike < strikes.size(); istrike++ )
			{
		
				ControlVarValues[istrike] += MlEqMaths::Max(m_controlVariateCurrentSpot[ipath]-strikes[istrike].m_strike,0.0);
			}
		}

		averageMonteCarloResults(ControlVarValues);



	for ( istrike = 0 ; istrike < results.getsize(); istrike++ )
	{
		results[istrike] = GetImpliedVol(startIndex,endIndex,ControlVarValues[istrike],strikes[istrike].m_strike);
	}


		for ( istrike = 0 ; istrike < strikes.size(); istrike++ )
		{
			ControlVarValues[istrike] -= Bs(m_termFwds[endIndex][0],m_controlVariateVols[endIndex-1],m_t[endIndex],strikes[istrike].m_strike,1.0,1.0);
		}

		for ( istrike = 0 ; istrike < strikes.size(); istrike++ )
		{
			 values[istrike] = values[istrike] -ControlVarValues[istrike];
		}


	}	
	else if ( m_controlVariatePrices.rows() )
	{
		if ( m_controlVariatePrices.cols() != strikes.size() ){
			throw("incorrect number of columns in m_controlVariatePrices");}

		CVector temp(ControlVarValues.getsize());
		for ( istrike = 0 ; istrike < strikes.size(); istrike++ )
		{
			ControlVarValues[istrike] = m_controlVariatePrices[endIndex-1][istrike];
		}


		for ( istrike = 0 ; istrike < strikes.size(); istrike++ )
		{
			ControlVarValues[istrike] -= Bs(m_termFwds[endIndex][0],m_controlVariateVols[endIndex-1],m_t[endIndex],strikes[istrike].m_strike,1.0,1.0);
		}

		for ( istrike = 0 ; istrike < strikes.size(); istrike++ )
		{
			 values[istrike] = values[istrike] -ControlVarValues[istrike];
		}


	}	

//	call implied volcalculator here

	if ( impliedVolFlag )
	{
		for ( istrike = 0 ; istrike < results.getsize(); istrike++ )
		{
			results[istrike] = GetImpliedVol(startIndex,endIndex,values[istrike],strikes[istrike].m_strike/factor);
		}
	}
	else
	{
		for ( istrike = 0 ; istrike < results.getsize(); istrike++ )
		{
			results[istrike] = values[istrike];
		}
	}
}	
	

/***************************************************************
**	Class  : MlEqMonteCarlo
**	Routine: generatePath
**	Returns: nothing
**	Action : executes one Monte Carlo timestep
****************************************************************/

void MlEqMonteCarlo::calculateOptions(CVector& results,MlEqAssetHandle asset,int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor)
{
	if ( m_nAssets == 0 && m_assets.size() == 0 )
	{
		// forward skew mc was initialized without an asset
		int iasset = 0;
		calculateOptions(results,startIndex,endIndex,Strikes,controlVariate,impliedVolFlag,factor,iasset);
		return;

	}
	
	int	iasset = findAsset(asset);
	calculateOptions(results,startIndex,endIndex,Strikes,controlVariate,impliedVolFlag,factor,iasset);
}

/***************************************************************
**	Class  : MlEqMonteCarlo
**	Routine: generatePath
**	Returns: nothing
**	Action : executes one Monte Carlo timestep
****************************************************************/

void MlEqMonteCarlo::calculateOptions(CVector& results,int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor,int iasset)
{	
	
	int ipath,istrike;
	results.resize(Strikes.size());
	CVector values(Strikes.size());
	CVector ControlVarValues(Strikes.size());
	
	vector< MlEqStrike >  strikes; 
	strikes.resize(Strikes.size());
	
	for ( istrike = 0 ; istrike < strikes.size(); istrike++ )
	{
		MlEqStrike::convertStrikes(strikes[istrike],*Strikes[istrike]);
	}


	double payoff;
	for ( ipath = 0 ; ipath < m_nPaths; ipath++ )
	{
		for ( istrike = 0 ; istrike < strikes.size(); istrike++ )
		{
			payoff = MlEqMaths::Max(factor*GetPathValue(ipath,endIndex,iasset)/GetPathValue(ipath,startIndex,iasset)-strikes[istrike].m_strike,0.0);
			completePayoff(payoff,ipath);
			values[istrike] += payoff;
		}
	}

	averageMonteCarloResults(values);
	for ( int i = 0; i < values.getsize(); i++ ){
		values[i] /= factor;
	}


//	call implied volcalculator here

	if ( impliedVolFlag )
	{
		for ( istrike = 0 ; istrike < results.getsize(); istrike++ )
		{
			results[istrike] = GetImpliedVol(startIndex,endIndex,values[istrike],strikes[istrike].m_strike/factor,iasset);
		}
	}
	else
	{
		for ( istrike = 0 ; istrike < results.getsize(); istrike++ )
		{
			results[istrike] = values[istrike];
		}
	}
}	
	


/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: getSlope
**	Returns: nothing
**	Action : executes one Monte Carlo timestep
****************************************************************/

int CForwardSkewMC::getSlope(double& slope,double & smile,double& maxDiff,double& minDiff,CVector& impliedVols,CVector& termVols,int atmIndex,CVector& accuracy)	
{
		int j;
		double x;

		slope = 0.0;
		smile = 0.0;
		int Suceed= 1;
//			calibrate to slope first
			maxDiff=-1e99;minDiff=1e99; 
			for ( j = 0 ; j < impliedVols.getsize(); j++ )
			{
				if ( j == atmIndex )continue;

				x =  (impliedVols[j]-impliedVols[atmIndex])
					-(termVols[j]-termVols[atmIndex]);

				if ( j < atmIndex )
				{
					slope += x;
				}
				else
				{
					slope -= x;
				}

				smile += x;

				maxDiff = MlEqMaths::Max(maxDiff ,x);
				minDiff = MlEqMaths::Min(minDiff ,x);
			}

			if ( MlEqMaths::Max(fabs(maxDiff),fabs(minDiff)) > accuracy[j] )
			{		 
				Suceed = 0;
			}
		

		return Suceed;		
}

/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: updateSlope
**	Returns: nothing
**	Action : executes one Monte Carlo timestep
****************************************************************/

int CForwardSkewMC::updateSlope(CVector& impliedVols,CVector& termVols,int atmIndex,CVector& accuracy)	
{

	int suceed;
	double slope,smile,maxDiff,minDiff;
	if ( suceed = getSlope(slope,smile,maxDiff,minDiff,impliedVols,termVols,atmIndex,accuracy) )
	{
		return suceed;
	}
	
	double compare = 0.0;
	for ( int i = 0 ; i < accuracy.getsize(); i++ )
	{
		compare += accuracy[i];
	}

	double slopeFac = slope/termVols[atmIndex];
	double smileFac = smile/termVols[atmIndex];
	double betaAdjust	= 1.0+0.5*slopeFac;
	double smileAdjust	= 1.0+0.5*smileFac;

	int idate;
	for ( idate = 0; idate < termVols.getsize(); idate++ )
	{		
		m_beta[idate]	*=	slopeFac ;		
		m_volvol[idate]	*=	smileAdjust;	
	}
	return 0;
}

/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: GlobalCalibrate
**	Returns: nothing
**	Action : executes one Monte Carlo timestep
****************************************************************/


void CForwardSkewMC::GlobalCalibrate(CMatrix& TermSkews,vector<vector<MlEqStrikeHandle> > & Strikes,
					 long& idum,CVector& accuracy)// calibrates paths from idate->idate+1
{	
	// start calibrating beta on the short end first	
	int i,idateCalib;

	CVector volAtmSave;
	volAtmSave = m_VolAtm;

	m_volvolMonitor.resize(2);	
	int calibGlobalFlag = 0;
	m_calibResults.resize(m_calibTimeIndex.getsize(),6);

//  find last slice	that is calibrated to beta

	int j=-1;	
	for ( i = 0 ; i < m_calibrationInfo.rows(); i++ )
	{
		if ( m_calibrationInfo[i][1] == 1 ){
			j = i;
		}
	}

	if ( j == -1 ){
		throw("no beta calibration points have been entered");
	}

	int idateCalibIndex = 0;
	idateCalib = m_calibTimeIndex[idateCalibIndex];
	calibGlobalFlag = 1;

	m_calibResults.resize(m_calibTimeIndex.getsize(),6);
	calibrateTimeSlice(idateCalibIndex,m_calibVols[idateCalibIndex],m_pCalibStrikes[idateCalib-1],idum,accuracy,calibGlobalFlag,m_calibrationInfo);

	for (  i = 1; i <= j; i++ )
	{
		calibGlobalFlag = 0;
		calibrateTimeSlice(i,m_calibVols[i],m_pCalibStrikes[m_calibTimeIndex[i]-1],idum,accuracy,calibGlobalFlag,m_calibrationInfo);
	}


	double avBeta=0;
	int k = 0;
	for ( i = 0 ; i <= j; i++ )
	{	
		if ( m_calibrationInfo[i][1] == 1 )
		{
			k++;
			idateCalib = m_calibTimeIndex[i];
			double xbeta = m_beta[idateCalib]/m_VolAtm[idateCalib-1]; //sos
			avBeta += xbeta;
		}
	}
	
	if ( k == 0 ){
		throw("no beta calibration point has been added in global calibration");
	}

	avBeta /= (double)k;

	// reset beta to average beta and restart
	for ( i = 0 ; i < m_beta.getsize(); i++ ){
		m_beta[i] = avBeta*m_VolAtm[i];
	}
	


//  find first calibration start date index for mreversion

	int lastbetaIndex = j;
	int firstMeanRevIndex = m_calibrationInfo.rows()-1;
	for ( i = 0 ; i < m_calibrationInfo.rows(); i++ )
	{
		if ( m_calibrationInfo[i][3] == 1 )
		{
			firstMeanRevIndex = i;
			break;
		}
	}

	if ( lastbetaIndex >= firstMeanRevIndex ){
		throw("inconsistent calibration set up");
	}

	if ( firstMeanRevIndex == m_calibrationInfo.rows() ){
		throw("no meanreversion calib entered");
	}

//  restart from beginning with average beta

	CMatrix calibInfo(m_calibrationInfo.rows(),4);

	for ( i = 0 ; i < m_calibrationInfo.rows(); i++ )
	{
		calibInfo[i][0] = 1.0;
		calibInfo[i][1] = 0.0;
		calibInfo[i][2] = 0.0;
		calibInfo[i][3] = 0.0;
	}
	
	for ( i = 0 ; i < firstMeanRevIndex; i++ )
	{	
		calibGlobalFlag = 0;
		calibrateTimeSlice(i,m_calibVols[i],m_pCalibStrikes[m_calibTimeIndex[i]-1],idum,accuracy,calibGlobalFlag,calibInfo);
	}

//  start calibrating to meanreversion 

// find last kappa calibration point:

	int lastMeanRevIndex=firstMeanRevIndex;
	for ( i = firstMeanRevIndex+1 ; i < m_calibrationInfo.rows(); i++ )
	{
		if ( m_calibrationInfo[i][3] == 1 ){
			lastMeanRevIndex = i;
		}
	}

	if ( lastMeanRevIndex != firstMeanRevIndex ){
		throw("only calibration to one meanreversionrate implemented yet");
	}

//	for ( i = 0 ; i < m_calibrationInfo.rows(); i++ )
//	{
//		calibInfo[i][0] = 1.0;
//		calibInfo[i][1] = 0.0;
//		calibInfo[i][2] = 0.0;
//		calibInfo[i][3] = 1.0;
//	}


	calibGlobalFlag = 1;
	int suceedType  = 3;
	int maxiterlocal=12;
	
	CVector Accuracy(m_calibVols.cols());

	for ( i = 0 ; i < Accuracy.getsize(); i++ ){
		Accuracy[i] = 0.01;
	}

	idateCalib = m_calibTimeIndex[lastMeanRevIndex]-1;//sos
	if ( idateCalib < 0 ){
		throw("calibration to meanreversion must occur at a later date");
	}

	calibrateTimeSlice(lastMeanRevIndex,m_calibVols[lastMeanRevIndex],m_pCalibStrikes[idateCalib],idum,Accuracy,calibGlobalFlag,m_calibrationInfo,maxiterlocal,suceedType);	


//  continue to end


	suceedType = 1;
	calibGlobalFlag = 0;
	for ( j = lastMeanRevIndex+1  ; j < m_calibTimeIndex.getsize(); j++ )
	{
		idateCalib = m_calibTimeIndex[j];
		calibrateTimeSlice(j,m_calibVols[j],m_pCalibStrikes[idateCalib-1],idum,accuracy,calibGlobalFlag,m_calibrationInfo);
	}


	for ( i = 0 ; i < m_calibTimeIndex.getsize(); i++ )
	{
		idateCalib = m_calibTimeIndex[i];

		m_calibResults[i][0] = idateCalib;
		m_calibResults[i][1] = m_VolAtm[idateCalib-1]/volAtmSave[idateCalib-1];
		m_calibResults[i][2] = m_beta[idateCalib-1]/m_VolAtm[idateCalib-1];
		m_calibResults[i][3] = m_volvol[idateCalib-1];
		m_calibResults[i][4] = m_meanReversionRate[idateCalib-1];
		m_calibResults[i][5] = 0;;

	}

}			
	

	
/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: calibrateTimeSlice
**	Returns: nothing
**	Action : executes one Monte Carlo timestep
****************************************************************/
	
	
int CForwardSkewMC::calibrateTimeSlice(int idateCalibIndex,CVector& TermSkews,	vector< MlEqStrikeHandle > & Strikes,
					 long& idum,CVector& accuracy,int calibGlobalFlag,CMatrix& calibrationInfo,int maxiter,int succeedType,int firstIdateStart)// calibrates paths from idate->idate+1
{	

// function calibrates paths from idate->idate+1
// if calibrationInfo columns are all set to zero routine will just evolve
// from timestep idate = firstidateStart; idate < idateCalib; idate++ )
// if maxiter == 0 routine will also check objective function for initial gues but not perform any iterations
// function returns value of objective function; if negative calibration has failed

	
	int idateCalib	= m_calibTimeIndex[idateCalibIndex];
	
// routine calibrates patharray[idateCalib]

	int fillFlag	= 1;
	int useRandom   = 0;
	int idate;

	if ( idateCalibIndex == 0 ){
		initializeRandomGenerator(m_numberOfFactors);
	}


	int idateStart=0;
	int firstidateStart = 0;
	if ( firstIdateStart != -1 ){
		firstidateStart = firstIdateStart;
	}

	if ( !calibGlobalFlag && idateCalibIndex>0 )
	{
		idateStart = m_calibTimeIndex[idateCalibIndex-1];
		firstidateStart = idateStart;
	}
	
	if ( idateStart == 0 ){
		m_volvolMonitor.resize(2);
	}

    if ( m_newcurrentVols.getsize() == 0 )
	{
		m_newcurrentVols.resize(m_nPaths);//[ipath]
		m_currentVols.resize(m_nPaths);//[ipath]
	}
	
	for ( idate = firstidateStart; idate < idateCalib; idate++ )
	{	
		
		GenerateRandomNumbers(idate);
		int calibflag = 2;
		if ( idate == 0 )
		{
			m_pcurrentVols		= &m_currentVols;
			m_pnewcurrentVols	= &m_newcurrentVols;
			m_psavecurrentVols	= &m_currentVols;
		}
		else 
		{
//			no calibration done here

			if ( &(*m_psavecurrentVols)	== &m_newcurrentVols )
			{
				if ( idate > idateStart ){
					m_pcurrentVols		=& m_currentVols;
				}
				m_pnewcurrentVols	= &m_currentVols;
			}
			else if ( &(*m_psavecurrentVols)	== &m_currentVols )
			{
				if ( idate > idateStart ){
					m_pcurrentVols	= &m_newcurrentVols;
				}
				m_pnewcurrentVols	= &m_newcurrentVols;
			}
			else{
				throw("incorrect pointzer arithmatic encountered in calibration routine");
			}
		}

		generatePath(idate,calibflag,idum);
	}

	// check if any calibration needs to be performed at all

	int continueCalibration=0;
	for (int i = 0 ; i < calibrationInfo.cols(); i++ )
	{
		if ( calibrationInfo[idateCalibIndex][i] )
		{
			continueCalibration = 1;
			break;
		}
	}
	if ( continueCalibration == 0 ){// no vcalibration required in this case
		return 0;
	}
	
	CVector impliedVols;
	int impliedvolflag = 1;
	calculateOptions(impliedVols,0,idateCalib,Strikes,m_controlVariate,impliedvolflag,m_initialSpot);


	int suceed = 1;
//	int maxiter = 12;
	double x,volatm;

	int atmIndex;
	if ( !MlEqVolatilityStructure::findATMFwdVol(atmIndex,Strikes,m_termFwds[idateCalib][0]) ){
		throw("please add atm forward vol as calibration point");
	}
	volatm = TermSkews[atmIndex];
	
	fitPolynomial fittingCurve,mktCurve;
	
	CMatrix NonZeroPartialDerivativesSpecification;
	double FinalTolerance		= 1e-12;
	double StoppingTolerance	= 1e-12;
	double InitialTolerance		= 2*FinalTolerance;
	int outputFlag				= 0;
	
	CVector initialGuess(3);
	initialGuess[0]	= TermSkews[atmIndex];
	initialGuess[1]	= -2.0;
	initialGuess[2]	= 0.0;
	
	CMatrix xValBounds(3,2);
	xValBounds[0][0] = 0.5*initialGuess[0];
	xValBounds[0][1] = 1.5*initialGuess[0];
	xValBounds[1][0] = -5.0;
	xValBounds[1][1] = 5.0;;
	xValBounds[2][0] = -5.0;
	xValBounds[2][1] = 0.1;
	
	vector < MlEqNormalizedStrike > normStrikes;
	normStrikes.resize(Strikes.size());
	
	double mat	= m_t[idateCalib];
	double fwd	= 1.0;
	for ( int i = 0 ; i < idateCalib-1; i++ ){
		fwd *= m_fwds[i][0];
	}
	double strike	= 1.0;

	const std::vector<long>& mcDates	= GetMCDates();
	
	for (int istrike = 0 ; istrike < Strikes.size(); istrike++ )
	{
		normStrikes[istrike].initialize(strike,mcDates[idateCalib],volatm,m_termFwds[idateCalib][0],false,m_hDate);//fwdsos 
		MlEqStrike::convertStrikes(normStrikes[istrike],*Strikes[istrike]);
	}

	
	CMatrix mktdata(TermSkews.getsize(),2);
	for ( int istrike = 0 ; istrike < Strikes.size(); istrike++ )
	{
		mktdata[istrike][0] = normStrikes[istrike].m_strike;
		mktdata[istrike][1] = TermSkews[istrike];
	}
	
	
	CMatrix ObjectiveBounds;
	mktCurve.initialize(initialGuess,mktdata,
					    xValBounds,ObjectiveBounds,
					    InitialTolerance,FinalTolerance,StoppingTolerance,
					    NonZeroPartialDerivativesSpecification,outputFlag);
	
	int maximizeFlag = 0;
	int returnCode;
	
	CVector mktParams(initialGuess);
	mktCurve.solve( mktParams, maximizeFlag, returnCode);
	
//	fittingCurve: initialize curve
	
	CMatrix fitdata(TermSkews.getsize(),2);
	
	for (int istrike = 0 ; istrike < Strikes.size(); istrike++ )
	{
		fitdata[istrike][0] = normStrikes[istrike].m_strike;
		fitdata[istrike][1] = impliedVols[istrike];
	}
	
	fittingCurve.initialize(initialGuess,fitdata,
					    xValBounds,ObjectiveBounds,
					    InitialTolerance,FinalTolerance,StoppingTolerance,
					    NonZeroPartialDerivativesSpecification,outputFlag);
	
	CVector fittingParams(initialGuess);
	fittingCurve.solve( fittingParams, maximizeFlag, returnCode);

	// check initial guess

	double fitxslope;
	suceed = 1;

	if ( succeedType == 1 )
	{
		for (int i = 0 ; i < accuracy.getsize(); i++ )
		{
			if ( calibrationInfo[idateCalibIndex][i] )
			{
				if ( ( x = fabs(fittingParams[i]-mktParams[i]) ) > accuracy[i] )
				{	
					suceed = 0;
					break;
				}	
			}		

		}
	}
	else if ( succeedType == 2 )
	{
		for (int i = 0 ; i < TermSkews.getsize(); i++ )
		{
			if ( ( x = fabs(impliedVols[i]-TermSkews[i]) ) > accuracy[i])
			{
				suceed = 0;
				break;
			}
		}
	}
	else
	{
		fitxslope=0.0;
		for (int i = 0 ; i < TermSkews.getsize(); i++ )
		{
			if ( i < atmIndex ){
				fitxslope += impliedVols[i]-TermSkews[i];
			}
			else{
				fitxslope -= impliedVols[i]-TermSkews[i];
			}

		}

		if ( (x=fabs(fitxslope)) > accuracy[0] ){
			suceed = 0;
		}		
	}


	if ( suceed == 1 )
	{
		m_psavecurrentVols  = m_pnewcurrentVols;
		m_pcurrentVols		= m_pnewcurrentVols;

		m_calibResults[idateCalibIndex][0] = idateCalib-1;
		m_calibResults[idateCalibIndex][1] = 1.0;
		m_calibResults[idateCalibIndex][2] = m_beta[idateCalib-1]/m_VolAtm[idateCalib-1];
		m_calibResults[idateCalibIndex][3] = m_volvol[idateCalib-1];
		m_calibResults[idateCalibIndex][4] = m_meanReversionRate[idateCalib-1];
		m_calibResults[idateCalibIndex][5] = 0;

		m_volvolMonitor[0] += mktParams[2];
		m_volvolMonitor[1] += fittingParams[2];

		return x;
	}

	if ( maxiter ==0  ){
		return 0;
	}

	double beta0	=  0.0;
	double beta1	=  m_beta[idateCalib-1]/m_VolAtm[idateCalib-1];
	double slope0	=  0.0;
	double slope1	=  fittingParams[1];	
	double beta2,slope,factor,mrev2,mrslope;

	if ( m_rev0 < 0 && calibrationInfo[idateCalibIndex][3] )
	{
		m_rev0		= 1.06667*m_meanReversionRate[idateCalib-1];
		m_rev1		= m_meanReversionRate[idateCalib-1];
		m_xslope0	= 1.2*fitxslope;
	}
	if ( calibrationInfo[idateCalibIndex][3] ){
		m_xslope1	= fitxslope;
	}
	

	double modelvolatm = impliedVols[atmIndex];
	double f = 1.0;
	double volvolslope0=0.0;
	double volvolslope1=fittingParams[2];
	double volvol2,volvol0 = 0.0;
	double volvol1 = m_volvol[idateCalib-1];
	double /*f0,f1,f2,*/fslope0,fslope1;
	int iter;

	for ( iter = 0; iter < maxiter; iter++ )
	{	

//      update level

		if ( calibrationInfo[idateCalibIndex][0] )
		{
			factor = m_dt[idateCalib-1]/m_t[idateCalib]*m_VolAtm[idateCalib-1]/volatm;
			factor = 1.0-(fittingParams[0]-mktParams[0])/m_VolAtm[idateCalib-1]/factor;

			if ( factor < 0.0 ){
				factor = 0.5;
			}

			if ( factor > 4.0 ){
				factor = 4.0;
			}

			f *= factor;
			multiplyForwardVols(idateCalib-1,factor);

		}

//		update skew

		if ( calibrationInfo[idateCalibIndex][1] )
		{
			slope = (slope1-slope0)/(beta1-beta0);
			beta2 = -(slope1-mktParams[1])/slope + beta1;
		
			if ( calibGlobalFlag )
			{

				if ( fabs(beta2)< 1e-3 && iter == 0  ){
					throw( "zero beta encountered during calibration. Have U set initial guess to zero?");
				}

				for (int i = 0 ; i < m_beta.getsize(); i++ ){
					m_beta[i] = beta2*m_VolAtm[i];
				}
			}
			else{
				m_beta[idateCalib-1] = beta2*m_VolAtm[idateCalib-1];
			}
		}
		else if ( calibrationInfo[idateCalibIndex][0] )
		{//sos
			// update betas anyway

			m_beta[idateCalib-1] *= factor;

		}

			
//      update volvol

		if ( calibrationInfo[idateCalibIndex][2] )
		{

			slope = (volvolslope1-volvolslope0)/(volvol1-volvol0);
			volvol2 = -(volvolslope1-mktParams[2])/slope+volvol1;

			if ( volvol2 < 0.0 )
			{
				volvol2 = volvol1;
//				throw("negative vol of vol encountered during calibration");
			}
			
			if ( calibGlobalFlag )
			{
				for ( int j = 0 ; j < m_volvol.getsize(); j++ ){
					m_volvol[j] = volvol2;
				}
			}
			else{
				m_volvol[idateCalib-1] = volvol2;
			}
		}

		
//		update meanreversion


		if ( calibrationInfo[idateCalibIndex][3] )
		{
		
			mrslope = (m_xslope1-m_xslope0)/(m_rev1-m_rev0);
			mrev2 = -(m_xslope1-0.0)/mrslope + m_rev1;
		
			if ( calibGlobalFlag )
			{
				for (int i = 0 ; i < m_meanReversionRate.getsize(); i++ ){
					m_meanReversionRate[i] = mrev2;
				}
			}
			else{
				m_meanReversionRate[idateCalib-1] = mrev2;
			}

			if ( m_mrUpdateOnly > 0)
			{
				m_rev0			= m_rev1;
				m_rev1			= mrev2;
				m_xslope0		= fitxslope;
				// m_xslope1 was set earlier

				return 1e10;
			}
		}

		fillFlag	= 0;
		useRandom	= 1;

		for ( idate = idateStart; idate < idateCalib; idate++ )
		{

			if ( idate == idateStart )
			{
				m_pcurrentVols = m_psavecurrentVols;
				if ( &(*m_pcurrentVols)	== &m_newcurrentVols ){
					m_pnewcurrentVols	= &m_currentVols;
				}
				else if ( &(*m_pcurrentVols)	== &m_currentVols ){
					m_pnewcurrentVols	= &m_newcurrentVols;
				}
				else{
					throw("incorrect pointzer arithmatic encountered in calibration routine");
				}

			}
			else if ( idate > idateStart )
			{
				if ( &(*m_psavecurrentVols)	== &m_newcurrentVols )
				{
					m_pcurrentVols		= &m_currentVols;
					m_pnewcurrentVols	= &m_currentVols;
				}
				else if ( &(*m_psavecurrentVols) == &m_currentVols )
				{
					m_pcurrentVols		= &m_newcurrentVols;
					m_pnewcurrentVols	= &m_newcurrentVols;
				}
				else{
					throw("incorrect pointzer arithmatic encountered in calibration routine");
				}
			}
			int calibflag = 0;
			m_controlVariateCalcFlag = 0;
			generatePath(idate,calibflag,idum);
		}
		
		calculateOptions(impliedVols,0,idateCalib,Strikes,m_controlVariate,impliedvolflag,m_initialSpot);
		for (int istrike = 0 ; istrike < Strikes.size(); istrike++ ){
			fitdata[istrike][1] = impliedVols[istrike];
		}
		
		fittingCurve.initialize(initialGuess,fitdata,
							xValBounds,ObjectiveBounds,
							InitialTolerance,FinalTolerance,StoppingTolerance,
							NonZeroPartialDerivativesSpecification,outputFlag);
	
		fittingCurve.solve( fittingParams, maximizeFlag, returnCode);

		modelvolatm = impliedVols[atmIndex];

		suceed = 1;
		if ( succeedType == 1 )
		{
			for (int i = 0 ; i < accuracy.getsize(); i++ )
			{
				if ( calibrationInfo[idateCalibIndex][i] )
				{
					if ( ( x = fabs(fittingParams[i]-mktParams[i]) ) > accuracy[i] )
					{	
						suceed = 0;
						break;
					}	
				}	

			}
		}
		else if ( succeedType == 2 )
		{
			for (int i = 0 ; i < TermSkews.getsize(); i++ )
			{
				double x = fabs(impliedVols[i]-TermSkews[i]); 
				if ( x > accuracy[i] )
				{
					suceed = 0;
					break;
				}	
			}
		}
		else
		{
			fitxslope=0.0;
			for (int i = 0 ; i < TermSkews.getsize(); i++ )
			{
				if ( i < atmIndex ){
					fitxslope += impliedVols[i]-TermSkews[i];
				}
				else{
					fitxslope -= impliedVols[i]-TermSkews[i];
				}
			}

			if ( (x = fabs(fitxslope) ) > accuracy[0] ){
				suceed = 0;
			}

		}

		if ( suceed == 1 )
		{

			m_psavecurrentVols  = m_pnewcurrentVols;
			m_pcurrentVols		= m_pnewcurrentVols;

			m_calibResults[idateCalibIndex][0] = idateCalib;
			m_calibResults[idateCalibIndex][1] = f;
			m_calibResults[idateCalibIndex][2] = m_beta[idateCalib-1]/m_VolAtm[idateCalib-1];
			m_calibResults[idateCalibIndex][3] = m_volvol[idateCalib-1];
			m_calibResults[idateCalibIndex][4] = m_meanReversionRate[idateCalib-1];
			m_calibResults[idateCalibIndex][5] = iter;

			m_controlVariateCalcFlag = m_controlVariate;

			m_volvolMonitor[0] += mktParams[2];
			m_volvolMonitor[1] += fittingParams[2];

			return x;
		}

		if ( calibrationInfo[idateCalibIndex][0] )
		{
			//f0				= f1;
			//f1				= f2;
			fslope0			= fslope1;
			fslope1			= fittingParams[0];
		}

		if ( calibrationInfo[idateCalibIndex][1] )
		{
			beta0			= beta1;
			beta1			= beta2;
			slope0			= slope1;
			slope1			= fittingParams[1];
		}

		if ( calibrationInfo[idateCalibIndex][2] )
		{
			volvol0			= volvol1;
			volvol1			= volvol2;
			volvolslope0	= volvolslope1;
			volvolslope1	= fittingParams[2];
		}

		if ( calibrationInfo[idateCalibIndex][3] )
		{
			m_rev0			= m_rev1;
			m_rev1			= mrev2;
			m_xslope0		= m_xslope1;
			m_xslope1		= fitxslope;
		}
			
//		update atm vols
//		update skew
//		update volvol		
		
	}	
		
	if ( iter == maxiter )
	{		
			int calibFail = 1;

//			keep going either way: do more work here

			m_psavecurrentVols  = m_pnewcurrentVols;
			m_pcurrentVols		= m_pnewcurrentVols;

			m_calibResults[idateCalibIndex][0] = idateCalib;
			m_calibResults[idateCalibIndex][1] = f;
			m_calibResults[idateCalibIndex][2] = m_beta[idateCalib-1]/m_VolAtm[idateCalib-1];
			m_calibResults[idateCalibIndex][3] = m_volvol[idateCalib-1];
			m_calibResults[idateCalibIndex][4] = m_meanReversionRate[idateCalib-1];
			m_calibResults[idateCalibIndex][5] = iter;

			m_controlVariateCalcFlag = m_controlVariate;

			m_volvolMonitor[0] += mktParams[2];
			m_volvolMonitor[1] += fittingParams[2];

			return x;
	}

	return x;
}



/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: generatePath
**	Returns: nothing
**	Action : executes one Monte Carlo timestep
****************************************************************/


void CForwardSkewMC::generatePath(int idate,int calibflag,long& idum)// generates paths from idate->idate+1
{		

	if ( m_nDates == 0 ){
		return;
	}

	double rand,volrand,volold,volnew,volfwd;
	double y,z,localVolFactor,beta;		
	int ipath;

	int n_offset = m_numberOfFactors*idate;
		
//	start generating the next volstates first
	
	if ( idate == 0 )
	{
		m_initialShortVol = m_bareForwardVolcurve[0].getValue(0.0);
		(*m_pcurrentVols)[0] = m_initialShortVol;
	}

	double	averageVol = 0.0;		
	double  volatm0 = (*m_pcurrentVols)[0];
		
	if ( idate > 0 ){
		volfwd = m_VolAtm[idate]/m_VolAtm[idate-1];
	}else{
		volfwd = 1.0;
	}
		
	m_currentVariance = 0.0;
	m_currentAverage = 0.0;
		
	volold = volatm0;

	double maxlocalvolfactor = -1e99;	
	double x,avbeta=0.0;
	for ( ipath = 0 ; ipath < m_nPaths; ipath++ )
	{			
		if ( idate > 0 ){	
			volold = (*m_pcurrentVols)[ipath];
		}	
				
		volrand =  GetNextRandom(ipath,n_offset);				
		double volvol = m_volvol[idate];
		getLocalVolFactor(localVolFactor,beta,idate,ipath,volold);		
		maxlocalvolfactor = MlEqMaths::Max(maxlocalvolfactor,localVolFactor);
		
		y	=   exp(-0.5*volvol*volvol*m_dt[idate]+volvol*m_sqrtdt[idate]*volrand);	
		z   =   m_meanReversionRate[idate]+pow((volold/m_meanReversionLevel[idate]-1.0),2.0);
		z	=   exp(-z*m_dt[idate]);
		x	=   volold*y*localVolFactor-(volold*y*localVolFactor-m_meanReversionLevel[idate])*(1.0-z);
		
		volnew  = volfwd*x;

		if ( volnew < 0.05 ){
			volnew = 0.05;
		}

		averageVol+= volnew;
		
		if ( calibflag == 1){	
			m_currentAverage	+= localVolFactor;
			m_currentVariance   += localVolFactor*localVolFactor;
		}	
		
		if ( volold < 1e-6){
			throw( "zero vol generated" );
		}

		(*m_pnewcurrentVols)[ipath] = volnew;	

		if ( m_saveVolofVolInfo ){
			avbeta+= beta;			
		}
	}		
	
	averageVol /= m_nPaths;
	if ( m_saveVolofVolInfo )
	{
		avbeta /= m_nPaths;
		if ( idate == 0 ){
			m_averageBeta[idate] = 1.0;
		}else{
			m_averageBeta[idate] = avbeta;
		}
	}

	if ( calibflag == 1)
	{	
		m_currentVariance	/= m_nPaths;
		m_currentAverage	/= m_nPaths;
		m_currentVariance	-= m_currentAverage*m_currentAverage;
	}	
			
	double volmultiplyer	 = m_VolAtm[idate];
	double volatm			 = m_VolAtm[idate];
		
	for ( ipath = 0; ipath < m_nPaths; ipath++ ){	
		(*m_pnewcurrentVols)[ipath] *= volmultiplyer/averageVol;
	}	
	
	double factor = 1.0;
	double vol;
		
	calibrate(idate,*m_pnewcurrentVols,factor,calibflag);
	
	iVector indexVec(m_nPaths);	
	for ( ipath = 0; ipath < m_nPaths; ipath++ )
	{		
		vol = (*m_pnewcurrentVols)[ipath];
		vol *= m_volmultiplyer;
		
		int index=0;
		MlEqMaths::locate(m_discreteVolArray,vol,index);
		index++;
		if ( index == m_discreteVol.getsize() )
			index--;
		
		indexVec[ipath] = index;
	}	
	
	if ( m_saveVolGrid == 1)
	{	
		m_volMap[idate] = indexVec;
		m_discreteVolMatrix[idate]= m_discreteVol;
	}	
	else if ( m_saveVolGrid == 2){	
		m_continuousVolMatrix[idate] = (*m_pnewcurrentVols);
	}	

	
//  start creating spot states
			
	initTable(m_numberGridpoints,idate);
	
	double newLogSpot,prevSpot;			
	double termvolvol = 0.0;
	double spotvolvol = 0.0,z1=0.0;
	
	setUpNextTimeStep(idate);
	n_offset++;

	for ( ipath = 0 ; ipath < m_nPaths; ipath++ )
	{
		rand		=	GetNextRandom(ipath,n_offset);		
		vol			=	m_discreteVol[indexVec[ipath]];
		newLogSpot	=	getNewLogSpotIncrement(rand,indexVec[ipath],idate,ipath);
		prevSpot	=	(*m_pPath_array)[ipath][m_iasset][idate];


		(*m_pPath_array)[ipath][m_iasset][idate+1] = prevSpot*exp(newLogSpot);			
		saveVolofVolInformation(idate,ipath,vol,termvolvol,z,spotvolvol,z1);

	}			

	centerForwards(idate+1,m_iasset,m_termFwds[idate+1][0]);

	calculateVolofVolInformation(idate,termvolvol,z, spotvolvol, z1);	
	generateControlVariatePath(idate,calibflag,idum);
}	

/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: centerForwards
**	Returns: nothing
**	Action : centers forward values numerically
****************************************************************/

void MlEqMonteCarlo::centerForwards(int idate,int iasset,double fwd)
{
	if ( !m_centerFwds ){
		return;
	}

	double fac = 0.0;
	for ( int ipath = 0 ; ipath < m_nPaths; ipath++ ){
		fac += (*m_pPath_array)[ipath][iasset][idate];
	}

	fac /= m_nPaths;
	fac = fwd/fac;

	for ( int ipath = 0 ; ipath < m_nPaths; ipath++ ){
		(*m_pPath_array)[ipath][iasset][idate] *= fac;
	}
}

/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: generateControlVariatePath
**	Returns: nothing
**	Action : executes one Monte Carlo timestep
****************************************************************/

void CForwardSkewMC::generateControlVariatePath(int idate,int calibflag,long& idum)// generates paths from idate->idate+1
{		

	if ( m_controlVariate == 0 || m_controlVariateCalcFlag == 0 )
		return;

	double rand,vol;
	int ipath;
	int iRandom=0;
	int n_offset = m_numberOfFactors*idate;
	double newSpot,prevSpot;			
	
	n_offset++;
	for ( ipath = 0 ; ipath < m_nPaths; ipath++ )
	{
		iRandom		=	ipath;
		rand		=	GetNextRandom(iRandom,n_offset);//MlEqMaths::gasdev(&idum);////MlEqMaths::gasdev(&idum);////;//GetNextRandom(iRandom);//MlEqMaths::gasdev(&idum);//GetNextRandom(iRandom);//randoms[iRandom++];// m_randomNumbers->throwNext(fillflag,useRandoms);//generateRandom(idate,idum,fillFlag);		
		iRandom++;

		prevSpot	=	m_controlVariateCurrentSpot[ipath];
		vol			=	m_VolAtm[idate];//m_controlVariateVols[idate-1];
		newSpot		=	prevSpot*m_fwds[idate][0]*exp(-0.5*vol*vol*m_dt[idate]+vol*m_sqrtdt[idate]*rand);

		m_controlVariateCurrentSpot[ipath] = newSpot;
	}			

}	



/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: calculateVolofVolInformation
**	Returns: nothing
**	Action : calculates spot vol of vol and term vol of vols empirically
****************************************************************/


void CForwardSkewMC::calculateVolofVolInformation(int idate,double termvolvol,double z,double spotvolvol,double z1)
{
	if ( !m_saveVolofVolInfo )
		return;

		termvolvol	/= m_nPaths;
		z			/= m_nPaths;

		spotvolvol  /= m_nPaths;
		z1			/= m_nPaths;

		double t = 0.0;
		for ( int i = 0 ; i <= idate;i++ ){
			t += m_dt[idate];
		}

		z = log(z/(termvolvol*termvolvol));
		if ( fabs( z ) < 1e-3 ){
			z = 0.0;
		}
		else if ( z < 0.0 ){
			throw(  "negative variance for vol of vol encountered" );
		}

		m_termvolvol[idate] = sqrt(z/t);

		z1 = z1-spotvolvol*spotvolvol;

		if ( fabs( z1 ) < 1e-3 ){
			z1 = 0.0;
		}
		else if ( z1 < 0.0 ){
			throw(  "negative variance for vol of vol encountered" );
		}

		m_spotvolvol[idate] = sqrt(z1)/m_sqrtdt[idate];
		
}	
			

/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: saveVolofVolInformation
**	Returns: nothing
**	Action : 
****************************************************************/


void CForwardSkewMC::saveVolofVolInformation(int idate,int ipath,double vol,double& termvolvol,double&z,double& spotvolvol,double&z1)
{

	if ( !m_saveVolofVolInfo )
		return;

	double prevvol;

	if ( idate == 0 ){
		prevvol = m_initialShortVol;
	}
	else{
		prevvol = m_prevVol[ipath];
	}

	double x = log(vol/prevvol);
	termvolvol	+= vol;
	z			+= vol*vol;
	spotvolvol  +=	x;
	z1			+= x*x;

	m_prevVol[ipath] = vol;		
	m_saveDiscreteVols[idate] =  m_discreteVol[idate];

}

/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: getNewLogSpot
**	Returns: double
**	Action : returns new log of new spot value
****************************************************************/


double CForwardSkewMC::getNewLogSpotIncrement(double omega,int volindex,int idate,int ipath)
{
	double logFwdFactor = getlogFwdFactors(idate,ipath) ;
	
	if ( m_hermiteSkew.getsize() == 0 || idate == 0 )
	{
		double guess = getGuess(omega,volindex);
		return guess+logFwdFactor;
	}

	double logspot = m_currentHermiteSkews[volindex].getStateVar(omega)+logFwdFactor;
	return logspot;
}

		

/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: getGuess
**	Returns: double
**	Action : returns cummulative probability for an input logK
****************************************************************/


double CForwardSkewMC::getGuess(double omega,int volindex)
{
	int index=0;
	MlEqMaths::locate(m_omegaGridPoints[volindex],omega,index);
	index++;
	if ( index == m_omegaGridPoints[volindex].getsize() )
		index--;

	double res;
	res = m_spotGrid[index][volindex];//sos

	if ( index < m_spotGrid.rows()-1 )
	{
		res = res + ( m_spotGrid[index+1][volindex]-res )/( m_omegaGridPoints[volindex][index+1]
				- m_omegaGridPoints[volindex][index])*(omega-m_omegaGridPoints[volindex][index]);
	}
	else
	{
		res = res + ( m_spotGrid[index-1][volindex]-res )/( m_omegaGridPoints[volindex][index-1]
				- m_omegaGridPoints[volindex][index])*(omega-m_omegaGridPoints[volindex][index]);
	}
	
	return res;
}	



/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: calibrate
**	Returns: double
**	Action : calibrates forward vols to the cliquet market
****************************************************************/

void CForwardSkewMC::calibrate(int idate,CVector& currentAtmVols,double factor,int calibflag)
{

// calibflag = 

// 0: no calibration
// 1: calibration to the cliquet market
// 2: calibration to first slice only
// 3: calibration to first slice and at the money termstructure

	CVector tmpVols(currentAtmVols.getsize()+1);;
	
	for ( int i=0; i< currentAtmVols.getsize(); i++ )
	{
		tmpVols[i+1] = currentAtmVols[i];
	}


	hpsort((long)m_nPaths, (double*)tmpVols);

	m_discreteVolProb.resize(m_numberVolStates);
	m_discreteVol.resize(m_numberVolStates);
	m_discreteVolArray.resize(m_numberVolStates);

	double    cL			= m_ForwardVolcurve[idate].cL();
	double    cR			= m_ForwardVolcurve[idate].cR();
	int       addTanhWing	= m_ForwardVolcurve[idate].addTanhWing();
	int		  ypower		= m_ForwardVolcurve[idate].yPower();
	
	
	CVector strikes(m_vols[idate].getsize());
		
	for (int i =0 ; i < strikes.getsize(); i++ )
	{
		//the following are calibration strikes in terms of forward
		strikes[i] = exp(m_normalizedStrikes[idate][i]*m_VolAtm[idate]*m_sqrtdt[idate]);
	}
	
	
	CVector values(m_normalizedStrikes[idate].getsize());
	
	double volmultiplyer = 1.0;
	
	m_bareVolAtm  = m_bareForwardVolcurve[idate].getValue(0.0);
	
	calibrateHelper(values,tmpVols,idate,volmultiplyer,strikes,calibflag);
	
	if ( calibflag == 0 || ( calibflag == 2 && idate >  0 ) || ( calibflag == 3 && idate >  0 ) )
	{	
		return;
	};	
	
	CVector bareVols(m_numberVolStates);
	bareVols = m_vols[idate];
	CVector xstrikes(values.getsize());
	
	int max_tries = 200;
	
	double accuracy = 0.0005;
	
	int iter,flag,j;
	double volatm,volold,volnew;
	
	
	for ( iter = 0; iter < max_tries; iter++ )
	{		
		
		flag = 0;
		for (int i = 0; i < strikes.getsize(); i++ )
		{		
			if ( fabs(values[i]) > accuracy )
			{	
				flag = 1;
				break;
			}	
		}		
			
		if ( flag == 1 )
		{	
			for ( j = 0 ; j < values.getsize(); j++ )
			{
				bareVols[j] += values[j]; 
			}
			
			int method = 0;	
			volatm =  MlEqMaths::linearInterp(m_normalizedStrikes[idate],bareVols,0.0,method);
			
			for (int i = 0 ; i < strikes.getsize(); i++ )
			{
				xstrikes[i] = log(strikes[i])/(volatm*m_sqrtdt[idate]);
			}
			
			volold = m_bareForwardVolcurve[idate].getValue(0.0);
			
			m_bareForwardVolcurve[idate].setInputData( xstrikes, bareVols, ypower );
			m_bareForwardVolcurve[idate].setAlgPars(addTanhWing, cL, cR );
			m_bareForwardVolcurve[idate].buildCurve( );
			
			m_bareVolAtm = m_bareForwardVolcurve[idate].getValue(0.0);
			
			volnew = m_bareForwardVolcurve[idate].getValue(0.0);
			
			volmultiplyer *= volnew/volold;
			
			m_bareVolAtm = volatm;
			
			calibrateHelper(values,tmpVols,idate,volmultiplyer,strikes,calibflag);
						
		}	
		else	
			break;
		
	}		
		

	if ( iter == max_tries )
	{
		throw(  "calibration failed" );
	}

	m_volmultiplyer = volmultiplyer;
	m_bareVolAtm  = m_bareForwardVolcurve[idate].getValue(0.0);
	
		
}



/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: calibrateHelper
**	Returns: double
**	Action : calculates current vol minus mkt vols
****************************************************************/


void CForwardSkewMC::calibrateHelper(CVector& values,CVector& simulatedVols,int idate,double volmultiplyer,CVector& strikes,int calibflag)
{
		
	int nbins = m_numberVolStates*100;

	int jpath = m_nPaths/nbins;
	
	CVector discreteVolProb(nbins);	
	CVector discreteVol(nbins);

	CVector xdata(nbins);
	CVector ydata(nbins);
		
	int index;
	int i;
	
	double test = 0.0;
	

	for ( i = 1 ; i < nbins; i++ )
	{	
		index = i*jpath;
			
		xdata[i] = (double)index/(double)m_nPaths;
		ydata[i-1] = simulatedVols[index]*volmultiplyer;
		
		test += xdata[i]-xdata[i-1];
		
		discreteVolProb[i-1] = xdata[i]-xdata[i-1];
	}	
	
	discreteVolProb[i-1] = 1.0-xdata[i-1];
	ydata[i-1] = simulatedVols[simulatedVols.getsize()-1]*volmultiplyer;
		
	test = 0.0;
	for ( i = 0 ; i < discreteVolProb.getsize(); i++ )
	{	
		test+=discreteVolProb[i];
	}

	test=0.0;
	for ( i = 0 ; i < nbins; i++ )
	{	
		if ( i == 0 )
		{
			discreteVol[i] = ydata[i]; 
		}
		else
		{
			discreteVol[i] = (ydata[i]+ydata[i-1])/2.0;
		}
		
		test += discreteVolProb[i]*discreteVol[i];
	}	
		
		
//  corsen grid:	

	int inverseCoarseness = nbins/m_numberVolStates;
	
	double bucketProb;
	double bucketVol,z;
	int k =0;
	int j=0;
	test = 0.0;

	bucketProb	= 0.0;
	bucketVol	= 0.0;	

	for ( i = 0 ; i < nbins; i++ )
	{	
		
		bucketProb	+=	discreteVolProb[i];
		bucketVol	+=	discreteVolProb[i]*discreteVol[i];
		k++;
			
		if ( k == inverseCoarseness || i == nbins-1 )
		{

			m_discreteVolProb[j] = bucketProb;
			z = bucketVol/bucketProb;
			m_discreteVol[j]  = z;

			test += m_discreteVol[j]*m_discreteVolProb[j];

			m_discreteVolArray[j] = ydata[i-1];
			j++;
			k=0;
			bucketProb	= 0.0;
			bucketVol	= 0.0;			

		}
	}	
	

////////////////////////////
	

  if ( idate>0 && calibflag ==2 || idate>0 && calibflag == 3 || calibflag == 0 )
    return;

//  start sweep calibration
		
	double lower_x,upper_x;
	
	int max_tries= 10;
	
	double accuracy = 0.0005;

	double rvols,val;

	lower_x = 0.01;
	upper_x = 1.5;
	
	double acc = 1e-15;//sos1e-9;
	
	for ( i = 0; i < strikes.getsize(); i++ )
	{		
				
		val = calibrateFcn(strikes[i],idate);

		rvols =  MlEqBSImpliedVol(
							 val,	
							 1,		
							 m_dt[idate],		
							 strikes[i],	
							 1.0,
							 1 ,					
							 eROOTFIND_BRENT_GROWRANGE,	
							 acc,		
							 lower_x,	
							 upper_x);


		values[i] = m_vols[idate][i]-rvols;

	}

}


/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: calibrateHelper
**	Returns: double
**	Action : calculates current vol minus mkt vols
****************************************************************/

double CForwardSkewMC::calibrateFcn(double strike,int idate)
{		
	
//  assumes strike in terms of the forward
	
	double fwd = 1.0;//sos m_fwdFactors[idate];
	double mat = m_dt[idate];
	double val = 0.0;
	double test=0.0;
	double sqdt = sqrt(mat);
	double vol,x,yvol;	
	double volatm = m_bareForwardVolcurve[idate].getValue(0.0);
	double excessVolOfVol; 

	for ( int istate = 0 ; istate < m_numberVolStates; istate++ )
	{					
		double prob		= m_discreteVolProb[istate];
		vol				= m_discreteVol[istate];
		x				= log(strike)/(vol*sqdt);
		excessVolOfVol  = getExcessVolOfVol(idate,x);

//		yvol = sqrt(m_bareForwardVolcurve[idate](x))/volatm*vol;sos
		yvol			= 	getNewVolSkew(m_bareForwardVolcurve[idate], x, vol, volatm,excessVolOfVol);
		test			+= prob*vol;
		val				+= prob*Bs(fwd,yvol,mat,strike,1.0,1.0);			
	}			
	return val;
}

/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: GetNextRandom
**	Returns: double
**	Action : generates MonteCarlo paths
****************************************************************/

double MlEqMonteCarlo::GetNextRandom(int& iRandom,int& n_offset)
{	
	if ( m_randomsAreExternal  ){
		return (*m_pallMultiRandoms)[iRandom][m_iasset][n_offset];//[ipath][iasset][m_numberOfFactors*idate]
	}

	if ( m_allRandoms.rows() == 0 ){
		return m_oneStepRandoms[n_offset++];
	}
	else{
		return m_allRandoms[iRandom][n_offset];
	}

}

/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: setExternalRandoms
**	Returns: nothing
**	Action : don't do this at home
****************************************************************/

void  MlEqMonteCarlo::setExternalRandoms(int iasset,GMatrix < CVector >&	m_allMultiRandoms,cTensor& Path_array)
{

//  multiasset forward skew montecaltos can be constructed as a collection of single
//  CForwardSkewMC like this one. In order to do that externally created (correlated)
//  random numbers need to be passed into here

	m_iasset = iasset;
	m_pallMultiRandoms = &m_allMultiRandoms;//don't try this at home
	m_randomsAreExternal = 1;
	m_pPath_array = &Path_array;// don't do this at home
}


/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: GenerateRandomNumbers
**	Returns: void
**	Action : generates randomNumbers
****************************************************************/


void  MlEqMonteCarlo::GenerateRandomNumbers(int idate)
{	
	if ( m_randomsAreExternal  ){
		return;
	}

/*	int i;	
	if ( m_randomNumberFlag == 0 )
	{	
		int fillflag=0,useRandom=0;
		if ( idate == 0 ){
			m_oneStepRandoms.resize(m_numberOfFactors*m_nPaths);
		}

		for ( i = 0 ; i < m_oneStepRandoms.getsize(); i++ ){	
			m_oneStepRandoms[i] = m_randomGenerator->throwNext(fillflag,useRandom);
		}	
		return;
	}	
*/
	
//	generate all random numbers here

	if ( idate == 0 && m_allRandoms.rows() == 0 )
	{
		m_allRandoms.resize(m_nPaths,m_numberOfFactors*m_nDates);
		for (int ipath = 0 ; ipath <m_nPaths; ipath++ ){	
			m_randomGenerator->generateRandoms(m_allRandoms[ipath],ipath);
		}	
	}	

	return;
}		

/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: GenerateRandomNumbers
**	Returns: void
**	Action : generates correlated randomNumbers
****************************************************************/


void  MlEqMonteCarlo::GenerateRandomNumbers(Cholesky& cholesky)
{

	int iasset,ipath;	
	
//	generate all random numbers here

	if ( m_allMultiRandoms.rows() == 0 )
	{
		m_allMultiRandoms.resize(m_nPaths,m_nAssets);
		for ( ipath = 0 ; ipath < m_nPaths; ipath++ )
		{
			for ( iasset = 0 ; iasset < m_nAssets; iasset++ ){
				m_allMultiRandoms[ipath][iasset].resize(m_nDates*m_numberOfFactorsPerAsset[iasset]);
			}
		}

		CVector rndTemp;
		int ndim = 0;
		for ( iasset = 0 ; iasset < m_nAssets; iasset++ ){
			ndim += m_numberOfFactorsPerAsset[iasset];
		}
//		ndim *= (m_nDates-1);
		ndim *= m_nDates;

		for (ipath = 0 ; ipath <m_nPaths; ipath++ ){	
			m_randomGenerator->generateRandoms(m_allMultiRandoms[ipath],m_numberOfFactorsPerAsset,cholesky,ipath,rndTemp);
		}	

	}	

	return;
}		




/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: generatePath
**	Returns: double
**	Action : generates MonteCarlo paths
****************************************************************/

void CForwardSkewMC::GeneratePath(int calibflag)
{
	if ( m_nDates == 0 ){
		return;
	}

	int idate;	
	initializeRandomGenerator(m_numberOfFactors);
	long idum = m_seed;
	for ( idate = 0; idate < m_nDates; idate++ )
	{
		GenerateRandomNumbers(idate);
		initIRStateVars(idate);
		if ( idate == 0 )
		{
			m_currentVols.resize(m_nPaths);
			m_pcurrentVols		= &m_currentVols;
			m_pnewcurrentVols	= &m_currentVols;
		}
		generatePath(idate,calibflag,idum);
	}

	m_pcurrentVols		= NULL;
	m_pnewcurrentVols	= NULL;
	m_psavecurrentVols	= NULL;
	m_currentVols.resize(0);
	m_allRandoms.resize(0,0);
	m_oneStepRandoms.resize(0);
}



/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: GetImpliedVol
**	Returns: double
**	Action : gets impied vol form option price
****************************************************************/

double MlEqMonteCarlo::GetImpliedVol(int iStart,int iEnd,double optionPrice,double strike,int iasset)
{
	int cp = 1;
	int rootfind_flag = 8;
	double accuracy   =1e-8;//1e-5; 
	double lower_x    =0.0; 
	double upper_x    =1.5;
	double discount_factor = 1.0;
	double fwd,dt;
	double result;
	
	dt = 0.0;
	fwd = 1.0;
	for ( int i = iStart; i < iEnd; i++ )
	{
		dt		+=	m_dt[i];
		fwd		*=	m_fwds[i][iasset];
	}
	

	result =  MlEqBSImpliedVol(optionPrice,fwd,dt,strike,discount_factor,cp,
								rootfind_flag ,accuracy ,lower_x ,upper_x);
	
	return result;
}

/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: GetImpliedVol
**	Returns: double
**	Action : gets impied vol form option price
****************************************************************/

void MlEqMonteCarlo::averageMonteCarloResults(CVector& results)
{
	if ( m_nDates == 0 ){
		return;
	}
	m_randomGenerator->average(results,m_nPaths);
}



/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: GetImpliedVol
**	Returns: double
**	Action : gets implied vol form option price
****************************************************************/

void MlEqMonteCarlo::averageMonteCarloResults(double& result)
{
	CVector val(1);
	val[0] = result;
	m_randomGenerator->average(val,m_nPaths);
	result = val[0];
}

/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: GetImpliedVol
**	Returns: double
**	Action : gets implied vol form option price
****************************************************************/

double CForwardSkewMC::getlogFwdFactors(int idate,int ipath)
{
	return m_logFwdFactors[idate];
}


/***************************************************************
**	Class  : CForwardSkewMC 
**	Routine: GetImpliedVol
**	Returns: double
**	Action : gets impied vol form option price
****************************************************************/

void MlEqMonteCarlo::completePayoff(double& payoff,int ipath)
{
	m_randomGenerator->completeValue(payoff,ipath);
}


///////////////////////////////////////////////////////////////

void QuickMultiMc::choleskySetup()
{	
	int idate,iasset,jasset;
	double x;
	int nassets = m_nassets;
	
	CMatrix covariance(nassets,nassets);
	CVector diag(nassets);
	m_choleskyVols.resize(m_mDates,nassets,nassets);
		
	for ( idate = 0 ; idate < m_mDates; idate++ )
	{		
		//	calculate forward covariance		
		for ( iasset = 0 ; iasset < nassets; iasset++ )
			for ( jasset = 0 ; jasset < nassets; jasset++ )
			{	
				//set to 0
				m_choleskyVols[idate][iasset][jasset] = 0.0;
		
				double term = m_times[idate+1]-m_times[idate];
		
				double iVol = m_vols[iasset];
				double jVol = m_vols[jasset];
				double corr = getCorrelation(iasset,jasset);
				x = iVol * jVol * corr * term;
				
				// implementation only for constant correlations
				if ( iasset == jasset ){
					diag[iasset] = x;
				}
				covariance[iasset][jasset] = x ;
			}
				
		// start with cholesky decomposition		
		MlEqMaths::choldc(covariance,diag);
						
		// remember lower triangle of covariance has been changed !
		for ( iasset = 0 ; iasset < nassets; iasset++ )
		{				
			for ( jasset = 0 ; jasset < nassets; jasset++ )
			{				
				if ( iasset == jasset ){
					m_choleskyVols[idate][iasset][jasset]  += diag[iasset];
				}
				else if ( jasset < iasset ){
					m_choleskyVols[idate][iasset][jasset]  += covariance[iasset][jasset];
				}
			}	
		}		
	}			
}				

void QuickMultiMc::initialize(CVector& times,CMatrix& inputs,CVector& params)
{
	int i;

	m_times = times;

	for ( i = 0 ; i < inputs.rows(); i++ )
	{
		m_fwdRates[i] = inputs[i][0];
		m_blends[i]   = inputs[i][1];
		m_vols[i]	   = inputs[i][2];//[iasset]
	}

	m_nDates = times.getsize();
	m_mDates = m_nDates-1;
	m_dt.resize(m_mDates);

	for ( i = 0 ; i < m_mDates; i++ )
	{
		m_dt[i] = m_times[i+1]-times[i];
	}

	m_nassets						= m_fwdRates.getsize();
	m_seed							= params[0];
	m_nPayoffDates					= params[1];
	m_numberStepsBetweenSettings	= params[2];	;
	m_nPaths						= m_nPayoffDates*m_numberStepsBetweenSettings;

}

void QuickMultiMc::getNewSpots(CVector& newSpots,int idate,CVector& currentSpots,CVector& randoms)
{
	int iasset,ifactor;
	double x,vol;
	
	int nassets = m_nassets;
	for ( iasset = 0; iasset < nassets; iasset++ )
	{	
		x = 0.0;
		for ( ifactor = 0; ifactor < nassets; ifactor++ )
		{	
			double cvols = m_choleskyVols[idate][iasset][ifactor]; //(idate,iasset,ifactor);
			double random = randoms[ifactor];
			x += cvols * random;
		}
		
		vol = m_vols[iasset];
		newSpots[iasset] = currentSpots[iasset]*exp(m_fwdRates[iasset]*m_dt[idate])
						   *exp(-0.5*vol*vol*m_dt[idate] + x );
				
	}									

}



void QuickMultiMc::simulate(CVector& results,CVector& params)
{

	CVector currentSpots(m_nassets),newSpots(m_nassets),randoms(m_nassets);

	CMatrix pathArray(m_nDates,m_nPaths);

	int iasset,idate,i,ipath;
	for ( iasset = 0 ; iasset < m_nassets; iasset++ )
	{
		currentSpots[iasset]	= 100.0;
		pathArray[0][iasset]	= 100.0;
	}

	long idum = m_seed;

	double value = 0.0;

	for ( ipath = 0 ; ipath < m_nPaths; ipath++ )
	{
		for ( idate = 0 ; idate < m_nDates; idate++ )
		{
			for ( i = 0 ; i < m_nassets; i++ )
			{
				randoms[i] = MlEqMaths::gasdev(&idum);
			}

			getNewSpots(newSpots,idate,currentSpots,randoms);

			pathArray[idate] = newSpots;
		}

		value += payoff(pathArray,params);

		currentSpots = newSpots;

	}

	value /= value/m_nPaths;
	results[0] = value;
}



double QuickMultiMc::payoff(CMatrix& pathArray,CVector& params)
{

	int idate,iasset,istep;
	double barrier = params[0];
	double Coupon = params[1];

	
	int hitbarr;
	CVector barrHistory(m_nPayoffDates);

	int k = 1;
	for ( idate = 0 ; idate < m_nPayoffDates; idate++ )
	{
		hitbarr = 0;
		for ( istep = 0 ; istep < m_numberStepsBetweenSettings; istep++ )
		{
			if ( hitbarr == 0 )
			{
				for ( iasset = 0 ; iasset > m_nassets; iasset++ )
				{
					if ( pathArray[k][iasset] < barrier )
					{
						hitbarr = 1;
						break;
					}
				}
			}

			k++;
		}

		if ( hitbarr == 1 )
		{
			barrHistory[idate] = 1;
		}

	}

	double payout = 0.0;
	for ( idate = 0 ; idate < m_nPayoffDates; idate++ )
	{
		double levFactor = barrHistory[idate];

		for ( int jdate = idate-1; jdate >= 0 ; jdate--)
		{
			if ( barrHistory[jdate] )
			{
				break;
			}
			else
			{
				levFactor++;
			}
		}

		payout += levFactor*Coupon;

	}

	return payout;	

}


/////////////////////////////////////////////////////////////////////////////
//
//
//				CMonteCarloHelper
//
//
/////////////////////////////////////////////////////////////////////////////


/***************************************************************
**	Class  : CMonteCarloHelper 
**	Routine: operator=
**	Returns: nothing
**	Action : assignment operator
****************************************************************/

const CMonteCarloHelper& CMonteCarloHelper::operator=(const CMonteCarloHelper& mch)
{
	m_PayoffDates = mch.m_PayoffDates;
	m_bUseToday = mch.m_bUseToday;
	m_vectorFixedSpots = mch.m_vectorFixedSpots;
	m_vectorMCDates = mch.m_vectorMCDates;
	m_mapFromPayoff	= mch.m_mapFromPayoff;

	return *this;
}


/***************************************************************
**	Class  : CMonteCarloHelper 
**	Routine: CMonteCarloHelper
**	Returns: nothing
**	Action : constructor
****************************************************************/


void CMonteCarloHelper::initialize(const std::vector<long>& Dates, const CMatrix& Fixings, long nToday)
//	Dates - the product payoff is a function of these dates
//	Fixings - these are the spot values corresponding to elements in Dates before nToday
//	Today's date
{		

//	m_PayoffDates	= Dates;	

//	if ( m_PayoffDates.size() == 0 ){
//		throw("no payoff dates have been entered into monteCarlo helper");}

	m_vectorMCDates.clear();	
	// nToday is always needed and is always the first element in m_vectorMCDates
	m_vectorMCDates.push_back(nToday);
	
	// We need to determine the number of elements in Dates that we need.
	// This is the number of elements in Dates occurring after nToday.
	m_bUseToday = false;
	long nLow = 1;
	long nHigh = Dates.size() - 1;
	long nElements = -1;
	while (nLow <= nHigh){
		long nMid = (nLow + nHigh) / 2;
		if (nToday >= Dates[nMid - 1] && nToday < Dates[nMid]){
			// done.
			nElements = Dates.size() - nMid;
			if (nToday == Dates[nMid - 1]) m_bUseToday = true;			// this is the direct hit condition			
			break;
		} else if (nToday >= Dates[nMid]){
			nLow = nMid + 1;
		} else {
			nHigh = nMid - 1;
		}
	}
	if (nElements < 0){
		nElements = Dates[Dates.size() - 1] <= nToday ? 0 : Dates.size();
		if (nToday == Dates[Dates.size() - 1]) m_bUseToday = true;		// this is the direct hit condition on the last element
	}

	// form m_vectorMCDates	
	for (long nElement = 1; nElement <= nElements; nElement++){
		m_vectorMCDates.push_back(Dates[Dates.size() - nElements + nElement - 1]);
	}

	// get m_vectorFixedSpots

	int i;
	vector < vector < double > > fixedSpots;
//	STLVectorVectorTransposeFromCMatrix(fixedSpots,Fixings);
	STLVectorVectorFromCMatrix(fixedSpots,Fixings);



	for ( i = 0 ; i < fixedSpots.size(); i++ ){
		fixedSpots[i].resize(Dates.size() - nElements - (m_bUseToday ? 1 : 0));
	}

	CMatrixFromSTLVectorVector(m_vectorFixedSpots,fixedSpots);

	int iDate = 0;
		// put fixings in first
	while (iDate < m_vectorFixedSpots.cols()){
		iDate++;
	}


	mapFromPayoff(Dates, Fixings, nToday,Dates);


}

/***************************************************************
**	Class  : CMonteCarloHelper 
**	Routine: CMonteCarloHelper
**	Returns: nothing
**	Action : constructor
****************************************************************/

void CMonteCarloHelper::initialize(const std::vector<long>& Dates, const CMatrix& Fixings, long nToday, const std::vector<long>& PayoffDates)
{

	initialize(Dates,Fixings,nToday);
	mapFromPayoff(Dates, Fixings, nToday,PayoffDates);

}

/***************************************************************
**	Class  : CMonteCarloHelper 
**	Routine: CMonteCarloHelper
**	Returns: nothing
**	Action : constructor
****************************************************************/

void CMonteCarloHelper::mapFromPayoff(const std::vector<long>& Dates, const CMatrix& Fixings, long nToday, const std::vector<long>& PayoffDates)
{

	m_PayoffDates = PayoffDates;

	if ( m_PayoffDates.size() == 0 ){
		throw("no payoff dates have been entered into monteCarlo helper");}


	int iDate = 0;
	int iasset = 0;
		// put fixings in first
	while (iDate < m_vectorFixedSpots.cols()){
		iDate++;
	}

	m_mapFromPayoff.resize(PayoffDates.size()-iDate);

	for ( int i = iDate ; i < PayoffDates.size(); i++ )
	{
		for ( int j = 0; j < m_vectorMCDates.size(); j++ )
		{
			if ( PayoffDates[i] == m_vectorMCDates[j] ){
				m_mapFromPayoff[i-iDate] = j;
				break;
			}

			if ( j == m_vectorMCDates.size() ){
				throw("indexing error in seeting up MC dates");
			}
		}
	}

}

/***************************************************************
**	Class  : CMonteCarloHelper 
**	Routine: GetMCDates
**	Returns: nothing
**	Action : 
****************************************************************/


int CMonteCarloHelper::mapFromPayoff(int idate)
{
// idate is a simulation date index
	return m_mapFromPayoff[idate];
}

/***************************************************************
**	Class  : CMonteCarloHelper 
**	Routine: GetMCDates
**	Returns: nothing
**	Action : 
****************************************************************/

CMonteCarloHelper::CMonteCarloHelper(const std::vector<long>& Dates, const CMatrix& Fixings, long nToday)
{
	initialize(Dates, Fixings,nToday);
}

/***************************************************************
**	Class  : CMonteCarloHelper 
**	Routine: GetMCDates
**	Returns: nothing
**	Action : 
****************************************************************/

const std::vector<long>& CMonteCarloHelper::GetMCDates(void) 
{
	return m_vectorMCDates;
}	


/***************************************************************
**	Class  : CMonteCarloHelper 
**	Routine: GetMCDates
**	Returns: nothing
**	Action : 
****************************************************************/

long CMonteCarloHelper::GetMCDate(int idate) 
{
	if ( idate < 0 || idate>= m_vectorMCDates.size() ){
		throw("incorrect date index invoked");
	}
	return m_vectorMCDates[idate];
}	

/***************************************************************
**	Class  : MlEqMonteCarlo
**	Routine: GetMCDates
**	Returns: nothing
**	Action : 
****************************************************************/

const std::vector<long>&	MlEqMonteCarlo::GetMCDates(void) 
{
	return m_mMCHelper.GetMCDates();
}

/***************************************************************
**	Class  : MlEqMonteCarlo
**	Routine: GetMCDates
**	Returns: nothing
**	Action : 
****************************************************************/

long	MlEqMonteCarlo::GetMCDate(int idate) 
{
	return m_mMCHelper.GetMCDate(idate);
}


/***************************************************************
**	Class  : CMonteCarloHelper 
**	Routine: GetMCDates
**	Returns: nothing
**	Action : 
****************************************************************/

const std::vector<long>& CMonteCarloHelper::GetPayoffDates(void) const
{
	return m_PayoffDates;
}	


/***************************************************************
**	Class  : CMonteCarloHelper 
**	Routine: postInitialize
**	Returns: nothing
**	Action : 
****************************************************************/

void CMonteCarloHelper::postInitialize(std::vector<long>&	PayoffDates)
{
	m_PayoffDates = PayoffDates;
}



/***************************************************************
**	Class  : CMonteCarloHelper 
**	Routine: CreatePayoffPath
**	Returns: nothing
**	Action : 
****************************************************************/
	
	
bool CMonteCarloHelper::CreatePayoffPath(const CMatrix& MCSpots, CMatrix& PayoffSpots,CVector& payoutPathDiscounts,CVector& simulatedDiscounts)
{	
	for (int iasset = 0; iasset < MCSpots.rows(); iasset++){
		int iDate = 0;
		// put fixings in first
		while (iDate < m_vectorFixedSpots.cols())
		{
			PayoffSpots[iasset][iDate]			= m_vectorFixedSpots[iasset][iDate];
			payoutPathDiscounts[iDate]  = 0.0;// if setting is in the past: don't need discounts any more
			iDate++;
		}
		
/*		// now put in MCSpots
		int nIndex = 0;
		if (!m_bUseToday) nIndex++;		// i.e. miss out the first monte carlo spot
		for (; nIndex < MCSpots.cols(); nIndex++)
		{
			PayoffSpots[iasset][iDate]		= MCSpots[iasset][nIndex];
			payoutPathDiscounts[iDate]		= simulatedDiscounts[nIndex];
			iDate++;

		}
*/

		for (int n  = 0; n < PayoffSpots.cols()-iDate; n++)
		{
			PayoffSpots[iasset][n+iDate]	= MCSpots[iasset][m_mapFromPayoff[n]];
			payoutPathDiscounts[n+iDate]	= simulatedDiscounts[m_mapFromPayoff[n]];
		}

	}
	
	return true;
}	
	

/***************************************************************
**	Class  : CMonteCarloHelper 
**	Routine: CreatePayoffPath
**	Returns: nothing
**	Action : 
****************************************************************/
	
	
bool CMonteCarloHelper::CreatePayoffPath(const CMatrix& MCSpots, CMatrix& PayoffSpots)
{	
	for (int iasset = 0; iasset < MCSpots.rows(); iasset++){
		int iDate = 0;
		// put fixings in first
		while (iDate < m_vectorFixedSpots.cols()){
			PayoffSpots[iasset][iDate]			= m_vectorFixedSpots[iasset][iDate];
			iDate++;
		}
		
		// now put in MCSpots
/*		int nIndex = 0;
		if (!m_bUseToday) nIndex++;		// i.e. miss out the first monte carlo spot
		for (; nIndex < MCSpots.cols(); nIndex++){
			PayoffSpots[iasset][iDate++]		= MCSpots[iasset][nIndex];

		}
*/

		for (int n  = 0; n < PayoffSpots.cols()-iDate; n++){
			PayoffSpots[iasset][n+iDate]	= MCSpots[iasset][m_mapFromPayoff[n]];
		}


	}
	
	return true;
}	
	
/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: simulate
**	Returns: nothing
**	Action : 
****************************************************************/

void MlEqMonteCarlo::simulate(CMatrix& result,product& deriv)
{
//	generatePaths();

	CMatrix payoffs;
	CMatrix pathArray;
	CVector discounts;

	deriv.setUp(payoffs,*this);

	CMatrix payoffPathArray(m_nAssets, m_mMCHelper.GetPayoffDates().size());	
	CVector payoffDiscounts(m_mMCHelper.GetPayoffDates().size());

//	std::ofstream lolo("C:\\stuff\\MC.txt");

	if ( m_nDates == 0 )
	{
		m_nPaths = 1;
	}
	for (int ipath = 0; ipath < m_nPaths; ipath++ )
	{	
		m_currentPath	= ipath;
		pathArray		= GetPathArray(ipath);
		discounts		= GetDiscounts(ipath);

		createPayoutPath(payoffPathArray,pathArray,payoffDiscounts,discounts);

		deriv.payout(payoffs,payoffPathArray,payoffDiscounts,*this);		
		deriv.accruePayout(deriv.m_results,payoffs);

/*		if( (ipath+1) % 1000 == 0)
		{
			lolo << ipath+1 << '\t' ;
			
			int r = deriv.m_results.rows();
			int c = deriv.m_results.cols();
			for (int i = 0 ; i < r; i++ )
			{
				for (int j = 0 ; j < c; j++ ){
				lolo << deriv.m_results[i][j] / double(ipath) << '\t';
				}
			}
			lolo << '\n';
		}
*/	}	

	deriv.createResult(deriv.m_results,*this);

	result = deriv.m_results;
}


/*void CLocalVolMC::simulate(CMatrix& result,product& deriv)
{

	CMatrix payoffs;
	CMatrix pathArray;

	deriv.setUp(payoffs,*this);

	CMatrix payoffPathArray(m_nAssets, m_nDates);	
	CVector payoffDiscounts = m_discount;

	if ( m_nDates == 0 )
	{
		m_nPaths = 1;
	}
	for (int ipath = 0; ipath < m_nPaths; ++ipath )
	{	
		m_currentPath	= ipath;
		pathArray		= GetPathArray(ipath);

		deriv.payout(payoffs,pathArray,payoffDiscounts,*this);		
		deriv.accruePayout(deriv.m_results,payoffs);
	}	

	deriv.createResult(deriv.m_results,*this);

	result = deriv.m_results;
}

*/


/***************************************************************
**	Class  : MlEqMonteCarlo
**	Routine: simulate
**	Returns: nothing
**	Action : 
****************************************************************/

void MlEqMonteCarlo::basicInitializeMcHelper(product& deriv,CMatrix& fixings)
{
	long valDate = m_hDate->GetDate();
	vector<long > vDates(deriv.m_payoffDates.getsize());
	for ( int i = 0 ; i < vDates.size();i++){
		vDates[i] = deriv.m_payoffDates[i];
	}

	m_mMCHelper.initialize(vDates,fixings, valDate);

}

/***************************************************************
**	Class  : MlEqMonteCarlo
**	Routine: simulate
**	Returns: nothing
**	Action : 
****************************************************************/

double MlEqMonteCarlo::AmericanDigital(CVector& prices,int ipath,int idate,double DownBarrier,double finalSpot,
									 int iasset,int method,double stressForwardBarrierVol,
									 double stressForwardBarrierVolSlope,int npoints,double blend)
{
	throw("function not implemented");
	return -1;
}


/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: simulate
**	Returns: nothing
**	Action : 
****************************************************************/

void MlEqMonteCarlo::initializeMcHelper(product& deriv,CMatrix& fixings)
{
/*	long valDate = m_hDate->GetDate();
	vector<long > vDates(deriv.m_payoffDates.getsize());
	for ( int i = 0 ; i < vDates.size();i++){
		vDates[i] = deriv.m_payoffDates[i];
	}

	m_mMCHelper.initialize(vDates,fixings, valDate);
*/

	vector<long > additionalSimDates;
	long valDate = m_hDate->GetDate();
	
	if (valDate >= deriv.m_payoffDates[0] ){
		basicInitializeMcHelper(deriv,fixings);
		return;
	}

	if(deriv.m_payoffDates.getsize()==1){
		basicInitializeMcHelper(deriv,fixings);
		return;
	}


	if((deriv.m_payoffDates[1]-deriv.m_payoffDates[0])>(deriv.m_payoffDates[0]-valDate)){
	// fist period is less than the interval
		basicInitializeMcHelper(deriv,fixings);
		return;
	}


	
	long interval=(deriv.m_payoffDates[1]-deriv.m_payoffDates[0]); // interval in days
		
	int NumberOfAdditionalDates=(deriv.m_payoffDates[0]-valDate)/interval;//+1?
		
	additionalSimDates.resize(NumberOfAdditionalDates);
	
	long  idate=NumberOfAdditionalDates-1;
	long  cycleDate(deriv.m_payoffDates[0]);
	
	
	while(idate>=0)// >=0?
	{	
		cycleDate					-=	interval;
		additionalSimDates[idate--]	=	cycleDate;
	}
	

	initializeMcHelper(deriv,fixings,additionalSimDates);

}



/***************************************************************
**	Class  : MlEqMonteCarlo
**	Routine: simulate
**	Returns: nothing
**	Action : 
****************************************************************/

void MlEqMonteCarlo::initializeMcHelper(product& deriv,CMatrix& fixings,const vector<long >& additionalSimDates)
{
	if ( additionalSimDates.size() == 0 ){
		initializeMcHelper(deriv,fixings);
		return;
	}

	long valDate = m_hDate->GetDate();
	vector<long > vDates(deriv.m_payoffDates.getsize());

	for ( int i = 0 ; i < vDates.size();i++){
		vDates[i] = deriv.m_payoffDates[i];
	}

	std::vector<long> outPutSet;

	merge(outPutSet,vDates,additionalSimDates);

	m_mMCHelper.initialize(outPutSet,fixings, valDate,vDates);


	m_mMCHelper.postInitialize(vDates);

}



	
/////////////////////////////////////////////////////////////////////////////
//
//
//				Generic MonteCarlo
//
//
/////////////////////////////////////////////////////////////////////////////


/***************************************************************
**	Class  : MlEqMonteCarlo
**	Routine: createPayoutPath
**	Returns: nothing
**	Action : 
****************************************************************/



void MlEqMonteCarlo::createPayoutPath(CMatrix& payoutPathArray,CMatrix& simulatedPathArray)
{
	m_mMCHelper.CreatePayoffPath(simulatedPathArray,payoutPathArray);
}


/***************************************************************
**	Class  : MlEqMonteCarlo
**	Routine: createPayoutPath
**	Returns: nothing
**	Action : 
****************************************************************/



void MlEqMonteCarlo::createPayoutPath(CMatrix& payoutPathArray,CMatrix& simulatedPathArray,CVector& payoutPathDiscounts,CVector& simulatedDiscounts)
{
	m_mMCHelper.CreatePayoffPath(simulatedPathArray,payoutPathArray,payoutPathDiscounts,simulatedDiscounts);
}

/***************************************************************
**	Class  : MlEqMonteCarlo
**	Routine: initializeRandomGenerator
**	Returns: double
**	Action : generates MonteCarlo paths
****************************************************************/

void MlEqMonteCarlo::initializeRandomGenerator(int dimensionalityFactor)
{	
	if ( m_randomsAreExternal ){
		return;
	}

	int numberScenariosToBeStored=0,numberRandomsPerScenario=0;
	if ( m_randomNumberFlag == 0 )
	{
		long idum = m_seed;
		m_randomGenerator = new	randomGenerator(idum,numberScenariosToBeStored,
													numberRandomsPerScenario);
		return;
	}
	else if ( m_randomNumberFlag == 1 )
	{
		int dimension	= dimensionalityFactor*m_nDates;
		int skip		= 1000;
		int leap		= 1;

		m_randomGenerator = new Sobol(dimension,skip,leap,numberScenariosToBeStored,numberRandomsPerScenario);
	}
	else if ( m_randomNumberFlag >= 2 && m_randomNumberFlag <= 5)
	{
			int ndim = dimensionalityFactor*m_nDates;
			CUnitcube* D = new CUnitcube(ndim);
			//	select alpha
			//	CNiederreiteralpha *alpha = new CNiederreiteralpha;

			CAlpha* alpha;
			if ( m_randomNumberFlag == 2 ){
				alpha = new CPrimerootalpha;
			}
			else if ( m_randomNumberFlag == 3 ){
				alpha = new CNiederreiteralpha;
			}
			else if ( m_randomNumberFlag == 4 ){
				alpha = new CBakeralpha;
			}
			else if ( m_randomNumberFlag == 5 ){
				alpha = new CPrimerootalpha;
			}
			else{
				throw("randomNumberFlag must be between 1 and 5");
			}

			CArithmeticmeanweights* weights = new CArithmeticmeanweights;// weight
			int d = 1;// d >= r 
	//		CIdentityperiodization *periodization=new CIdentityperiodization;
			CBakerperiodization* periodization=new CBakerperiodization;
//			CTANHperiodization* periodization=new CTANHperiodization(1.0,1);
//			CDEperiodization* periodization=new CDEperiodization(1.0);
			CNTparameters* par = new CNTparameters(alpha,weights,periodization);
			CNTintegrator* NTintegrator= new CNTintegrator;

			m_randomGenerator = new Diophantine(ndim,m_nPaths,alpha ,
								weights,D,periodization,NTintegrator,par);
	}
	else if( m_randomNumberFlag == 6 )
	{
		int dimension	= dimensionalityFactor*m_nDates;
		m_randomGenerator = new CMersenneTwister(dimension);
	}

}

/***************************************************************
**	Class  : MlEqMonteCarlo
**	Routine: GetPathArray
**	Returns: nothing
**	Action : 
****************************************************************/

long MlEqMonteCarlo::GetNumberOfFutureDates()
{ 
//	int n = m_mMCHelper.GetPayoffDates().size() ;	
	long n = m_nDates;
	return n;
}


/***************************************************************
**	Class  : MlEqMonteCarlo
**	Routine: GetPathArray
**	Returns: nothing
**	Action : 
****************************************************************/


CMatrix& MlEqMonteCarlo::GetPathArray(int ipath)
{

	return (*m_pPath_array)[ipath];
}

/***************************************************************
**	Class  : MlEqMonteCarlo
**	Routine: Getdt
**	Returns: nothing
**	Action : Getdt[0] = t1-t0
****************************************************************/


double MlEqMonteCarlo::Getdt(int idate)
{
	return m_dt[idate];
}


/***************************************************************
**	Class  : MlEqMonteCarlo
**	Routine: Getdt
**	Returns: nothing
**	Action : Getdt[0] = t1-t0
****************************************************************/


double MlEqMonteCarlo::GetTime(int idate)
{
	return m_t[idate];
}


/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: Getdt
**	Returns: nothing
**	Action : Getdt[0] = t1-t0
****************************************************************/

double	MlEqMonteCarlo::GetTimeAtPayoffDate(int ipayoffDate)
{
	int n = m_mMCHelper.mapFromPayoff(ipayoffDate);
	return m_t[n];
}


/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: GetPathValue
**	Returns: double
**	Action : 
****************************************************************/


double	MlEqMonteCarlo::GetPathValue(int ipath,int idate,int iasset)
{
	return (*m_pPath_array)[ipath][iasset][idate];
}


/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: GetPathValue
**	Returns: double
**	Action : 
****************************************************************/

void MlEqMonteCarlo::initMC(const vector< long >& simDates)
{	
	CMatrix Fixings;
	m_mMCHelper.initialize(simDates, Fixings, m_hDate->GetDate());
}



/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: PutPathValue
**	Returns: double
**	Action : 
****************************************************************/

void MlEqMonteCarlo::PutPathValue(int ipath,int idate, int iasset, double fValue)
{
	//m_Path_array[ipath][iasset][idate] = fValue;
	throw "MlEqMonteCarlo::PutPathValue is not supported in this release";
}


/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: GetPathValue
**	Returns: double
**	Action : 
****************************************************************/


void CForwardSkewMC::generatePaths()
{
	if ( m_nDates == 0 ){
		return;
	}

	long idum = m_seed;
	GeneratePath(m_calibflag,idum,m_accuracy);
}

	
	

/***************************************************************
**	Class   : MlEqHullAndWhite
**	Function: B
**	Returns : double
**	Comment : this function usses indexStart and indexEnd
****************************************************************/

double MlEqHullAndWhite::B(double t,double T,int indexStart,int indexEnd)
{
	double res;
	int nsize = m_calibrationDates.getsize();
	int i = MlEqMaths::Min(indexStart,nsize-1);

	
	if ( indexStart < nsize-1 ){
		res = m_gamma_ti[i]/m_lambda*(exp(-t*m_lambda)-exp(-m_times[indexStart+1]*m_lambda));
	}
	else{
		res = m_gamma_ti[i]/m_lambda*(exp(-t*m_lambda)-exp(-T*m_lambda));
		return res;
	}

	for ( int i = indexStart+1 ; i < indexEnd; i++ ){
		res += m_gamma_ti[i]/m_lambda*(exp(-t*m_lambda)-exp(-m_times[i+1]*m_lambda));
	}

	if ( T > m_times[indexEnd] )
	{
		i = MlEqMaths::Min(indexEnd,nsize-1);
		res += m_gamma_ti[i]/m_lambda*(exp(-m_times[indexEnd]*m_lambda)-exp(-T*m_lambda));
	}

	return res;
}
/***************************************************************
**	Class   : MlEqHullAndWhite 
**	Function: B
**	Returns : double
**	Comment : this function also returns indexStart and indexEnd
****************************************************************/

double MlEqHullAndWhite::B(int& indexStart,int& indexEnd,long tDate,long TDate)
{

	if ( tDate > TDate ){
		throw("rateSettingDate cannot be greater than rateMaturityDate");
	}

	int nsize = m_calibrationDates.getsize();
	indexStart=0;
	MlEqMaths::Locate(m_calibrationDates,tDate,indexStart);
	indexEnd=0;
	MlEqMaths::Locate(m_calibrationDates,TDate,indexEnd);

	double t = m_hDate->GetYearFraction(tDate);
	double T = m_hDate->GetYearFraction(TDate);

	return B(t,T,indexStart,indexEnd);
}

/***************************************************************
**	Class   : MlEqHullAndWhite 
**	Function: B
**	Returns : double
**	Comment : this function also returns indexStart and indexEnd
****************************************************************/

void MlEqHullAndWhite::Locate(int& indexStart,int& indexEnd,long tDate,long TDate)
{

	if ( tDate > TDate ){
		throw("rateSettingDate cannot be greater than rateMaturityDate");
	}

	int nsize = m_calibrationDates.getsize();
	indexStart=0;
	MlEqMaths::Locate(m_calibrationDates,tDate,indexStart);
	indexEnd=0;
	MlEqMaths::Locate(m_calibrationDates,TDate,indexEnd);

}


/***************************************************************
**	Class   : MlEqHullAndWhite 
**	Function: B
**	Returns : double
**	Comment : 
****************************************************************/


double MlEqHullAndWhite::B(long tDate,long TDate)
{
	int indexStart,indexEnd;
	if ( tDate < m_hDate->GetDate() ){
		throw("B(.,.) function is tried to calculate starting from the past");
	}

	return B(indexStart,indexEnd,tDate,TDate);
}



/***************************************************************
**	Class   : MlEqHullAndWhite 
**	Function: integral_g
**	Returns : double
**	Comment : 
****************************************************************/


double MlEqHullAndWhite::integral_g(double t,double T,int indexStart,int indexEnd)
{

	double res;
	int nsize = m_calibrationDates.getsize();
	int i = MlEqMaths::Min(indexStart,nsize-1);

	if ( indexStart < nsize-1 ){
		res = m_beta_ti[i]/m_lambda*(exp(m_times[indexStart+1]*m_lambda)-exp(t*m_lambda));
	}
	else{
		res = m_beta_ti[i]/m_lambda*(exp(T*m_lambda)-exp(t*m_lambda));
		return res;
	}

	for ( int i = indexStart+1 ; i < indexEnd; i++ ){
		res += m_beta_ti[i]/m_lambda*(exp(m_times[i+1]*m_lambda)-exp(m_times[i]*m_lambda));
	}

	if ( T > m_times[indexEnd] )
	{
		i = MlEqMaths::Min(indexEnd,nsize-1);
		res += m_beta_ti[i]/m_lambda*(exp(T*m_lambda)-exp(m_times[indexEnd]*m_lambda));
	}

	return res;
}


/***************************************************************
**	Class   : MlEqHullAndWhite 
**	Function: integral_g
**	Returns : double
**	Comment : 
****************************************************************/

double MlEqHullAndWhite::integral_g(long t,long T)
{

	if ( t< m_hDate->GetDate() ){
		throw("integral_g(.,.) function is tried to calculate starting from the past");
	}

	int indexStart,indexEnd;
	return integral_g(indexStart,indexEnd,t,T);
}

/***************************************************************
**	Class   : MlEqHullAndWhite 
**	Function: integral_g
**	Returns : double
**	Comment : this function also returns indexStart and indexEnd
****************************************************************/


double MlEqHullAndWhite::integral_g(int& indexStart,int& indexEnd,long tDate,long TDate)
{

	if ( tDate > TDate ){
		throw("rateSettingDate cannot be greater than rateMaturityDate");
	}

	if ( tDate == TDate ){
		return 0.0;
	}

	int nsize = m_calibrationDates.getsize();
	indexStart=0;
	MlEqMaths::Locate(m_calibrationDates,tDate,indexStart);

	indexEnd=0;
	MlEqMaths::Locate(m_calibrationDates,TDate,indexEnd);

	double t = m_hDate->GetYearFraction(tDate);
	double T = m_hDate->GetYearFraction(TDate);

	return integral_g(t,T,indexStart,indexEnd);
}


/***************************************************************
**	Class   : MlEqHullAndWhite 
**	Function: integral_g2
**	Returns : double
**	Comment : this fuction uses indexStart and indexEnd
****************************************************************/

double MlEqHullAndWhite::integral_g2(double t,double T,int indexStart,int indexEnd)
{

	double res;
	int nsize = m_calibrationDates.getsize();
	int i = MlEqMaths::Min(indexStart,nsize-1);

	if ( indexStart < nsize-1 ){
		res = pow(m_beta_ti[i],2.0)/(2.0*m_lambda)*(exp(2.0*m_lambda*m_times[indexStart+1])-exp(2.0*m_lambda*t));
	}
	else{
		res = pow(m_beta_ti[i],2.0)/(2.0*m_lambda)*(exp(2.0*m_lambda*T)-exp(2.0*m_lambda*t));
		return res;
	}

	for ( int i = indexStart+1 ; i < indexEnd; i++ ){
		res += pow(m_beta_ti[i],2.0)/(2.0*m_lambda)*(exp(2.0*m_times[i+1]*m_lambda)-exp(2.0*m_times[i]*m_lambda));
	}

	if ( T > m_times[indexEnd] )
	{
		i = MlEqMaths::Min(indexEnd,nsize-1);
		res += pow(m_beta_ti[i],2.0)/(2.0*m_lambda)*(exp(2.0*T*m_lambda)-exp(2.0*m_times[indexEnd]*m_lambda));
	}
	return res;
}

/***************************************************************
**	Class   : MlEqHullAndWhite 
**	Function: integral_g2
**	Returns : double
**	Comment : this function also returns indexStart and indexEnd
****************************************************************/

double MlEqHullAndWhite::integral_g2(int& indexStart,int& indexEnd,long tDate,long TDate)
{

	if ( tDate > TDate ){
		throw("rateSettingDate cannot be greater than rateMaturityDate");
	}

	if ( tDate == TDate ){
		return 0.0;
	}

	double t	= m_hDate->GetYearFraction(tDate);
	double T	= m_hDate->GetYearFraction(TDate);
	int nsize	= m_calibrationDates.getsize();

	indexStart=0;
	MlEqMaths::Locate(m_calibrationDates,tDate,indexStart);

	indexEnd=0;
	MlEqMaths::Locate(m_calibrationDates,TDate,indexEnd);

	return integral_g2(t,T,indexStart,indexEnd);

}


/***************************************************************
**	Class   : MlEqHullAndWhite 
**	Function: integral_g2
**	Returns : double
**	Comment : 
****************************************************************/

double MlEqHullAndWhite::integral_g2(long t,long T)
{

	int indexStart,indexEnd;
	if ( t < m_hDate->GetDate() ){
		throw("integral_g2(.,.) function is tried to calculate starting from the past");
	}

	return integral_g2(indexStart,indexEnd,t,T);
}


/***************************************************************
**	Class   : MlEqHullAndWhite
**	Function: P
**	Returns : double
**	Comment : bond under Q_T
****************************************************************/

double MlEqHullAndWhite::P(long t, long T,double stateVar)
{
	int indexStart,indexEnd;

	double fwdDisc,B,var;
	return P(indexStart,indexEnd,fwdDisc,B,var,t,T,stateVar);
}

/***************************************************************
**	Class   : MlEqHullAndWhite
**	Function: P
**	Returns : double
**	Comment : this function also returns indexStart and indexEnd
****************************************************************/

double MlEqHullAndWhite::P(int& indexStart,int& indexEnd,double& fwdDisc,double& b,double& var,
						   long t, long T,double stateVar)
{
	fwdDisc = m_discountCurve->GetDiscountFactor(m_hDate->GetDate(),T)/
					 m_discountCurve->GetDiscountFactor(m_hDate->GetDate(),t);

	b		= 	B(indexStart,indexEnd,t,T);

	var		=	integral_g2(m_hDate->GetDate() ,t);

	double z	=	fwdDisc*exp(-0.5*b*b*var + b*stateVar);

	return z;
}


/***************************************************************
**	Class   : MlEqHullAndWhite
**	Function: P
**	Returns : double
**	Comment : this function uses indexStart and indexEnd
****************************************************************/

double MlEqHullAndWhite::P(long t, long T,double stateVar,int indexStart,int indexEnd,double fwdDisc,double B,double var)
{
//	double fwdDisc = m_discountCurve->GetDiscountFactor(m_hDate->GetDate(),T)/
//					 m_discountCurve->GetDiscountFactor(m_hDate->GetDate(),t);


	double dT		= m_hDate->GetYearFraction(T);
	double dt		= m_hDate->GetYearFraction(t);

//	double b    =   B(dt,dT,indexStart,indexEnd);
//	double var	=	integral_g2(0,dt,0,indexStart);//sos1integral_g2(t,T,indexStart,indexEnd);
	double z	=	fwdDisc*exp(-0.5*B*B*var + B*stateVar);

	return z;
}

/***************************************************************
**	Class   : MlEqHullAndWhite
**	Function: Z
**	Returns : double
**	Comment : 
****************************************************************/

double MlEqHullAndWhite::Z(long t,long T)
{
	double fwdDisc = m_discountCurve->GetDiscountFactor(m_hDate->GetDate(),T)/
					 m_discountCurve->GetDiscountFactor(m_hDate->GetDate(),t);

	return fwdDisc;
}

/***************************************************************
**	Class   : MlEqHullAndWhite
**	Function: CapletPrice
**	Returns : double
**	Comment : 
****************************************************************/

double MlEqHullAndWhite::CapletPrice(double B, long RateSettingDate,long RateMaturityDate,double strike,DayCountConventionEnum & dayCount )
{

	if ( RateSettingDate > RateMaturityDate || RateSettingDate < m_hDate->GetDate()){
		throw(" input error in caplet calculation");
	}

	MlEqDate dateToDouble(m_hDate->GetDate(),dayCount);

	double ZT		= Z(dateToDouble.GetDate(),RateMaturityDate);
	double Zt		= Z(dateToDouble.GetDate(),RateSettingDate);
	double delta	= dateToDouble.GetYearFraction(RateMaturityDate)-dateToDouble.GetYearFraction(RateSettingDate);
	double fwd		= ZT/Zt;
	double sig		= pow(B,2.0)*integral_g2(dateToDouble.GetDate(),RateSettingDate);
	sig				= sqrt(sig);
	strike			= 1.0/(1.0+delta*strike);
	double d1		= log(fwd/strike)/sig + sig/2.0;
	double d2		= d1-sig;

	double cplPrice = Zt*normal(-d2)-ZT/strike*normal(-d1);
	return cplPrice;
}


/***************************************************************
**	Class   : MlEqHullAndWhite
**	Function: CapletPrice
**	Returns : double
**	Comment : 
****************************************************************/

double MlEqHullAndWhite::CapletPrice(long RateSettingDate,long RateMaturityDate,double strike,DayCountConventionEnum & dayCount )
{
	double b = B(RateSettingDate,RateMaturityDate);
	return CapletPrice(b,RateSettingDate, RateMaturityDate, strike, dayCount );
}

/***************************************************************
**	Class   : MlEqHullAndWhite
**	Function: CapletPrice
**	Returns : double
**	Comment : 
****************************************************************/

double MlEqHullAndWhite::CapletPrice(int icalibDate,double strike)
{

//	calculates caplet with ratesetting at m_calibrationDates[icalibdate] and ratematurity m_calibrationDates[icalibdate+1]

	if ( icalibDate > m_calibrationDates.getsize()-2 ){
		throw("incorrect icalibdate entered in calibration caplet calculation");
	}

	int idate		= icalibDate;
	double ZT		= Z(m_hDate->GetDate(),m_calibrationDates[idate+1]);
	double Zt		= Z(m_hDate->GetDate(),m_calibrationDates[idate]);
	double T		= m_hDate->GetYearFraction(m_calibrationDates[idate+1]);
	double t		= m_hDate->GetYearFraction(m_calibrationDates[idate]);
	double delta	= T-t;
	double fwd		= ZT/Zt;
	double sig		= pow(B(t,T,idate,idate+1),2.0)*integral_g2(0,T,0,idate+1);//sos1integral_g2(t,T,idate,idate+1);
	sig				= sqrt(sig);
	strike			= 1.0/(1.0+delta*strike);
	double d1		= log(fwd/strike)/sig + sig/2.0;
	double d2		= d1-sig;

	double cplPrice = Zt*normal(-d2)-ZT/strike*normal(-d1);
	return cplPrice;
}

/***************************************************************
**	Class   : MlEqHullAndWhite
**	Function: SwaptionPrice
**	Returns : double
**	Comment : 
****************************************************************/

double MlEqHullAndWhite::SwaptionPrice(GVector<long>& SwapFixingDates,double swaptionStrike,DayCountConventionEnum & dayCount,double integ_g2 ,int callPut)
{

//  SwapFixingDates[0] corresponds to the Swaption Maturity

//  callPut = -1 : receiver swaption
//  callPut =  1 : payer swaption

	double cp;
	if ( callPut == 1 ){
		cp = -1.0;// notice formula below describes receiver swaption
	}else if ( callPut == -1 ){
		cp = 1.0;
	}else{
		throw("incorrect callPut input chosen");
	}

	MlEqDate dateToDouble(m_hDate->GetDate(),dayCount);

	double Z0 = Z(m_hDate->GetDate(),SwapFixingDates[0]);
	double Zi,d1,d2,p,delta,strike;
	double swaptionVal=0.0;

//  check wether today's day is swaption maturity	

	if ( m_hDate->GetDate() == SwapFixingDates[0] )
	{
		for ( int i = 1 ; i < SwapFixingDates.getsize(); i++ )
		{
			Zi		= Z(m_hDate->GetDate(),SwapFixingDates[i]);
			delta	= dateToDouble.GetYearFraction(SwapFixingDates[i])-dateToDouble.GetYearFraction(SwapFixingDates[i-1]);
			strike = swaptionStrike*delta;

			if ( i == SwapFixingDates.getsize()-1){
				strike += 1.0;
			}
			swaptionVal += strike*Zi;
		}

		swaptionVal = MlEqMaths::Max((double)cp*(swaptionVal-1.0),0.0);
		return swaptionVal;
	}


	swaptionRootDataContainer c;
	c.pHullWhite		= this;
	c.swaptionStrike	= swaptionStrike;
	c.SwapFixingDates	= SwapFixingDates;
	c.dayCount			= dayCount;

	int retVal;
	double stateVar		= -1e99;
	int rootfind_flag	= eROOTFIND_BRENT_GROWRANGE;
	double sig			= sqrt(integ_g2);
	double lower_x		= -4.0*sig;
	double upper_x		= 4.0*sig;
	double accuracy		= 1e-6;
	int max_tries		= 50;

	retVal =  rootfind_solve(rootfind_flag,swaptionRoot, 
							  lower_x,upper_x,
							  accuracy, max_tries, 
							  0.0,&c, &stateVar,NULL );


	for ( int i = 1 ; i < SwapFixingDates.getsize(); i++ )
	{
		Zi		= Z(m_hDate->GetDate(),SwapFixingDates[i]);
		p		= P(SwapFixingDates[0], SwapFixingDates[i],stateVar);
		sig		= pow(B(SwapFixingDates[0], SwapFixingDates[i]),2.0)*integ_g2;
		sig		= sqrt(sig);
		d1		= log(Zi/Z0/p)/sig + sig/2.0;
		d2		= d1-sig;
		delta	= dateToDouble.GetYearFraction(SwapFixingDates[i])-dateToDouble.GetYearFraction(SwapFixingDates[i-1]);
		strike = swaptionStrike*delta;

		if ( i == SwapFixingDates.getsize()-1){
			strike += 1.0;
		}

		swaptionVal +=	strike*cp*(normal(cp*d1)*Zi-p*normal(cp*d2)*Z0);
	}

	return swaptionVal;
}

/***************************************************************
**	Class   : MlEqHullAndWhite
**	Function: SwaptionPrice
**	Returns : double
**	Comment : 
****************************************************************/

double MlEqHullAndWhite::SwaptionPrice(GVector<long>& SwapFixingDates,double swaptionStrike,DayCountConventionEnum & dayCount,int callPut )
{
	double integ_g2 = integral_g2(m_hDate->GetDate(),SwapFixingDates[0]);

	double test = BlackScholesSwaptionPrice(SwapFixingDates,swaptionStrike,0.2,dayCount ,callPut);
	return SwaptionPrice(SwapFixingDates,swaptionStrike,dayCount,integ_g2,callPut );
}


/***************************************************************
**	Class   : MlEqHullAndWhite
**	Function: BlackScholesSwaptionPrice
**	Returns : double
**	Comment : 
****************************************************************/

double MlEqHullAndWhite::BlackScholesSwaptionPrice(GVector<long>& SwapFixingDates,double swaptionStrike,double vol,DayCountConventionEnum & dayCount ,int callPut)
{
	return  BlackScholesSwaptionPricer(m_hDate->GetDate(),SwapFixingDates,swaptionStrike,vol,dayCount ,*m_discountCurve, callPut);
}

/***************************************************************
**	Class   : MlEqHullAndWhite
**	Function: BlackScholesCapletPrice
**	Returns : double
**	Comment : 
****************************************************************/

double MlEqHullAndWhite::BlackScholesCapletPrice(long RateSettingDate,long RateMaturityDate,double strike,double vol,DayCountConventionEnum & dayCount,int callput=1)
{
	return BlackScholesCapletPricer(m_hDate->GetDate(), RateSettingDate,RateMaturityDate,strike,vol,dayCount,*m_discountCurve,callput);
}

/***************************************************************
**	Class   : MlEqHullAndWhite
**	Function: CapletPrice
**	Returns : double
**	Comment : 
****************************************************************/
/*
double MlEqHullAndWhite::SwaptionPrice(GVector<long>& SwapFixingDates,double swaptionStrike,DayCountConventionEnum & dayCount,int callPut )
{

//  SwapFixingDates[0] corresponds to the Swaption Maturity

	double cp;
	if ( callPut == 1 ){
		cp = 1.0;
	}else if ( callPut == -1 ){
		cp = -1.0;
	}else{
		throw("incorrect callPut input chosen");
	}

	MlEqDate dateToDouble(m_hDate->GetRefDate(),dayCount);

	double Z0 = Z(m_hDate->GetDate(),SwapFixingDates[0]);
	double Zi,d1,d2,p,delta,strike;
	double swaptionVal=0.0;

// check wether today's day is swaption maturity	

	if ( m_hDate->GetDate() == SwapFixingDates[0] )
	{
		for ( int i = 1 ; i < SwapFixingDates.getsize(); i++ )
		{
			Zi		= Z(m_hDate->GetDate(),SwapFixingDates[i]);
			delta	= m_hDate->DateToDouble(SwapFixingDates[i])-m_hDate->DateToDouble(SwapFixingDates[i-1]);
			strike = swaptionStrike*delta;

			if ( i == SwapFixingDates.getsize()-1){
				strike += 1.0;
			}
			swaptionVal += strike*Zi;
		}

		swaptionVal = MlEqMaths::Max((double)cp*(swaptionVal-1.0),0.0);
		return swaptionVal;
	}


	swaptionRootDataContainer c;
	c.pHullWhite		= this;
	c.swaptionStrike	= swaptionStrike;
	c.SwapFixingDates	= SwapFixingDates;
	c.dayCount			= dayCount;

	int retVal;
	double stateVar = -1e99;
	
	int rootfind_flag	= eROOTFIND_BRENT_GROWRANGE;
	double sig			= sqrt(integral_g2(m_hDate->GetDate(),SwapFixingDates[0]));
	double lower_x		= -4.0*sig;
	double upper_x		= 4.0*sig;
	double accuracy		= 1e-6;
	int max_tries		= 50;


	retVal =  rootfind_solve(rootfind_flag,swaptionRoot, 
							  lower_x,upper_x,
							  accuracy, max_tries, 
							  0.0,&c, &stateVar,NULL );

	


	for ( int i = 1 ; i < SwapFixingDates.getsize(); i++ )
	{
		Zi		= Z(m_hDate->GetDate(),SwapFixingDates[i]);
		p		= P(SwapFixingDates[0], SwapFixingDates[i],stateVar);

		sig		= pow(B(SwapFixingDates[0], SwapFixingDates[i]),2.0)*integral_g2(m_hDate->GetDate(), SwapFixingDates[0]);
		sig		= sqrt(sig);
		d1		= log(Zi/Z0/p)/sig + sig/2.0;
		d2		= d1-sig;
		delta	= m_hDate->DateToDouble(SwapFixingDates[i])-m_hDate->DateToDouble(SwapFixingDates[i-1]);
		strike = swaptionStrike*delta;

		if ( i == SwapFixingDates.getsize()-1){
			strike += 1.0;
		}

		swaptionVal +=	strike*cp*(normal(cp*d1)*Zi-p*normal(cp*d2)*Z0);

	}
	return swaptionVal;
}	
*/

/***************************************************************
**	Class   : MlEqHullAndWhite
**	Function: initialize
**	Returns : nothing
**	Comment : 
****************************************************************/


void MlEqHullAndWhite::initialize(GVector<long>& Dates,CVector& beta,CVector& gamma,double lambda,MlEqZeroCurveHandle	discountCurve,	CVector& capletVols,CVector& capletStrikes,CVector& swaptionVols,CVector& swaptionStrikes,CMatrix& swaptionDates,DayCountConventionEnum & dayCount)
{

	m_calibrationDates	= Dates;
	m_beta_ti			= beta;
	m_gamma_ti			= gamma;
	m_lambda			= lambda;
	m_discountCurve		= discountCurve;	

	m_times.resize(m_calibrationDates.getsize());

	for ( int i = 0 ; i < m_calibrationDates.getsize(); i++ ){
		m_times[i] = m_hDate->GetYearFraction(m_calibrationDates[i]); 
	}

	// test calibration

//	CVector capletVols(m_calibrationDates.getsize()-2);

//	calibrateToCaplets(capletVols,capletStrikes);
	calibrateToCaplets(capletVols,capletStrikes,dayCount);

	GVector < GVector<long> > calibSwaptionDates;
	calibSwaptionDates.resize(swaptionVols.getsize());

	for (int i = 0 ; i < swaptionDates.rows(); i++ )
	{
		int k = 0;
		for ( int j = 0 ; j < swaptionDates.cols(); j++ )
		{
			if(swaptionDates[i][j] > 0 ){
				k++;
			}
			else{
				break;
			}


//			swaptionDates[i][j] > 0 : k++ ? break;
		}

		calibSwaptionDates[i].resize(k);
		for ( int j = 0 ; j < k; j++ ){
			calibSwaptionDates[i][j] = swaptionDates[i][j];
		}

	}


	calibrateToSwaptions(swaptionVols,swaptionStrikes,calibSwaptionDates,capletVols,capletStrikes,dayCount);
	SwaptionPrice(calibSwaptionDates[0],swaptionStrikes[0],dayCount,1);

}


/***************************************************************
**	Class   : none
**	Function: swaptionRoot
**	Returns : bool
**	Comment : 
****************************************************************/
	
bool swaptionRoot(double x,void* vp,double* f)
{		
	swaptionRootDataContainer* p = (swaptionRootDataContainer*)vp;
	
	MlEqHullAndWhiteHandle pHullWhite = p->pHullWhite;
	double swaptionStrike	= p->swaptionStrike;
	double strike			= 1.0;
	GVector<long>&	fixings	= p->SwapFixingDates;
	
	MlEqDate dateToDouble(pHullWhite->m_hDate->GetDate(),p->dayCount);
	
	double delta;
	double res = -strike;
	int i;
	for (  i = 1 ; i < fixings.getsize()-1; i++ )
	{	
		delta	=  dateToDouble.GetYearFraction(fixings[i]) - dateToDouble.GetYearFraction(fixings[i-1]);
		res		+= swaptionStrike*delta*pHullWhite->P(fixings[0],fixings[i],x);
	}	
	
//	int i	=	fixings.getsize();
	delta	=	dateToDouble.GetYearFraction(fixings[i]) - dateToDouble.GetYearFraction(fixings[i-1]);
	res		+=	(1.0+swaptionStrike*delta)*pHullWhite->P(fixings[0],fixings[i],x);
	*f		=	res;
	
	return true;
}	
	
	
	
/***************************************************************
**	Class   : MlEqHullAndWhiteCalibration 
**	Function: initialize
**	Returns : void
**	Comment : 
****************************************************************/
	
void MlEqHullAndWhite::initialize(CVector& calibrationDates, long rateSetting,long rateMaturity)
{	
	m_calibrationDates.resize(calibrationDates.getsize());
	for ( int i = 0 ; i < calibrationDates.getsize();i++){
		m_calibrationDates[i] = calibrationDates[i];
	}
		
	m_times.resize(m_calibrationDates.getsize());
	m_gamma_ti.resize(m_calibrationDates.getsize());
	for (int i = 0 ; i < m_calibrationDates.getsize(); i++ ){
		m_times[i]		=	m_hDate->GetYearFraction(m_calibrationDates[i]);
		m_gamma_ti[i]	=	1.0;
	}
	
	m_lambda = 1.0/25.0;
	double res = B(rateSetting,rateMaturity);
}	
	

/***************************************************************
**	Class   : none
**	Function: capletRoot
**	Returns : bool
**	Comment : 
****************************************************************/
	
bool capletRoot(double x,void* vp,double* f)
{	
	double res;

	swaptionRootDataContainer* p = (swaptionRootDataContainer*)vp;
	
	MlEqHullAndWhiteHandle pHullWhite	= p->pHullWhite;


	long RateSettingDate	= p->RateSettingDate;	
	long RateMaturityDate	= p->RateMaturityDate;
	double capletStrike		= p->capletStrike;
	double price			= p->price;	

	res		=	pHullWhite->CapletPrice(x,RateSettingDate,RateMaturityDate,capletStrike,p->dayCount );
	res		-=	price;
	*f		=	res;
	
	return true;
}	

	
/***************************************************************
**	Class   : MlEqHullAndWhite 
**	Function: B
**	Returns : double
**	Comment : this function usses indexStart and indexEnd
****************************************************************/
	
void MlEqHullAndWhite::calibrateToCaplets(CVector& capletVols,CVector& capletStrikes,DayCountConventionEnum & dayCount)
{	
	
	if ( ( capletVols.getsize() != m_calibrationDates.getsize()-2 ) ||
		 ( m_calibrationDates[0] != m_hDate->GetDate() ) )
	{
		throw("input error in setting up HullWhite: check whether first calibration date coincides with refDate and that number of caplet vols provided is correct");
	}
	
	swaptionRootDataContainer c;

	c.pHullWhite		= this;
	c.dayCount			= dayCount;
	
	int retVal;
	double B			= -1e99;	
	int rootfind_flag	= eROOTFIND_BRENT_NOGROW;
	double lower_x		= 1e-6;
	double upper_x;
	double accuracy		= 1e-6;
	int max_tries		= 50;
	int callput			= 1;

	for ( int i = 0 ; i < capletVols.getsize(); i++ )
	{

		upper_x				=	sqrt(m_hDate->GetYearFraction(m_calibrationDates[i+1]));

		c.RateSettingDate	=	m_calibrationDates[i+1];
		c.RateMaturityDate	=	m_calibrationDates[i+2];
		c.capletStrike		=	capletStrikes[i],	
		c.price				=	BlackScholesCapletPrice(m_calibrationDates[i+1],m_calibrationDates[i+2],capletStrikes[i],capletVols[i],dayCount,callput);

		retVal =  rootfind_solve(rootfind_flag,capletRoot, 
								  lower_x,upper_x,
								  accuracy, max_tries, 
								  0.0,&c, &B,NULL );


		m_gamma_ti[i+1]	=	B*m_lambda/(exp(-m_times[i+1]*m_lambda)-exp(-m_times[i+2]*m_lambda)); 
	}

}		
		
		
/***************************************************************
**	Class   : none
**	Function: swaptionCalibRoot
**	Returns : bool
**	Comment : 
****************************************************************/
	
int	 swaptionCalibRoot (const CVector& x, void* vp, int calc_for_grad, CVector& f)
{	
	swaptionRootDataContainer* p = (swaptionRootDataContainer*)vp;
	
	MlEqHullAndWhiteHandle pHullWhite = p->pHullWhite;
	GVector <GVector<long> >&	fixings	= p->AllSwapFixingDates;
	int callPut = 1;
	
	for ( int i = 0 ; i < x.getsize(); i++ ){
		pHullWhite->m_beta_ti[i] = exp(x[i]);
	}

	pHullWhite->calibrateToCaplets(p->capletVols,p->capletStrikes,p->dayCount);

	for ( int i = 0 ; i < fixings.getsize(); i++ ){
		f[i] =  pHullWhite->SwaptionPrice(fixings[i],p->swaptionStrikes[i],p->dayCount,callPut);
	}
	for ( int i = 0 ; i < fixings.getsize(); i++ ){
		f[i] -= p->swaptionVals[i];
	}
	

	return 0;
}	
	
extern void rootfind_underdetermined_solve(const CVector& initial_x, const CVector& bump_x, const CVector& tolerances, 
								   p_dev_underdet_func fn, void* vp, int max_tries, 
								   int max_restarts, const CMatrix& weights, CVector& found_x);

/***************************************************************
**	Class   : MlEqHullAndWhite 
**	Function: B
**	Returns : double
**	Comment : this function usses indexStart and indexEnd
****************************************************************/

void MlEqHullAndWhite::calibrateToSwaptions(CVector& swaptionVols,CVector& swaptionStrikes,GVector < GVector< long > > & swaptionDates,CVector& capletVols,CVector& capletStrikes,DayCountConventionEnum & dayCount)
{
	if ( m_gamma_ti.getsize() <= swaptionVols.getsize() ){
		throw("too many swaptions entered");
	}

	if ( ( swaptionVols.getsize() != swaptionStrikes.getsize() ) ||
		 ( swaptionVols.getsize() != swaptionDates.getsize() )){
		throw("input error: check arraysize of swaption vols and swaption strikes as well as number of swaptions entered");
	}

	int nvariables = m_beta_ti.getsize();

//	for ( int i = 0; i < m_beta_ti.getsize(); i++ ){
//		m_beta_ti[i] = 1.0;
//	}

	swaptionRootDataContainer c;
	c.pHullWhite				= this;
	c.dayCount					= dayCount;
	c.AllSwapFixingDates		= swaptionDates;
	c.swaptionStrikes			= swaptionStrikes;
	c.capletVols				= capletVols;
	c.capletStrikes				= capletStrikes;

	CVector swaptionVals;
	int callPut = 1;
	swaptionVals.resize(swaptionVols.getsize());
	for (int i = 0 ; i < swaptionVals.getsize(); i++ ){
		swaptionVals[i]	= BlackScholesSwaptionPrice(swaptionDates[i],swaptionStrikes[i],swaptionVols[i],dayCount,callPut);
	}
	c.swaptionVals = swaptionVals;

	int max_tries = 500;
	int max_restarts = 50;
	CVector tolerance(swaptionVols.getsize());

	for ( int i = 0; i < tolerance.getsize(); i++ ){
		tolerance[i] = 1e-6;
	}

	CVector bump(nvariables);
	for ( int i = 0; i < nvariables; i++ ){
		bump[i] = 0.008;
	}

	
	CMatrix weights;
	CVector found_x;	

	CVector initial_x(m_beta_ti.getsize());

	for ( int i = 0; i < m_beta_ti.getsize(); i++ ){
		initial_x[i] = log(m_beta_ti[i]);
	}

	rootfind_underdetermined_solve(initial_x, bump, tolerance, 
								   swaptionCalibRoot, &c, max_tries, 
								   max_restarts, weights,found_x);


}


/***************************************************************
**	Class  : CStochIRForwardSkewMC 
**	Routine: getDiscountFactor
**	Returns: double
**	Action : 
****************************************************************/


double CStochIRForwardSkewMC::getDiscountFactor(int idateAsOf,int ipath,long To)
{
	int indexStart,indexEnd;
	double fwdDisc,B, var;
	return getDiscountFactor_5(indexStart,indexEnd,fwdDisc,B, var,idateAsOf,ipath,To);
}	


/***************************************************************
**	Class  : CStochIRForwardSkewMC 
**	Routine: getDiscountFactor
**	Returns: double
**	Action : 
****************************************************************/


double CStochIRForwardSkewMC::getDiscountFactor_5(int& indexStart,int& indexEnd,double& fwdDisc,double& B,double& var,
												  int idateAsOf,int ipath,long To)
{	
	int idate				=	idateAsOf;
	if ( idate == 0 ){

		pHullWhite->Locate(indexStart,indexEnd,m_hDate->GetDate(),To);

		return pHullWhite->Z(m_hDate->GetDate(),To);
	}

	const std::vector<long>& mcDates	=	m_mMCHelper.GetMCDates() ;
	double stateVar						=	getStateVar(ipath,idate);
	double disc	= pHullWhite->P(indexStart,indexEnd,fwdDisc,B,var,mcDates[idate],To,stateVar);
	return disc;
}	

/***************************************************************
**	Class  : CStochIRForwardSkewMC 
**	Routine: getDiscountFactor
**	Returns: double
**	Action : 
****************************************************************/


double CStochIRForwardSkewMC::getDiscountFactor_5(int idateAsOf,int ipath,int indexStart,int indexEnd,double fwdDisc,double B,double var)
{	
	if ( idateAsOf > indexStart ){
		throw("indexing problem in getDiscountFactor_5");
	}

	int idate				=	idateAsOf;
	const std::vector<long>& mcDates	=	m_mMCHelper.GetMCDates() ;

	if ( idate == 0 ){
		
		return pHullWhite->Z(m_hDate->GetDate(),mcDates[indexEnd]);
	}

	double stateVar	 =	getStateVar(ipath,idate);

	double disc = pHullWhite->P(mcDates[indexStart], mcDates[indexEnd],stateVar,idateAsOf,indexEnd,fwdDisc,B,var);

	return disc;
}	

/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: GetDiscounts
**	Returns : CVector
**	Comment : calculates stochastic discounts along a monteCarlo Path 
****************************************************************/

const CVector& CStochIRForwardSkewMC::GetDiscounts(int ipath)
{
	
//  calculates stochastic discounts along a MonteCarlo Path
	
	int idate;
	const std::vector<long>& mcDates = m_mMCHelper.GetMCDates() ;
	
	if ( m_stochasticDiscountsAlongPath.getsize() == 0 ){
		m_stochasticDiscountsAlongPath.resize(m_nDates+1);
	}
	
	if ( ipath == 0 && m_saveStochasticDiscounts ){
		m_stochasticDiscounts.resize(m_nPaths,m_nDates);
	}

	m_stochasticDiscountsAlongPath[0] = 1.0;
	for ( idate = 1 ; idate < m_nDates; idate++ )
	{
		m_stochasticDiscountsAlongPath[idate] = 
				m_stochasticDiscountsAlongPath[idate-1] * 
				getDiscountFactor(idate-1,ipath,mcDates[idate]);
	}

	
	return m_stochasticDiscountsAlongPath;
}	

CStochIRForwardSkewMC::~CStochIRForwardSkewMC()
{
	m_stochasticDiscounts.resize(0,0);
}


/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: GetDiscounts
**	Returns : CVector
**	Comment : calculates stochastic discounts along a monteCarlo Path 
****************************************************************/


double	CStochIRForwardSkewMC::GetDiscountToDate(int ipath,long To,int idateAsOf)
{

	double disc = getDiscountFactor(idateAsOf,ipath,To);
	return disc;
}


/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: GetDiscounts
**	Returns : CVector
**	Comment : calculates stochastic discounts along a monteCarlo Path 
****************************************************************/

double	CStochIRForwardSkewMC::GetStochasticDiscount(int ipath,int idate)
{
	
// this is the stochastic discount to maturity mcDates[idate]

	if ( m_saveStochasticDiscounts ){

		double res = m_stochasticDiscounts[ipath][idate];
		return res;
	}


	const std::vector<long>& mcDates = m_mMCHelper.GetMCDates() ;	
	double res = 1.0;

	for ( int i = 0 ; i < idate; i++ )
	{
		res *= getDiscountFactor(i,ipath,mcDates[i+1]);

		if ( i > 0 ){
			res *= m_convexityAdj[i-1];
		}
	}

	return res;
}

/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: GetDiscounts
**	Returns : CVector
**	Comment : calculates stochastic discounts along a monteCarlo Path 
****************************************************************/

double	CStochIRForwardSkewMC::GetDiscount(int ipath,int idate,int idateAsOf)
{

	if ( idate < idateAsOf ){
		throw("indexing error in CStochIRForwardSkewMC::GetDiscount");
	}

	const std::vector<long>& mcDates = m_mMCHelper.GetMCDates() ;

	double fwdDisc=0,B=0,var=0;
	int indexStart=0,indexEnd=0;

	double res = getDiscountFactor_5(indexStart,indexEnd,fwdDisc,B,var,
						idateAsOf,ipath,mcDates[idate]);

	return res;

}
/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: GetDiscounts
**	Returns : CVector
**	Comment : calculates stochastic discounts along a monteCarlo Path 
****************************************************************/

double	CStochIRForwardSkewMC::GetDiscount(int ipath,int idate,int jdate,int idateAsOf)
{
	if ( idate < idateAsOf ){
		throw("indexing error in --Discount");
	}

	const std::vector<long>& mcDates = m_mMCHelper.GetMCDates() ;

	double disc		=  GetDiscountToDate(ipath,mcDates[jdate],idateAsOf);
//	disc			/= GetDiscountToDate(ipath,mcDates[idate],idateAsOf);

	return disc;
}

	
/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: calculateOneStepDiscounts
**	Returns : nothing
**	Comment : calculates discounts factors from idate-> idate+1
****************************************************************/
	
void CStochIRForwardSkewMC::initIRStateVars(int idate)
{	
	if ( idate == 0 )
	{
		for (int idate = 0 ; idate < m_nDates; idate++ ){
			createStateVariables(idate);
		}
	}


}	
	
/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: calculateOneStepDiscounts
**	Returns : nothing
**	Comment : calculates discounts factors from idate-> idate+1
****************************************************************/
	
void CStochIRForwardSkewMC::calculateOneStepDiscounts(int idate)
{	
	
//  calculates discounts factors from idate-> idate+1
	
	int idateAsOf = idate;
	const std::vector<long>& mcDates = m_mMCHelper.GetMCDates() ;

	if ( m_oneStepDiscounts.getsize() == 0 ){
		m_oneStepDiscounts.resize(m_nPaths);
	}
	
	double fwdDisc,B,var;
	int indexStart=0,indexEnd=0;
	int ipath = 0;

	m_oneStepDiscounts[0] = getDiscountFactor_5(indexStart,indexEnd,fwdDisc,B,var,idateAsOf,ipath,mcDates[idate+1]);
	for (ipath = 1 ; ipath < m_nPaths; ipath++ )
	{
		m_oneStepDiscounts[ipath] = getDiscountFactor_5(idateAsOf,ipath,indexStart,indexEnd,fwdDisc,B,var);
	}


}


	
/***************************************************************
**	Class  : CStochIRForwardSkewMC 
**	Routine: getStateVar
**	Returns: double
**	Action : 
****************************************************************/

double CStochIRForwardSkewMC::getStateVar(int ipath,int idate)
{
	if ( idate == 0 ){
		return 0.0;
	}
	return m_StateVar[ipath][idate-1];//xsos
}


/***************************************************************
**	Class  : CStochIRForwardSkewMC 
**	Routine: createStateVariables
**	Returns: double
**	Action : 
****************************************************************/


void CStochIRForwardSkewMC::createStateVariables()
{
	for ( int idate = 0 ; idate < m_nDates; idate++ ){ 
		createStateVariables(idate);
	}
}
	
/***************************************************************
**	Class  : CStochIRForwardSkewMC 
**	Routine: createStateVariables
**	Returns: void
**	Action : 
****************************************************************/


void CStochIRForwardSkewMC::createStateVariables(int idate)
{

	// creates statevariable S_{idate+1}
	// m_StateVar[ipath][idate] stores state Variable S_{idate+1}

	if ( m_StateVar.rows() == 0 ){
		m_StateVar.resize(m_nPaths,m_nDates);
	}
	
	if ( m_convexityAdj.getsize() == 0 ){
		m_convexityAdj.resize(m_nDates);
		m_bareConvexityAdj.resize(m_nDates);
//		m_bareConvexityAdj.resize(m_nDates);
	}

	const std::vector<long>& mcDates = m_mMCHelper.GetMCDates() ;
	
//  m_numberOfFactors*idate   : random for stochastic vol
//  m_numberOfFactors*idate+1 : random for spot
//  m_numberOfFactors*idate+2 : random for stochastic interest rates
	
	int n_offset = m_numberOfFactors*idate+1;// random number for spot
	int offset	 = n_offset+1;// random numbers for stochastic interest rates
	
	double rand,stateVar,convexityAdj;
	
	int nn_offset = m_numberOfFactors*(idate+1)+1;// random number for spot

	convexityAdj = -pHullWhite->B(mcDates[idate],mcDates[idate+1])*
					pHullWhite->integral_g2(mcDates[0],mcDates[idate]);
		
	m_convexityAdj[idate] = convexityAdj*pHullWhite->B(mcDates[idate],mcDates[idate+1]);
	m_bareConvexityAdj[idate] = convexityAdj;

	if ( idate > 0 ){
		m_convexityAdj[idate] += m_convexityAdj[idate-1];
	}

	if ( idate == m_nDates-1 ){
		for ( int i = 0; i < m_nDates; i++ ){
			m_convexityAdj[i] = exp(m_convexityAdj[i]);
		}
	}

	for ( int ipath = 0 ; ipath < m_nPaths; ipath++ )
	{	
//		rand = m_corr*GetNextRandom(ipath,n_offset)+m_xcorr*GetNextRandom(ipath,offset);	
//		rand = m_corr*GetNextRandom(ipath,nn_offset)+m_xcorr*GetNextRandom(ipath,offset);//xsossos	

		rand = GetNextRandom(ipath,offset);	

		stateVar = m_stateVarVols[idate]*rand + 0.0*convexityAdj;//xsos
		m_StateVar[ipath][idate] = stateVar;
		if ( idate > 0 ){
			m_StateVar[ipath][idate] +=	m_StateVar[ipath][idate-1];
		}
	}	



}

/***************************************************************
**	Class  : CStochIRForwardSkewMC 
**	Routine: calculateConvextityAdjustment
**	Returns: double
**	Action : 
****************************************************************/

double	CStochIRForwardSkewMC::calculateConvextityAdjustment(int idate,long T)
{

	const std::vector<long>& mcDates = m_mMCHelper.GetMCDates() ;

	double convexAdj=0.0;
	for (int i = 0 ; i <= idate; i++ ){
		convexAdj += m_bareConvexityAdj[i]*pHullWhite->B(mcDates[i],T);
	}

	convexAdj = exp(convexAdj);
	return convexAdj;
}

/***************************************************************
**	Class  : CStochIRForwardSkewMC 
**	Routine: calculateConvextityAdjustment
**	Returns: double
**	Action : 
****************************************************************/

double	CStochIRForwardSkewMC::calculateConvextityAdjustment(int iDate,int jDate)
{

	if ( jDate == iDate +1 ){
		return 1.0;// no convexity adjustments needed in this case");
	}

	if ( iDate > jDate ){
		throw("startDate greater than final date");
	}
	if ( jDate > m_nDates ){
		throw("indexing problem in calculateConvextityAdjustment");
	}

	double res = 1.0;
	for ( int i = iDate ; i < iDate; i++ )
	{
		if ( i > 0 ){
			res *= m_convexityAdj[i-1];
		}

	}
	return res;
}


/***************************************************************
**	Class  : CStochIRForwardSkewMC 
**	Routine: getlogFwdFactors
**	Returns: double
**	Action : 
****************************************************************/

double CStochIRForwardSkewMC::getlogFwdFactors(int idate,int ipath)
{

// note that CForwardSkewMC::getlogFwdFactors(idate,ipath)  has already forward disc included

	double logFwdFactor = 0.0;

	if ( !m_minimalHybrid )
	{

		const std::vector<long>& mcDates = m_mMCHelper.GetMCDates() ;
		double fwd	= pHullWhite->Z(mcDates[idate],mcDates[idate+1]);
			   fwd	/=	m_oneStepDiscounts[ipath];

		if ( idate > 0 ){
//			fwd /= m_convexityAdj[idate-1];//sos

		}

		logFwdFactor = log(fwd);
	}

	logFwdFactor += CForwardSkewMC::getlogFwdFactors(idate,ipath) ;


	if ( m_saveStochasticDiscounts )
	{
		const std::vector<long>& mcDates = m_mMCHelper.GetMCDates() ;

		double res = m_oneStepDiscounts[ipath];//sos getDiscountFactor(idate,ipath,mcDates[idate+1]);

		if ( idate ){
				res *= m_convexityAdj[idate-1];
		}
		else{

			if ( ipath == 0 && idate == 0 ){
					m_stochasticDiscounts.resize(m_nPaths,mcDates.size());
			}

			m_stochasticDiscounts[ipath][idate] = 1.0;
		}

		m_stochasticDiscounts[ipath][idate+1] = m_stochasticDiscounts[ipath][idate]*res;
	}
	

	return logFwdFactor;
}


 
/***************************************************************
**	Class  : CStochIRForwardSkewMC 
**	Routine: setUpNextTimeStep
**	Returns: void
**	Action : 
****************************************************************/

void CStochIRForwardSkewMC::setUpNextTimeStep(int idateCurrent)
{
	if ( m_StateVar.rows() == 0 ){
		createStateVariables();
	}
	calculateOneStepDiscounts(idateCurrent);
}



/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: initializeNumberOfFactors
**	Returns : nothing
**	Comment : initializes forward skew MonteCarlo
****************************************************************/

void CStochIRForwardSkewMC::initializeNumberOfFactors()
{
	m_numberOfFactors = 3;

}

/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: initIR
**	Returns : nothing
**	Comment : 
****************************************************************/

void  CStochIRForwardSkewMC::initIR(double fCorrelation,bool isMinimalHybrid)
{
	if ( isMinimalHybrid ){
		m_centerFwds	= false;
	}else{
		m_centerFwds	= false;
	}

	m_corr			= fCorrelation;
	m_xcorr			= sqrt(1.0-m_corr*m_corr);
	m_minimalHybrid = isMinimalHybrid;

	if ( m_stateVarVols.getsize() != 0 ){
		return;
	}

	m_stateVarVols.resize(m_nDates);
	const std::vector<long>& mcDates	= m_mMCHelper.GetMCDates();;

	for (int i = 0 ; i < m_nDates; i++ )
	{
		m_stateVarVols[i]	= pHullWhite->integral_g2(mcDates[i],mcDates[i+1]);
		m_stateVarVols[i]	= sqrt(m_stateVarVols[i]);
	}
}



/***************************************************************
**	Class  : CStochIRForwardSkewMC 
**	Routine: generatePath
**	Returns: nothing
**	Action : executes one Monte Carlo timestep
****************************************************************/

void  CStochIRForwardSkewMC::calculateOptions(CVector& results,MlEqAssetHandle asset, int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor)
{	
	if ( m_assets.size() == 0 && m_nAssets == 1 ){
		// asset have not been set ( e.g. forward skews must have been entered explicitely
		int iasset = 0;
		calculateOptions(results,startIndex,endIndex,Strikes,controlVariate,impliedVolFlag,factor,iasset);
		return;
	}

	int iasset = findAsset(asset);
	if ( iasset != 0 ){
		throw("only one asset should be set in single factor monte carlo");
	}

	calculateOptions(results,startIndex,endIndex,Strikes,controlVariate,impliedVolFlag,factor,iasset);
}


/***************************************************************
**	Class  : CStochIRForwardSkewMC 
**	Routine: generatePath
**	Returns: nothing
**	Action : executes one Monte Carlo timestep
****************************************************************/

void CStochIRForwardSkewMC::calculateOptions(CVector& results,int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor,int iasset)
{	
	
	if ( iasset > 0 ){
		throw("iasset for a single factor monte carlo need to be set to zero");
	}

	int ipath,istrike;
	results.resize(Strikes.size());
	CVector values(Strikes.size());
	CVector ControlVarValues(Strikes.size());
	
	vector< MlEqStrike >  strikes; 
	strikes.resize(Strikes.size());
	
	for ( istrike = 0 ; istrike < strikes.size(); istrike++ )
	{
		MlEqStrike::convertStrikes(strikes[istrike],*Strikes[istrike]);
	}



	double payoff,disc;
	for ( ipath = 0 ; ipath < m_nPaths; ipath++ )
	{
		for ( istrike = 0 ; istrike < strikes.size(); istrike++ )
		{
			payoff = MlEqMaths::Max(factor*GetPathValue(ipath,endIndex)/GetPathValue(ipath,startIndex)-strikes[istrike].m_strike,0.0);
			completePayoff(payoff,ipath);


			disc = GetStochasticDiscount(ipath,endIndex)/GetStochasticDiscount(ipath,startIndex);

			payoff *= disc;
			values[istrike] += payoff;
		}
	}

	disc = pHullWhite->Z(GetMCDate(startIndex),GetMCDate(endIndex));


	averageMonteCarloResults(values);
	for ( int i = 0; i < values.getsize(); i++ ){
		values[i] /= factor*disc;
	}


//	call implied volcalculator here

	if ( impliedVolFlag )
	{
		for ( istrike = 0 ; istrike < results.getsize(); istrike++ )
		{
			results[istrike] = GetImpliedVol(startIndex,endIndex,values[istrike],strikes[istrike].m_strike/factor);
		}
	}
	else
	{
		for ( istrike = 0 ; istrike < results.getsize(); istrike++ )
		{
			results[istrike] = values[istrike];
		}
	}
}	
	


/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: Initialize
**	Returns : nothing
**	Comment : 
****************************************************************/

void  CStochIRForwardSkewMC::Initialize(double spot,const vector<long>& simdates,const CVector& fwds,const CVector& discounts,const CMatrix& vols,const CVector& beta,const CVector& volvol,const CVector& meanRevRate,const CVector& meanRevLevel,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes,CMatrix & controlVariatePrices, MlEqZeroCurveHandle	discountCurve,MlEq2DInterpolatorHandle pAccessVolInterp)
{
	m_saveStochasticDiscounts = true;
	CForwardSkewMC::Initialize(spot,simdates,fwds,discounts,vols,beta,volvol,meanRevRate,meanRevLevel,modelParameters,calibTimeIndex,calibVols,pCalibStrikes,controlVariatePrices, pAccessVolInterp);	

	throw("if he is still with the firm Cuin needs to deal with this now");
	initIR(0.4,false);		// ToDo - David
}

/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: Initialize
**	Returns : nothing
**	Comment : 
****************************************************************/

void CStochIRForwardSkewMC::Initialize(double spot,const vector<long>& simdates,const CVector& fwds,const CVector& discounts,const CMatrix& vols,const CMatrix& modelMarketParameters,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes,CMatrix & controlVariatePrices, MlEqZeroCurveHandle	discountCurve,MlEq2DInterpolatorHandle pAccessVolInterp)
{
	m_saveStochasticDiscounts = true;
	CForwardSkewMC::Initialize(spot,simdates,fwds,discounts,vols,modelMarketParameters,modelParameters,calibTimeIndex,calibVols,pCalibStrikes,controlVariatePrices, pAccessVolInterp);

	throw("if he is still with the firm Cuin needs to deal with this now");

	initIR(0.4,false);		// ToDo - David
}


/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: Initialize
**	Returns : nothing
**	Comment : 
****************************************************************/

void CStochIRForwardSkewMC::Initialize(MlEqAsset& asset,const vector<long>& simdates,MlEqStochBetaVolatilityStructure& betaVolStruct,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes, MlEqHullAndWhiteHandle&  phullWhite, MlEqHullAndWhiteHandle& pPayHullWhite)
{
	pHullWhite			= phullWhite;
	m_pPayHullWhite		= pPayHullWhite;
	m_saveStochasticDiscounts = true;

	CForwardSkewMC::Initialize(asset,simdates,betaVolStruct,modelParameters, calibTimeIndex,calibVols,pCalibStrikes);

	throw("if he is still with the firm Cuin needs to deal with this now");

	initIR(0.4,false);	// ToDo - David
}

/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: Initialize
**	Returns : nothing
**	Comment : 
****************************************************************/

void CStochIRForwardSkewMC::Initialize(MlEqAsset& asset,const vector<long>& simdates,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes, MlEqHullAndWhiteHandle&  phullWhite, MlEqHullAndWhiteHandle& pPayHullWhite, double fCorrelation,bool isMinimalHybrid)
{
	pHullWhite					= phullWhite;
	m_pPayHullWhite				= pPayHullWhite;
	m_saveStochasticDiscounts	= true;

	CForwardSkewMC::Initialize(asset,simdates,modelParameters,calibTimeIndex,calibVols,pCalibStrikes);
	initIR(fCorrelation,isMinimalHybrid);
}

/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: Initialize
**	Returns : nothing
**	Comment : 
****************************************************************/

void CStochIRForwardSkewMC::Initialize(MlEqAsset& asset,product& deriv,MlEqStochBetaVolatilityStructure& betaVolStruct,CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes, MlEqHullAndWhiteHandle& phullWhite, MlEqHullAndWhiteHandle& pPayHullWhite)
{
	const std::vector<long> additionalSimDates;
	Initialize(asset,deriv,additionalSimDates,betaVolStruct,fixings,modelParameters, calibTimeIndex,calibVols,pCalibStrikes,phullWhite,pPayHullWhite);

}

/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: Initialize
**	Returns : nothing
**	Comment : 
****************************************************************/


void CStochIRForwardSkewMC::Initialize(MlEqAsset& asset,product& deriv,const std::vector<long>& additionalSimDates,MlEqStochBetaVolatilityStructure& betaVolStruct,CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes, MlEqHullAndWhiteHandle& phullWhite, MlEqHullAndWhiteHandle& pPayHullWhite)
{

	pHullWhite			= phullWhite;
	m_pPayHullWhite		= pPayHullWhite;
	m_saveStochasticDiscounts	= true;

	CForwardSkewMC::Initialize(asset,deriv,additionalSimDates,betaVolStruct,fixings,modelParameters, calibTimeIndex,calibVols,pCalibStrikes);

	throw("Cuin needs to deal with this now");
	initIR(0.4,false); // ToDo - David

}


/***************************************************************
**	Class   : CStochIRForwardSkewMC 
**	Function: Initialize
**	Returns : nothing
**	Comment : 
****************************************************************/

CStochIRForwardSkewMC::CStochIRForwardSkewMC(MlEqConstDateHandle hDate, double spot,const vector<long>& simdates,const CVector& fwds,const CVector& discounts,const CMatrix& vols,const CVector& beta,const CVector& volvol,const CVector& meanRevRate,const CVector& meanRevLevel,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes,CMatrix & controlVariatePrices, MlEqZeroCurveHandle	discountCurve,MlEq2DInterpolatorHandle pAccessVolInterp) : CForwardSkewMC(hDate)
{
	Initialize(spot,simdates,fwds,discounts,vols,beta,volvol,meanRevRate,meanRevLevel,modelParameters,calibTimeIndex,calibVols, pCalibStrikes,controlVariatePrices, discountCurve,pAccessVolInterp);
}


/***************************************************************
**	Class  : CStochIRForwardSkewMC 
**	Routine: GenerateRandomNumbers
**	Returns: double
**	Action : generates MonteCarlo paths
****************************************************************/


void  CStochIRForwardSkewMC::GenerateRandomNumbers(int idate)
{	
	if ( m_randomsAreExternal ){
		return;
	}

	MlEqMonteCarlo::GenerateRandomNumbers(idate);

//  correlate interest rates and spot

	if ( idate == 0 )
	{
		int n_offset,offset,jdate;
		for ( jdate = 0; jdate < m_nDates; jdate++ )
		{
			n_offset =	m_numberOfFactors*jdate+1;
			offset	 =	n_offset+1;
			for (int ipath = 0 ; ipath <m_nPaths; ipath++ ){	
				m_allRandoms[ipath][offset] = m_corr*m_allRandoms[ipath][n_offset]+ m_xcorr*m_allRandoms[ipath][offset];
			}
		}
	}

	return;
}	

/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: Initialize
**	Returns: void
**	Action : 
****************************************************************/

void CMultiAssetForwardSkewMC::Initialize(vector < MlEqAssetHandle >& asset,const std::vector<long>& mcDates,vector < MlEqStochBetaVolatilityStructureHandle >& betaVolStruct,GVector<int>& rateIsStochastic, const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes)
{
	
	m_nAssets = asset.size();

	m_assets.resize(m_nAssets);
	for ( int i = 0 ; i < m_nAssets; i++ ){
		m_assets[i] = asset[i];
	}


	if ( m_nAssets == 0 ){
		throw( " no assets have been set in CMultiAssetForwardSkewMC::Initialize");
	}

	long nToday = asset[0]->GetDateHandle()->GetDate();


	if ( mcDates[0] == nToday ){
		m_nDates = mcDates.size()-1;
	}
	else if ( mcDates[0] > nToday ){
		m_nDates = mcDates.size();
	}


	MlEqZeroCurveHandle zc =	asset[0]->GetPayZeroCurve(true);
	std::string	payCurr = zc->GetName() ;
	for ( int iasset = 1 ; iasset < m_nAssets; iasset++ )
	{
		zc = asset[iasset]->GetPayZeroCurve(true);
		if ( zc->GetName() != payCurr ){
			throw("pay currencies must be identical for all assets");
		}
	}


	CMatrix fixings;
	m_mMCHelper.initialize(mcDates,fixings,nToday);


	Initialize(nToday,asset,rateIsStochastic,modelParameters, calibTimeIndex,calibVols,pCalibStrikes)	;

	m_fwdSkewMC.resize(m_nAssets);

	for ( int iasset = 0 ; iasset < m_nAssets; iasset++ )
	{
		// do a new here
		if ( !rateIsStochastic[iasset] )
		{
			m_fwdSkewMC[iasset] = new CForwardSkewMC(betaVolStruct[iasset]->getDateToDouble());
			m_fwdSkewMC[iasset]->setExternalRandoms(iasset,m_allMultiRandoms,m_path_array);
			m_fwdSkewMC[iasset]->Initialize(*asset[iasset],mcDates,*betaVolStruct[iasset],modelParameters, calibTimeIndex,calibVols,pCalibStrikes);
		}
	}

	if(	m_nDates == 0 )	return;

	m_fwds.resize(m_nDates,m_nAssets);
	for ( int idate = 0 ; idate < m_nDates; idate++ ){
		for( int iasset = 0 ; iasset < m_nAssets; iasset++ ){
			m_fwds[idate][iasset] = m_fwdSkewMC[iasset]->GetForwardForward(idate,0);
		}
	}

	m_termFwds.resize(m_nDates+1,m_nAssets);
	for ( int idate = 0 ; idate < m_nDates+1; idate++ ){
		for ( int iasset = 0 ; iasset < m_nAssets; iasset++ ){
			m_termFwds[idate][iasset] = GetForwardValue(idate,iasset);
		}
	}

	int iasset = 0;
	m_dt.resize(m_nDates);
	for ( int idate = 0 ; idate < m_nDates; idate++ ){
		m_dt[idate] = m_fwdSkewMC[iasset]->Getdt(idate);
	}
	m_sqrtdt.resize(m_nDates);
	for ( int idate = 0 ; idate < m_nDates; idate++ ){
		m_sqrtdt[idate] = sqrt(m_dt[idate]);
	}


	m_t.resize(m_nDates+1);
	m_t[0] = 0.0;

	for ( int idate = 1 ; idate < m_t.getsize(); idate++ ){
		m_t[idate] = m_t[idate-1]+m_dt[idate-1];
	}

}




/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: calculateOptions
**	Returns: void
**	Action : 
****************************************************************/

//void  CMultiAssetForwardSkewMC::calculateOptions(CVector& results,int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor,int iasset)
//{
//		m_fwdSkewMC[iasset]->calculateOptions(results,startIndex,endIndex,Strikes,controlVariate,impliedVolFlag,factor);
//}


/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: calculateOptions
**	Returns: void
**	Action : 
****************************************************************/


//void  CMultiAssetForwardSkewMC::calculateOptions(CVector& results,MlEqAssetHandle asset, int startIndex,int endIndex,vector< MlEqStrikeHandle > & Strikes,int controlVariate,int impliedVolFlag,double factor)
//{
//	int iasset = findAsset(asset);
//	calculateOptions(results,startIndex,endIndex,Strikes,controlVariate,impliedVolFlag,factor,iasset);
//}




/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: Initialize
**	Returns: void
**	Action : 
****************************************************************/

void CMultiAssetForwardSkewMC::Initialize(vector < MlEqAssetHandle >& asset,const std::vector<long>& mcDates,GVector<int>& rateIsStochastic, const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes)
{
	vector < MlEqStochBetaVolatilityStructureHandle > betaVolStruct(asset.size());

	for ( int i = 0 ; i < asset.size(); i++ )
	{
		MlEqVolatilityStructureHandle			pVol =   asset[i]->GetVolatilityStructure();
		MlEqStochBetaVolatilityStructureHandle	pBetaVol(dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*pVol));

		if (!pBetaVol){
			throw(" asset handle must contain stochastic beta volatility structure for monte carlo calculation");
		}

		betaVolStruct[i] = pBetaVol;
	}

//	GVector<int> rateIsStochastic(asset.size());
	Initialize(asset,mcDates,betaVolStruct,rateIsStochastic, modelParameters, calibTimeIndex,calibVols,pCalibStrikes);	
}

/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: Initialize
**	Returns: void
**	Action : 
****************************************************************/

void CMultiAssetForwardSkewMC::Initialize(vector < MlEqAssetHandle >& asset,product& deriv, CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes)
{
	const std::vector<long> additionalSimDates;
	Initialize(asset,deriv,additionalSimDates,fixings,modelParameters,calibTimeIndex,calibVols,pCalibStrikes);
}


/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: simulate
**	Returns: nothing
**	Action : 
****************************************************************/

double CMultiAssetForwardSkewMC::AmericanDigital(CVector& prices,int ipath,int idate,double DownBarrier,double finalSpot,
									 int iasset,int method,double stressForwardBarrierVol,
									 double stressForwardBarrierVolSlope,int npoints,double blend)
{
	
	double res =  m_fwdSkewMC[iasset]->AmericanDigital( prices, ipath, idate, DownBarrier, finalSpot,
									  iasset, method, stressForwardBarrierVol,
									  stressForwardBarrierVolSlope, npoints, blend);

	return res;
}




/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: Initialize
**	Returns: void
**	Action : 
****************************************************************/

void CMultiAssetForwardSkewMC::Initialize(vector < MlEqAssetHandle >& asset,product& deriv,const std::vector<long>& additionalSimDates, CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes)
{
	std::vector<MlEqStochBetaVolatilityStructureHandle> betaVolStruct(asset.size());

	for ( int i = 0 ; i < asset.size(); i++ )
	{
		MlEqVolatilityStructureHandle			pVol =   asset[i]->GetVolatilityStructure();
		MlEqStochBetaVolatilityStructureHandle	pBetaVol(dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*pVol));

		if (!pBetaVol){
			throw(" asset handle must contain stochastic beta volatility structure for monte carlo calculation");
		}

		betaVolStruct[i] = pBetaVol;
	}


	GVector < int> rateIsStochastic(asset.size());		
	for (int i = 0; i < rateIsStochastic.getsize(); i++){
		rateIsStochastic[i] = 0;
	}

	
	Initialize(asset,deriv,additionalSimDates,betaVolStruct,rateIsStochastic,fixings,modelParameters, calibTimeIndex, calibVols,pCalibStrikes);
}


/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: initializeCalibrationToBasket
**	Returns: void
**	Action : 
****************************************************************/
	
void CMultiAssetForwardSkewMC::initializeCalibrationToBasket()
{	
	if ( !m_calibToBasket ){
		return;
	}


	CVector assetWeights(m_nAssets);
	for ( int i = 0 ; i < m_nAssets; i++ ){
		assetWeights[i] = 1.0/m_nAssets;
	}
	m_basketWeights = assetWeights;

//  calibration to the last simulated date
	long maturityDate =	m_mMCHelper.GetMCDate(m_nDates);

	vector < vector < MlEqStrikeHandle > > pVolStrikes;

	pVolStrikes.resize(m_nAssets+1);
//  calibration to at the money spot:
	for ( int i = 0 ; i < m_nAssets+1; i++ ){
		pVolStrikes[i].resize(1);
	}

	int nToday = m_assets[0]->GetDateHandle()->GetDate();

	double basketSpot = 0.0;
	for ( int i = 0 ; i < m_nAssets; i++ ){
		basketSpot += assetWeights[i]*m_assets[i]->GetSpot(nToday);
	}

 
	double basketFwd = 0.0;
	for ( int i = 0 ; i < m_nAssets; i++ ){
		basketFwd += assetWeights[i]*m_assets[i]->GetQuantoForward(nToday,maturityDate,false); 
	}

	pVolStrikes[0][0] = new MlEqStrike(basketSpot);
	for ( int i = 0 ; i < m_nAssets; i++ ){
		pVolStrikes[i+1][0] = new MlEqStrike(m_assets[i]->GetSpot(nToday));
	}

	m_pBasketCalibrationStrike = pVolStrikes[0][0];

	CVector result;
	result = CalculateEuropeanThreeMomentBasket(m_assets,assetWeights,maturityDate,pVolStrikes);

	double mat = m_hDate->GetYearFraction(maturityDate);
	double impliedVol  =  MlEqBSImpliedVol(result[0],basketFwd,mat,m_pBasketCalibrationStrike->m_strike,1.0,1);

	m_calibBasketValue = impliedVol;
}	

/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: CalculateBasket
**	Returns: void
**	Action : 
****************************************************************/
	
void CMultiAssetForwardSkewMC::CalibrateToBasket()
{

	double impliedBskVol = CalculateBasket();

	// shift correlation here re-run initilaize function and iterate
	
//do some work here alex
}

/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: CalculateBasket
**	Returns: void
**	Action : 
****************************************************************/

double	CMultiAssetForwardSkewMC::GetForwardValue(int idate,int iasset)
{
	double res = m_fwdSkewMC[iasset]->GetForwardValue(idate,0);
	return res;
}



/***************************************************************
**	Class   : CForwardSkewMC 
**	Function: 
**	Returns : nothing
**	Comment :  calculates forward vols between setting iStart and iEnd
****************************************************************/
	
int CMultiAssetForwardSkewMC::GetNumberOfVolStates(int ipath,int idate,int iasset)
{
	int n = m_fwdSkewMC[iasset]->GetNumberOfVolStates(ipath,idate,0);
	return  n;
}


/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: CalculateBasket
**	Returns: void
**	Action : 
****************************************************************/

double CMultiAssetForwardSkewMC::GetCorrelation(MlEqAssetHandle& asset1,MlEqAssetHandle& asset2)
{
	MlEqCorrelationMatrixHandle pCorrel = asset1->GetCorrelationMatrix();
	double correl = pCorrel->GetCorrelation(asset1->GetName(),asset2->GetName());					

	if ( m_correlationShift.getsize() == 0 && m_correlationFactor.getsize() == 0 ){
		return correl;
	}

	correl = (correl+m_correlationShift[0])*m_correlationFactor[0];
	return correl;
}

/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: CalculateBasket
**	Returns: void
**	Action : 
****************************************************************/
	
double CMultiAssetForwardSkewMC::CalculateBasket()
{
	MlEqStrike bskStrike;
	MlEqStrike::convertStrikes(bskStrike,*m_pBasketCalibrationStrike);


	double bskPrice = 0.0;
	double bskVal;
	for ( int ipath = 0 ; ipath < m_nPaths; ipath++ )
	{
		bskVal = 0.0;
		for ( int iasset = 0 ; iasset < m_nAssets; iasset++ ){
			bskVal += m_basketWeights[iasset]*(*m_pPath_array)[ipath][iasset][m_nDates+1];
		}

		bskPrice += MlEqMaths::Max(bskVal-bskStrike.m_strike,0.0);
	}

	bskPrice /= m_nPaths;
	return bskPrice;
}


/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: Initialize
**	Returns: void
**	Action : 
****************************************************************/

void CMultiAssetForwardSkewMC::Initialize(vector < MlEqAssetHandle >& asset,product& deriv,const std::vector<long>& additionalSimDates,vector < MlEqStochBetaVolatilityStructureHandle >& betaVolStruct,GVector<int>& rateIsStochastic, CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes)
{		
	bool								bBaseHelperSet = false;
	
	m_nAssets = asset.size();

	m_assets.resize(m_nAssets);
	for ( int i = 0 ; i < m_nAssets; i++ ){
		m_assets[i] = asset[i];
	}

	MlEqZeroCurveHandle zc =	asset[0]->GetPayZeroCurve(true);
	std::string	payCurr = zc->GetName() ;
	for ( int iasset = 1 ; iasset < m_nAssets; iasset++ )
	{
		zc = asset[iasset]->GetPayZeroCurve(true);
		if ( zc->GetName() != payCurr ){
			throw("pay currencies must be identical for all assets");
		}
	}


	int nDates = deriv.m_payoffDates.getsize();

	long nToday = betaVolStruct[0]->getDateToDouble()->GetDate();
	std::vector<long> mDates(nDates);

	for ( int n = 0 ; n < nDates; n++ ){
		mDates[n] = deriv.m_payoffDates[n];
	}


	initializeMcHelper(deriv,fixings,additionalSimDates);
//	m_mMCHelper.initialize(mDates,fixings,nToday);sossos

	bBaseHelperSet = true;

	const std::vector<long>& mcDates = m_mMCHelper.GetMCDates();


	if ( mcDates[0] == nToday ){
		m_nDates = mcDates.size()-1;
	}
	else if ( mcDates[0] > nToday ){
		m_nDates = mcDates.size();
	}

	;

	Initialize(nToday,asset,rateIsStochastic,modelParameters, calibTimeIndex,calibVols,pCalibStrikes)	;

	m_fwdSkewMC.resize(m_nAssets);

	for ( int iasset = 0 ; iasset < m_nAssets; iasset++ )
	{
		// do a new here
		if ( !rateIsStochastic[iasset] )
		{
			m_fwdSkewMC[iasset] = new CForwardSkewMC(betaVolStruct[iasset]->getDateToDouble());
			m_fwdSkewMC[iasset]->setExternalRandoms(iasset, m_allMultiRandoms, m_path_array);
			m_fwdSkewMC[iasset]->Initialize(*(asset[iasset]), deriv,additionalSimDates, *(betaVolStruct[iasset]), fixings, modelParameters, calibTimeIndex, calibVols, pCalibStrikes);
			if (!bBaseHelperSet){
				m_mMCHelper = m_fwdSkewMC[iasset]->m_mMCHelper;
				bBaseHelperSet = true;
			}
		}
	}
	if (!bBaseHelperSet) throw "Base class Monte Carlo helper class has not been set up.";


	if(	m_nDates == 0 )	return;


	m_fwds.resize(m_nDates,m_nAssets);
	for ( int idate = 0 ; idate < m_nDates; idate++ ){
		for( int iasset = 0 ; iasset < m_nAssets; iasset++ ){
			m_fwds[idate][iasset] = m_fwdSkewMC[iasset]->GetForwardForward(idate,0);
		}
	}

	m_termFwds.resize(m_nDates+1,m_nAssets);
	for ( int idate = 0 ; idate < m_nDates+1; idate++ ){
		for ( int iasset = 0 ; iasset < m_nAssets; iasset++ ){
			m_termFwds[idate][iasset] = GetForwardValue(idate,iasset);
		}
	}

	int iasset = 0;
	m_dt.resize(m_nDates);
	for ( int idate = 0 ; idate < m_nDates; idate++ ){
		m_dt[idate] = m_fwdSkewMC[iasset]->Getdt(idate);
	}
	m_sqrtdt.resize(m_nDates);
	for ( int idate = 0 ; idate < m_nDates; idate++ ){
		m_sqrtdt[idate] = sqrt(m_dt[idate]);
	}


	m_t.resize(m_nDates+1);
	m_t[0] = 0.0;

	for ( int idate = 1 ; idate < m_t.getsize(); idate++ ){
		m_t[idate] = m_t[idate-1]+m_dt[idate-1];
	}

}


/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: Initialize
**	Returns: void
**	Action : 
****************************************************************/

void CMultiAssetForwardSkewMC::Initialize(vector < MlEqAssetHandle >& asset,product& deriv,vector < MlEqStochBetaVolatilityStructureHandle >& betaVolStruct,GVector<int>& rateIsStochastic, CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes)
{		
	const std::vector<long> additionalSimDates;
	Initialize(asset,deriv,betaVolStruct,rateIsStochastic, fixings,modelParameters, calibTimeIndex,calibVols,pCalibStrikes);
}


/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: Initialize
**	Returns: void
**	Action : 
****************************************************************/

void CMultiAssetForwardSkewMC::Initialize(long nToday,vector < MlEqAssetHandle >& asset,GVector<int>& rateIsStochastic, const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes)
{
	m_nAssets = asset.size();
	m_numberOfFactorsPerAsset.resize(m_nAssets);	
	
	GVector< long > index(m_nAssets);
	GVector< long > jndex(m_nAssets);
		
	int n,i;
	int dimensionalityFactor = 0;
	for (i = 0 ; i < m_nAssets; i++ )
	{	
		if( rateIsStochastic[i] )
		{			
			n = 3;
			dimensionalityFactor += n;
			m_numberOfFactorsPerAsset[i] = 3;

			if ( i == 0 ){
				index[i] = n-1;
			}else{
				index[i] = index[i-1] + n;
			}
		}	
		else
		{	
			n = 2;
			dimensionalityFactor += n;
			m_numberOfFactorsPerAsset[i] = 2;
			
			if ( i == 0 ){
				index[i] = n-1;
			}else{
				index[i] = index[i-1] + n;
			}
		}
	}

	jndex[0] = 0;
	for ( i = 1 ; i < index.getsize(); i++ ){
		jndex[i] = index[i-1]+1;
	}

	m_numberOfFactors	= dimensionalityFactor;
	m_nPaths			= modelParameters[5];
	m_randomNumberFlag	= modelParameters[12];

	int randomNumberFlag = modelParameters[12];
	if ( randomNumberFlag == 0 ){
		throw("random number flag must be greater than zero");
	}

//  define correlation
	int ifactor=0,jfactor=0,j;
	CMatrix correl(dimensionalityFactor, dimensionalityFactor);
	MlEqCorrelationMatrixHandle pCorrel = asset[0]->GetCorrelationMatrix();

	for ( ifactor = 0; ifactor < dimensionalityFactor; ifactor++ )
	{
		MlEqMaths::Locate(jndex,ifactor,i);
		for ( jfactor = 0; jfactor < dimensionalityFactor; jfactor++ )
		{	
			MlEqMaths::Locate(jndex,jfactor,j);			
			if ( ifactor - jndex[i] == 0 )
			{	
				if ( jfactor - jndex[j] == 0 ){
					// enter vol-vol correlation
					correl[ifactor][jfactor] = pCorrel->GetVolVolCorrelation(asset[i]->GetName(),asset[j]->GetName());
				}
				else if ( jfactor - jndex[j] == 1 ){

					// enter vol spot correlation = 0 
				}
				else if ( jfactor - jndex[j] == 2 ){
					// enter vol-rate correlation
					correl[ifactor][jfactor] = 0.0;
				}
			}
			else if ( ifactor - jndex[i] == 1 )
			{
				if ( jfactor - jndex[j] == 0 ){

					// enter spot-vol correlation = 0
				}
				else if ( jfactor - jndex[j] == 1 ){
					// enter spot spot correlation  
					correl[ifactor][jfactor] = GetCorrelation(asset[i],asset[j]);
					//pCorrel->GetCorrelation(asset[i]->GetName(),asset[j]->GetName());					
				}
				else if ( jfactor - jndex[j] == 2 ){
					// enter spot-rates correlation 
					correl[ifactor][jfactor] = 0.0;
				}
			}
			else if ( ifactor - jndex[i] == 2 )
			{
				if ( jfactor - jndex[j] == 0 ){
					// enter rates-vol correlation
					correl[ifactor][jfactor] = 0.0;
				}
				else if ( jfactor - jndex[j] == 1 ){
					// enter rates spot correlation  
					correl[ifactor][jfactor] = 0.0;
				}
				else if ( jfactor - jndex[j] == 2 ){
					// enter rates-rates correlation
					correl[ifactor][jfactor] = 0.0;
				}	
			}						
		}		
	}			

	
	
	for ( ifactor = 0; ifactor < dimensionalityFactor; ifactor++ ){
		correl[ifactor][ifactor] = 1.0;
	}

	initializeRandomGenerator(dimensionalityFactor);
	
	Cholesky cholesky(correl,m_nDates-1);
	GenerateRandomNumbers(cholesky);
	
	
	m_seed			= modelParameters[4];
	m_calibflag		= modelParameters[6];

	int fitsize = 2;
	CVector accuracy(fitsize);
	accuracy[0]				= 0.0005;
	accuracy[1]				= 0.0003;
	m_accuracy				= accuracy;
	m_randomsAreExternal	= 0;
	

	const std::vector<long> mcDates = m_mMCHelper.GetMCDates();

	if ( mcDates[0] == nToday ){
		n = mcDates.size();
	}
	else if ( mcDates[0] > nToday )
	{
		n = mcDates.size()+1;
	}

	m_path_array.resize(m_nPaths,m_nAssets,n);
	m_pPath_array = &m_path_array;

}

/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: GetPathValue
**	Returns: double
**	Action : 
****************************************************************/


void CMultiAssetForwardSkewMC ::generatePaths()
{
	if ( m_nDates == 0 ){
		return;
	}

	long idum = m_seed;
	GeneratePath(m_calibflag,idum,m_accuracy);
}


/***************************************************************
**	Class  : MlEqMonteCarlo 
**	Routine: generatePath
**	Returns: double
**	Action : generates MonteCarlo paths
****************************************************************/

void CMultiAssetForwardSkewMC::GeneratePath(int calibflag,long idum,CVector& accuracy)
{
	for ( int i = 0 ; i < m_nAssets; i++ ){
		m_fwdSkewMC[i]->GeneratePath(calibflag,idum,accuracy);
	}
}

/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: CMultiAssetForwardSkewMC
**	Returns: nothing
**	Action : 
****************************************************************/


CMultiAssetForwardSkewMC::CMultiAssetForwardSkewMC(MlEqConstDateHandle hDate, vector < MlEqAssetHandle >& asset,const std::vector<long>& mcDates,vector < MlEqStochBetaVolatilityStructureHandle >& betaVolStruct,GVector<int>& rateIsStochastic, const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes) : MlEqMonteCarlo(hDate)
{
	Initialize(asset,mcDates,betaVolStruct,rateIsStochastic,modelParameters, calibTimeIndex,calibVols,pCalibStrikes);
}


/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: CMultiAssetForwardSkewMC
**	Returns: nothing
**	Action : 
****************************************************************/


CMultiAssetForwardSkewMC::CMultiAssetForwardSkewMC(MlEqConstDateHandle hDate, vector < MlEqAssetHandle >& asset,const std::vector<long>& mcDates,GVector<int>& rateIsStochastic, const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes) : MlEqMonteCarlo(hDate)
{
	Initialize(asset,mcDates,rateIsStochastic,modelParameters,calibTimeIndex,calibVols,pCalibStrikes);
}

/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: CMultiAssetForwardSkewMC
**	Returns: nothing
**	Action : 
****************************************************************/


CMultiAssetForwardSkewMC::CMultiAssetForwardSkewMC(MlEqConstDateHandle hDate, vector < MlEqAssetHandle >& asset,product& deriv,CMatrix& fixings,const CVector& modelParameters,CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes) : MlEqMonteCarlo(hDate)
{
	Initialize(asset,deriv,fixings,modelParameters,calibTimeIndex,calibVols,pCalibStrikes);
}

/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: CMultiAssetForwardSkewMC
**	Returns: nothing
**	Action : 
****************************************************************/


CMultiAssetForwardSkewMC::CMultiAssetForwardSkewMC(MlEqConstDateHandle hDate, vector < MlEqAssetHandle >& asset,product& deriv,vector < MlEqStochBetaVolatilityStructureHandle >& betaVolStruct,GVector<int>& rateIsStochastic, CMatrix& fixings,const CVector& modelParameters, CMatrix& calibTimeIndex,CMatrix& calibVols,vector<vector<MlEqStrikeHandle> >& pCalibStrikes) : MlEqMonteCarlo(hDate)
{
	Initialize(asset,deriv,betaVolStruct,rateIsStochastic,fixings,modelParameters, calibTimeIndex,calibVols,pCalibStrikes);
}


/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: CMultiAssetForwardSkewMC
**	Returns: nothing
**	Action : 
****************************************************************/


const	CVector& CMultiAssetForwardSkewMC::GetDiscounts(int ipath)// gets discounts along one monteCarlo path
{
	return m_fwdSkewMC[0]->GetDiscounts(ipath);
}

/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: CMultiAssetForwardSkewMC
**	Returns: nothing
**	Action : 
****************************************************************/

double	CMultiAssetForwardSkewMC::GetDiscount(int ipath,int idate)
{
	return m_fwdSkewMC[0]->GetDiscount(ipath,idate);
}

/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: CMultiAssetForwardSkewMC
**	Returns: nothing
**	Action : 
****************************************************************/

double	CMultiAssetForwardSkewMC::GetDiscount(int ipath,int idate,int idateAsOf)
{
	return m_fwdSkewMC[0]->GetDiscount(ipath,idate,idateAsOf);
}

/***************************************************************
**	Class  : CMultiAssetForwardSkewMC 
**	Routine: Initialize
**	Returns: void
**	Action : 
****************************************************************/

double  CMultiAssetForwardSkewMC::GetBridgeVolatility(int ipath,int idate,double strike,int iasset)
{
	double res;
	res = m_fwdSkewMC[iasset]->GetBridgeVolatility(ipath,idate,strike);
	return res;
}


/***************************************************************
**	Class  : CSequentialBlackMC 
**	Routine: CSequentialBlackMC
**	Returns: nothing
**	Action : 
****************************************************************/

CMatrix& CSequentialBlackMC::GetPathArray(int ipath)
{


	m_randomGenerator->generateRandoms(m_randoms,m_numberOfFactors,m_cholesky,ipath,m_rndTemp);
			
	double x;
	for ( int iasset = 0 ; iasset < m_nAssets; iasset++ )
	{

		for ( int idate = 0 ; idate < m_nDates; idate++ )
		{
			x = m_vols[iasset][idate]*m_randoms[iasset][idate]*m_sqrtdt[idate];

			m_currentPathValues[iasset][idate+1] = 

				m_currentPathValues[iasset][idate]*m_fwdfwds[iasset][idate]*
				exp(-0.5*m_vols[iasset][idate]*m_vols[iasset][idate]*m_dt[idate]+x);

		}
	}

	return m_currentPathValues;
}


/***************************************************************
**	Class  : CSequentialBlackMC 
**	Routine: CSequentialBlackMC
**	Returns: nothing
**	Action : 
****************************************************************/

void CSequentialBlackMC::initialize(const vector<MlEqAssetHandle>& assets, product& deriv, CMatrix& fixings, CMatrix& correl,int npath,int randomNumberFlag)
{

	int i;

	m_nPaths = npath;
	m_nAssets = assets.size();

	if ( m_nAssets == 0 ){
		throw("no assets have been entered");
	}

//  check whether pay currencies are consistent

	MlEqZeroCurveHandle zc =	assets[0]->GetPayZeroCurve(true);
	std::string	payCurr = zc->GetName() ;
	for ( int iasset = 1 ; iasset < m_nAssets; iasset++ )
	{
		zc = assets[iasset]->GetPayZeroCurve(true);
		if ( zc->GetName() != payCurr ){
//			throw("pay currencies must be identical for all assets");sos
		}
	}

	MlEqVolatilityStructureHandle pVol =   assets[0]->GetVolatilityStructure();	

	initializeMcHelper(deriv,fixings);
	const std::vector<long>& mcDatesVec		= m_mMCHelper.GetMCDates();
	const std::vector<long>& payoffDatesVec	= m_mMCHelper.GetPayoffDates();



	GVector<long> mcDates(mcDatesVec.size());
	GVector<long>  payoffDates(payoffDatesVec.size());

	for ( int i = 0 ; i < mcDates.getsize(); i++ ){
		mcDates[i] = mcDatesVec[i];
	}

	for ( int i = 0 ; i < payoffDates.getsize(); i++ ){
		payoffDates[i] = payoffDatesVec[i];
	}

	if ( mcDates[0] != m_hDate->GetDate() ){
		throw("inconsistent today's date entered");
	}

	if ( mcDates[0] == payoffDates[0] ){
		m_nDates = mcDates.getsize()-1;
	}
	else{
		m_nDates = mcDates.getsize();
	}

	m_dt.resize(m_nDates);
	m_t.resize(m_nDates+1);
	m_sqrtdt.resize(m_nDates);

	m_t[0] = 0.0;
	for ( i = 1 ; i <= m_nDates; i++ ){
		m_t[i] =  m_hDate->GetYearFraction(mcDates[i]);		
	}

	for ( i = 0 ; i < m_nDates; i++ ){
		m_dt[i]		=  m_t[i+1]-m_t[i];	
		m_sqrtdt[i] =  sqrt(m_dt[i]); 
	}


	int numberOfPeriods = mcDates.getsize();
	m_cholesky.initialize(correl,numberOfPeriods);

	m_numberOfFactors.resize(m_nAssets);
	for ( i = 0 ; i < m_nAssets; i++ ){
		m_numberOfFactors[i] = 1;
	}

	m_randoms.resize(m_nAssets);
	for ( i = 0 ; i < m_nAssets; i++ ){
		m_randoms[i].resize(m_nDates);
	}

	m_spots.resize(m_nAssets);

	int nToday = m_hDate->GetDate();
	for ( i = 0 ; i < m_nAssets; i++ ){
		m_spots[i] = assets[i]->GetSpot(nToday);
	}

	double fwd;
	m_vols.resize(m_nAssets,m_nDates);
	for ( i = 0; i < m_nAssets; i++ )
	{
		fwd = assets[i]->GetForward(nToday,mcDates[i+1], false);  
		MlEqStrike strike(fwd);
		for ( int idate = 0 ; idate < m_nDates; idate++ ){
			m_vols[i][idate] = assets[i]->GetCompositeVolatility(strike,mcDates[i],mcDates[i+1]);
		}
	}

	m_currentPathValues.resize(m_nAssets,m_nDates+1);

	for (  i = 0 ; i < m_nAssets; i++ ){
		m_currentPathValues[i][0] =	 m_spots[i];
	}

	m_fwdfwds.resize(m_nAssets,m_nDates);

	for ( i = 0 ; i < m_nAssets; i++ )
	{
		m_fwdfwds[i][0] = assets[i]->GetForward(mcDates[1], false)/assets[i]->GetSpot(mcDates[0]);
		for ( int idate = 1 ; idate < m_nDates; idate++ ){
			m_fwdfwds[i][idate] = assets[i]->GetForward(mcDates[idate+1], false)/assets[i]->GetForward(mcDates[idate], false);
		}

	}

//  set up discounts


	CVector ctimes;
	int i0=0;
	if ( m_hDate->GetYearFraction(mcDates[0]) > 1e-3 ){
		ctimes.resize(mcDates.getsize());
		m_discount.resize(mcDates.getsize());
	}
	else{
		i0=1;
		ctimes.resize(mcDates.getsize()-1);
		m_discount.resize(mcDates.getsize()-1);
	}
	

	for ( i = 0 ; i < ctimes.getsize(); i++ )
	{
		ctimes[i]		=	m_hDate->GetYearFraction(mcDates[i+i0]);
		m_discount[i]	=	zc->GetDiscountFactor(m_hDate->GetDate(),mcDates[i+i0]);
	}


	m_discountPrepended.resize(m_nDates+1);
	for ( int idate = 0 ; idate < m_nDates; idate++ ){
		m_discountPrepended[idate+1] = m_discount[idate];
	}


	int dimensionalityFactor = m_nAssets;
	int ndim = dimensionalityFactor*m_nDates;
	int m_randomNumberFlag = randomNumberFlag;

	CUnitcube* D = new CUnitcube(ndim);
			//	select alpha
			//	CNiederreiteralpha *alpha = new CNiederreiteralpha;

	CAlpha* alpha;
	if ( m_randomNumberFlag == 2 ){
				alpha = new CPrimerootalpha;
	}
	else if ( m_randomNumberFlag == 3 ){
				alpha = new CNiederreiteralpha;
	}
	else if ( m_randomNumberFlag == 4 ){
				alpha = new CBakeralpha;
	}
	else if ( m_randomNumberFlag == 5 ){
				alpha = new CPrimerootalpha;
	}
	else{
				throw("randomNumberFlag must be between 1 and 5");
	}

	CArithmeticmeanweights* weights = new CArithmeticmeanweights;// weight
	int d = 1;// d >= r 
	CBakerperiodization* periodization=new CBakerperiodization;
	CNTparameters* par = new CNTparameters(alpha,weights,periodization);
	CNTintegrator* NTintegrator= new CNTintegrator;

	m_randomGenerator = new Diophantine(ndim,m_nPaths,alpha ,
								weights,D,periodization,NTintegrator,par);

}

/***************************************************************
**	Class  : CHermiteMC 
**	Routine: CHermiteMC
**	Returns: initialize
**	Action : 
****************************************************************/


void CHermiteMC::initialize(GVector < HermiteCoeff > & HermiteCoeff,const vector<MlEqAssetHandle>& assets, product& deriv, CMatrix& fixings, CMatrix& correl,int npath,int randomNumberFlag,int ngauss,int nstdev)
{
	CSequentialBlackMC::initialize(assets, deriv, fixings,correl,npath,randomNumberFlag);

	if ( m_nAssets > 1 ){
		throw("this case has not been impülemented yet horse");
	}

	m_HermiteCoeff = HermiteCoeff;
	createTransitionCoeff();

	m_ngauss = ngauss;
	m_gaussPoints.resize(ngauss);
	m_gaussWeights.resize(ngauss);
	
	double lower = -nstdev;
	double upper = nstdev;

	MlEqMaths::dGauleg(lower,upper,m_gaussPoints,m_gaussWeights,ngauss,true);			

	m_sq2pi = sqrt(2.0*3.141592654);

}

/***************************************************************
**	Class  : CHermiteMC 
**	Routine: CHermiteMC
**	Returns: initialize
**	Action : 
****************************************************************/


void CHermiteMC::initialize(const vector<MlEqAssetHandle>& assets, product& deriv, CMatrix& fixings, CMatrix& correl,int npath,int randomNumberFlag,int ngauss,int nstdev)
{

	CSequentialBlackMC::initialize(assets, deriv, fixings,correl,npath,randomNumberFlag);

	const std::vector<long>& mcDates = 	GetMCDates();

	MlEqVolatilityStructureHandle hVol = assets[0]->GetVolatilityStructure();

	vector <MlEqVolDataHandle >	hVolData = 	hVol->getMidVolData();
	const vector<DATE>& VolDates = hVol->getDates();

	if ( VolDates.size() != hVolData.size() ){
		throw("are number of vol dates out of line with number of voldataobjects?");
	}

	int nsize;
	GVector < MLHermiteVolDataHandle >  hermiteVolData(hVolData.size());
	for ( int i = 0 ; i < hVolData.size(); i++ )
	{
		hermiteVolData[i] = dynamic_cast<MLHermiteVolData*>(&*hVolData[i]);	
		if ( hermiteVolData[i]  == NULL ){
			throw(" volatility in asset was not made of hermite vol data objects");
		}

		const HermiteCoeff hermiteCoeff = hermiteVolData[i]->getHermiteData(0);
		if ( i ){
			if (  hermiteCoeff.hCoeff.getsize() != nsize ){
				throw("currently each volslice object needs to have the same number of entries: sorry");
			}
		}
		else{
			nsize = hermiteCoeff.hCoeff.getsize();
		}

	}

	CVector x(VolDates.size());
	for ( int i = 0 ; i < x.getsize(); i++ ){
		x[i] = VolDates[i];
	}

	CVector y(mcDates.size());
	GVector < HermiteCoeff >  HermiteCoeff(mcDates.size());

	MlEqConstDateHandle		hDate = assets[0]->GetDateHandle();

	MlEqDateHandle nToday = new MlEqDate(*hDate);
	for ( int jdate = 0 ; jdate < mcDates.size(); jdate++ )
	{

		HermiteCoeff[jdate].mat = hDate->GetYearFraction(mcDates[jdate]);
		HermiteCoeff[jdate].fwd = assets[0]->GetForward(hDate->GetDate(),mcDates[jdate],false);
		HermiteCoeff[jdate].nToday = nToday;
		HermiteCoeff[jdate].maturity = mcDates[jdate];
		HermiteCoeff[jdate].hCoeff.resize(nsize);
	};


	for ( int k = 0 ; k < nsize; k++ )
	{
		for ( int idate = 0 ; idate < hermiteVolData.getsize(); idate++ ){
			y[idate] = hermiteVolData[idate]->getHermiteData(0).hCoeff[k];
		}
		
		for ( int jdate = 0 ; jdate < mcDates.size(); jdate++ ){
			HermiteCoeff[k].hCoeff[k] = MlEqMaths::linearInterp(x, y,mcDates[jdate]);
		}
	}


	if ( m_nAssets > 1 ){
		throw("this case has not been impülemented yet horse");
	}

	m_HermiteCoeff = HermiteCoeff;
	createTransitionCoeff();

	m_ngauss = ngauss;
	m_gaussPoints.resize(ngauss);
	m_gaussWeights.resize(ngauss);
	
	double lower = -nstdev;
	double upper = nstdev;

	MlEqMaths::dGauleg(lower,upper,m_gaussPoints,m_gaussWeights,ngauss,true);			

	m_sq2pi = sqrt(2.0*3.141592654);

}


/***************************************************************
**	Class  : CHermiteMC 
**	Routine: CHermiteMC
**	Returns: initialize
**	Action : 
****************************************************************/

double CHermiteMC::Wiener(CVector& currentHermite,int idate,double omega,CVector& prevHermite)
{

	double nfac = 1.0;
	double z		;
	double dt	= m_dt[idate];
	double w	= omega*sqrt(m_dt[idate]);

	double prevVal=0.0,prevPrevVal=0.0,res = 0.0,currentVal;


	for ( int i = 0 ; i < m_hermite; i++ )
	{
		z = 0.0;
		for ( int k = 0 ; k < m_hermite-i; k++ ){
			z += m_b[idate][i][k]*prevHermite[k];
		}

		currentVal = MlEqMaths::Hermite(i,w,dt,prevVal,prevPrevVal);
		currentHermite[i]	= currentVal;

		z *= currentVal/nfac;
		res += z;

		prevPrevVal = prevVal;
		prevVal		= currentVal;
		nfac *= (i+1);
	}

	return res;
}

/***************************************************************
**	Class  : CHermiteMC 
**	Routine: CHermiteMC
**	Returns: initialize
**	Action : 
****************************************************************/

double CHermiteMC::normalizeWiener(int idate,CVector& prevHermite)
{

	double res = 0.0,z,wiener;
	CVector currentHermite(m_hermite);

	for ( int i = 0 ; i < m_ngauss; i++ )
	{
		z = m_gaussWeights[i]/m_sq2pi*exp(-0.5*m_gaussPoints[i]*m_gaussPoints[i]);

		wiener = Wiener(currentHermite,idate,m_gaussPoints[i],prevHermite);
		res += z*exp(wiener);
	}

	return res;
}

/***************************************************************
**	Class  : CHermiteMC 
**	Routine: CHermiteMC
**	Returns: initialize
**	Action : 
****************************************************************/

CMatrix& CHermiteMC::GetPathArray(int ipath)
{

	m_randomGenerator->generateRandoms(m_randoms,m_numberOfFactors,m_cholesky,ipath,m_rndTemp);

	int iasset = 0;			
	double wiener,omega,normalization;

	CVector prevHermite(m_hermite),currentHermite(m_hermite);

	for ( int idate = 0 ; idate < m_nDates; idate++ )
	{

		omega	= m_randoms[iasset][idate];
		wiener	= Wiener(currentHermite,idate,omega,prevHermite);

		normalization = normalizeWiener(idate,prevHermite);

		m_currentPathValues[iasset][idate+1] = 
		m_currentPathValues[iasset][idate]*m_fwdfwds[iasset][idate]*exp(wiener)/normalization;

		prevHermite = currentHermite;

	}

	return m_currentPathValues;

}

/***************************************************************
**	Class  : CHermiteMC 
**	Routine: CHermiteMC
**	Returns: createTransitionCoeff
**	Action : 
****************************************************************/


void CHermiteMC::createTransitionCoeff()
{
	// determine m_hermite
	if ( m_HermiteCoeff.getsize() != m_nDates ){
		throw("incorrect dimensioning of hermite coefficient array");
	}

	m_hermite = m_HermiteCoeff[0].hCoeff.getsize();
	for ( int n = 1; n < m_HermiteCoeff.getsize(); n++ )
	{
		int m = m_HermiteCoeff[n].hCoeff.getsize();
		if ( m > m_hermite ){
			m_hermite = m;
		}
	}

	double dt1,dt2;
	m_b.resize(m_nDates,m_hermite);
	for ( int idate = 0 ; idate < m_nDates; idate++ )
	{
		dt2 = m_dt[idate];
		for ( int n = 0 ; n < m_hermite; n++ )
		{
			m_b[idate][n].resize(m_hermite-n);

			if ( idate == 0 ){
				m_b[idate][n] = m_HermiteCoeff[idate].getHermiteCoeff(n);
			}
			else
			{
				dt1 = m_dt[idate-1];
				for ( int k = 0 ; k < m_hermite-n; k++ )
				{
					double z;
					if ( k == 0 )
					{
						double x = 2.0/(double)(n+1);
						z  = pow(m_HermiteCoeff[idate].getHermiteCoeff(n),x)*(dt1+dt2);
						z -= pow(m_b[idate-1][n][k],x)*dt1;

						if ( n%2 == 0 && z < 0 ){
							throw(" U are in trouble when creating forward hermite parameters");
						}
						z /= dt2;
						z = pow(z,(double)n/2.0);
						m_b[idate][n][k] = z;
					}
					else
					{
						double x= (double)n/((double)(n+k));
						z = pow(m_b[idate-1][n+k][0],x)*pow(m_b[idate][n+k][0],x);
						m_b[idate][n][k] = MlEqMaths::binomial(n+k,k)*z;
					}
				}
			}
		}
	}
}


/***************************************************************
**	Class  : CQuasiMC
**	
**	This sets up a 1 path multiasset Monte Carlo where the
**	path is completely determined by past spot fixings
**	
****************************************************************/

void CQuasiMC::Initialize(std::vector<MlEqSpotScheduleHandle>& ahSpotSchedules, const std::vector<long>& anDates)
{
	m_nDates = anDates.size() - 1; 	// -1 Is here for the first date is today's date kludge in other Monte Carlos.
	m_nPaths = 1;
	m_nAssets = ahSpotSchedules.size();
	m_path_array.resize(m_nPaths, m_nAssets, m_nDates + 1);
	m_pPath_array = &m_path_array;
	m_numberOfFactors = 0;
	
	for (int nSchedule = 0; nSchedule < ahSpotSchedules.size(); nSchedule++){
		std::vector<double>	afSpots;
		ahSpotSchedules[nSchedule]->GetValuesBeforeDate(m_hDate->GetDate(), anDates, &afSpots, ahSpotSchedules[nSchedule]->GetName());
		for (int nDate = 0; nDate < afSpots.size(); nDate++){
			m_path_array[0][nSchedule][nDate] = afSpots[nDate];
		}
	}

	m_discount.resize(m_nDates + 1, 1.0);
	m_discountPrepended.resize(m_nDates + 1, 1.0);
}







////////////////////////////////////////////////////////////////////
//
//
//				 start local vol 
//
//
////////////////////////////////////////////////////////////////////		
	



/***************************************************************
**	Class  : LocalVolGrid 
**	Routine: initialise
**	Returns: 
**	Action : 
****************************************************************/


void LocalVolGrid::initialise(MlEqAsset asset,GVector< GVector < MlEqStrikeHandle > >& calibStrikes,GVector<long> calibDates,GVector < MlEqInterpolatorHandle >& localVolGridInterp,long horizonDate,long stepsPerDayRate,int numberSpacialPoints,int ngauss)
{	
	
	m_mktVols			  = asset.GetVolatilityStructure();
	m_calibDates		  = calibDates;
    m_spot				  = m_mktVols->getSpot();
	m_horizonDate		  = horizonDate;
	m_calibStrikes		  = calibStrikes;
	m_stepsPerDayRate	  = stepsPerDayRate;
	m_localVolInterp	  = localVolGridInterp;
	m_ngauss			  = ngauss;


	MlEqConstDateHandle   startdate	= m_mktVols->getDateToDouble();
	
	m_calibTimes.resize(m_calibDates.getsize());	
	for ( int i = 0 ; i < m_calibTimes.getsize(); i++ ){
		m_calibTimes[i]   = startdate->GetYearFraction(m_calibDates[i]);
	}
	
	int nToday   = startdate->GetDate();

	if ( ngauss == 0 )
	{
		int nsteps = (double)(horizonDate-nToday )/((double)stepsPerDayRate);
		GVector < long > tempGrid(nsteps);

		tempGrid[0] = nToday+stepsPerDayRate;
		for ( int i = 1 ; i < nsteps; i++ ){
			tempGrid[i] = tempGrid[i-1] + stepsPerDayRate;
		}

		merge( m_dateGrid, tempGrid,calibDates ); 
		m_timeGrid.resize(m_dateGrid.getsize());

		for ( int i = 0 ; i < m_dateGrid.getsize(); i++ ){
			m_timeGrid[i] = startdate->GetYearFraction(m_dateGrid[i]);
		}
	}
	else
	{

		m_gaussWeight.resize(m_calibDates.getsize(),ngauss);
		m_gaussPoint.resize(m_calibDates.getsize(),ngauss);
		m_timeGrid.resize(m_calibDates.getsize()*ngauss);
		m_dateGrid.resize(m_calibDates.getsize()*ngauss);

		long nStart = nToday;
		long nEnd;
		int k = 0;
		for ( int idate = 0 ; idate < m_calibDates.getsize(); idate++ )
		{
			nEnd = calibDates[idate];

			double t = 0.0;
			if ( idate > 0 ){
				t = m_calibTimes[idate-1];
			}

			MlEqMaths::dGauleg(t,m_calibTimes[idate],m_gaussPoint[idate],m_gaussWeight[idate],ngauss,true);			

			for ( int i = 0 ; i < ngauss; i++ ){
				m_timeGrid[k] =  m_gaussPoint[idate][i];
				m_dateGrid[k] =  nToday+(double)i/(double)(ngauss-1)*(nEnd-nStart);
				k++;
			}
		}	
		
	}

	m_fwdValues.resize(m_timeGrid.getsize());
	for ( int i = 0; i < m_fwdValues.getsize(); i++ ){
		m_fwdValues[i] = asset.GetForward(m_dateGrid[i],false)/m_spot;
	}

}	




/***************************************************************
**	Class  : LocalVolGrid 
**	Routine: calculateMostLikelyPath
**	Returns: mostlikelyPath,mostlikelyVariance
**	Action : 
****************************************************************/


void LocalVolGrid::iterateMostLikelyPath(CMatrix&	mostlikelyPath,CMatrix&	mostlikelyVariance,CVector& impliedVol,
										 GVector < MlEqStrikeHandle > & calibStrikes,int idateCalib,int idateCalibStart,
										 CMatrix&  previous_mostlikelyPath,CMatrix&  previous_mostlikelyVariance,bool use_prev_most_likely_variance,bool use_prev_most_likely_path)
{


	iterateMostLikelyPath(mostlikelyPath,mostlikelyVariance,impliedVol,
						  calibStrikes,idateCalib,idateCalibStart,
					      previous_mostlikelyPath,previous_mostlikelyVariance,use_prev_most_likely_variance,
						  use_prev_most_likely_path,m_localVolInterp);


}

	
/***************************************************************
**	Class  : LocalVolGrid 
**	Routine: calculateMostLikelyPath
**	Returns: mostlikelyPath,mostlikelyVariance
**	Action : 
****************************************************************/


void LocalVolGrid::iterateMostLikelyPath(CMatrix&	mostlikelyPath,CMatrix&	mostlikelyVariance,CVector& impliedVol,
										 GVector < MlEqStrikeHandle > & calibStrikes,int idateCalib,int idateCalibStart,
										 CMatrix&  previous_mostlikelyPath,CMatrix&  previous_mostlikelyVariance,bool use_prev_most_likely_variance,bool use_prev_most_likely_path,GVector < MlEqInterpolatorHandle >& localVolInterp)
{

	if ( idateCalib <= 0 ){
		throw("idatecalib index must start from 1");
	}

	int nsteps = idateCalib-idateCalibStart;

	if ( nsteps <= 0 ){ 
		throw("number of steps must be greate or equal to zero");
	}

	if ( previous_mostlikelyPath.rows() == 0 )
	{
		previous_mostlikelyPath.resize(calibStrikes.getsize(),m_ngauss*nsteps);
		mostlikelyPath.resize(calibStrikes.getsize(),m_ngauss*nsteps);
	}
	else if ( previous_mostlikelyPath.rows() != mostlikelyPath.rows() || 
			  previous_mostlikelyPath.cols() != mostlikelyPath.cols()  )
	{
		throw("iconsistent set up in calculateMostLikelyPath");
	}
	
	
	if ( previous_mostlikelyVariance.rows() == 0 )
	{
		previous_mostlikelyVariance.resize(calibStrikes.getsize(),m_ngauss*nsteps);
		mostlikelyVariance.resize(calibStrikes.getsize(),m_ngauss*nsteps);
	}
	else if ( previous_mostlikelyVariance.rows() != mostlikelyVariance.rows() || 
				  previous_mostlikelyVariance.cols() != mostlikelyVariance.cols()  )
	{
		throw("iconsistent set up in calculateMostLikelyPath");
	}


	int ndate		= m_ngauss*nsteps;
	int nstrikes	= calibStrikes.getsize();
	
	impliedVol.resize(nstrikes);
	
	CVector varT(nstrikes);
	double dt;
	
	CMatrix volvars(nstrikes,ndate);
	double volvar,eps=1e-2;
	
	int kdateCalib = idateCalibStart;//MlEqMaths::Max(idateCalibStart-1,0);

    for ( int istrike = 0 ; istrike < nstrikes; istrike++ )
	{		
		for ( int idate = 0 ; idate < ndate; idate++ )
		{
			if ( idate%m_ngauss == 0 && idate > 0 ){
				kdateCalib++;
			}

			dt = m_timeGrid[kdateCalib+idate];
			if ( idate > 0 ){
				dt -= m_timeGrid[kdateCalib+idate-1];
			}

			double prevSpot;
			if ( use_prev_most_likely_path ){
				prevSpot = previous_mostlikelyPath[istrike][idate];
			}
			else{
				prevSpot = m_fwdValues[idate];
			}

			MlEqInterpolatorHandle pInterp = getObjectFromVector(localVolInterp,kdateCalib);
			double v  = pInterp->getValue(prevSpot);

			double conv = 0.0;			
			if ( use_prev_most_likely_variance )
			{
				double prevVar  = previous_mostlikelyVariance[istrike][idate];

				double vp = pInterp->getValue(prevSpot*(1.0+eps));
				double vm = pInterp->getValue(prevSpot*(1.0-eps));

				conv = 0.5*prevVar*(vp*vp+vm*vm-2.0*v*v)/pow(eps*prevSpot,2.0);
			}
			
			volvar = v*v + conv;
			volvars[istrike][idate]	=  volvar*dt;
			varT[istrike]			+= volvar*dt;						
		}	
	}		
			
    for ( int istrike = 0 ; istrike < calibStrikes.getsize(); istrike++ )
	{
		MlEqStrike strike;
		MlEqStrike::convertStrikes(strike,*calibStrikes[istrike]);


		int kdateCalib = idateCalibStart;//MlEqMaths::Max(idateCalibStart-1,0);


		double var = 0.0;
		for ( int idate = 0 ; idate < ndate; idate++ )
		{

			if ( idate%m_ngauss == 0 && idate > 0 ){
				kdateCalib++;
			}

			dt = m_timeGrid[kdateCalib+idate];
			if ( idate > 0 ){
				dt -= m_timeGrid[kdateCalib+idate-1];
			}

			var	  += volvars[istrike][idate];						
			double alpha_tT = var/varT[istrike];

			double fac = m_fwdValues[(kdateCalib+1)*idate]*
				pow(strike.m_strike/m_spot/m_fwdValues[(kdateCalib+1)*ndate-1],alpha_tT);

			mostlikelyPath[istrike][idate]	   = fac*exp(0.5*alpha_tT*(varT[istrike]-var));
			mostlikelyVariance[istrike][idate] = fac*fac*exp(2.0*alpha_tT*(varT[istrike]-var))
											     -pow(mostlikelyPath[istrike][idate],2.0);

			MlEqInterpolatorHandle pInterp = getObjectFromVector(localVolInterp,kdateCalib);

			double currSpot = mostlikelyPath[istrike][idate];
			double v  = pInterp->getValue(currSpot);
			double vp = pInterp->getValue(currSpot*(1.0+eps));
			double vm = pInterp->getValue(currSpot*(1.0-eps));

			double conv = 0.5*mostlikelyVariance[istrike][idate]*(vp*vp+vm*vm-2.0*v*v)/pow(eps*currSpot,2.0);
	
			impliedVol[istrike] += (v*v+conv)*m_gaussWeight[kdateCalib][idate];

		}
	}



    for ( int istrike = 0 ; istrike < calibStrikes.getsize(); istrike++ )
	{
		double t;
		
		int n = MlEqMaths::Max(idateCalibStart-1,0);

		if ( idateCalibStart > 0 ){
			t = m_calibTimes[idateCalib-1]-m_calibTimes[n];
		}
		else
		{
			t = m_calibTimes[idateCalib-1];
		}
		impliedVol[istrike] = sqrt(impliedVol[istrike]/t);
	}

}


	
/***************************************************************
**	Class  : LocalVolGrid 
**	Routine: localToImplied
**	Returns: impliedVols,mostlikelyPath
**	Action : 
****************************************************************/


void LocalVolGrid::localToImplied(CVector& impliedVols,CMatrix& mostlikelyPaths, int idateCalib,int idateCalibStart,double accuracy,bool returnMostlikelypath)
{
	if ( idateCalib > m_calibStrikes.getsize() ){
		throw("incorrect date index entered in localtoimplied function");
	}

	localToImplied(impliedVols,mostlikelyPaths, m_calibStrikes[idateCalib],idateCalib,idateCalibStart,accuracy,returnMostlikelypath);
}



/***************************************************************
**	Class  : LocalVolGrid 
**	Routine: localToImplied
**	Returns: impliedVols,mostlikelyPath
**	Action : 
****************************************************************/

void LocalVolGrid::localToImplied(CVector& impliedVols,CMatrix& mostlikelyPaths, GVector < MlEqStrikeHandle > & calibStrikes,int idateCalib,int idateCalibStart,double accuracy,bool returnMostlikelypath)
{
	localToImplied(impliedVols,mostlikelyPaths, calibStrikes,idateCalib,idateCalibStart,accuracy,returnMostlikelypath,m_localVolInterp);
}


/***************************************************************
**	Class  : LocalVolGrid 
**	Routine: localToImplied
**	Returns: impliedVols,mostlikelyPath
**	Action : 
****************************************************************/


void LocalVolGrid::localToImplied(CVector& impliedVols,CMatrix& mostlikelyPaths, GVector < MlEqStrikeHandle > & calibStrikes,int idateCalib,int idateCalibStart,double accuracy,bool returnMostlikelypath,GVector < MlEqInterpolatorHandle >& localVolInterp)
{

	if ( idateCalibStart >= idateCalib ){
		throw("idateCalibStart cannot be greater or equal idateCalibEnd");
	}

	int ndate		= m_ngauss;
	int nstrikes	= calibStrikes.getsize();

	CMatrix* p_mostlikelyPath;
	CMatrix* p_previous_mostlikelyPath;
	CMatrix* p_mostlikelyVariance;
	CMatrix* p_previous_mostlikelyVariance;

	CMatrix mostlikelyPath,previous_mostlikelyPath;
	CMatrix mostlikelyVariance,previous_mostlikelyVariance;

	p_mostlikelyPath				= &mostlikelyPath;
	p_previous_mostlikelyPath		= &previous_mostlikelyPath;
	p_mostlikelyVariance			= &mostlikelyVariance;
	p_previous_mostlikelyVariance	= &previous_mostlikelyVariance;

	CVector targetVals(nstrikes);

	bool use_most_likely_variance	= false;
	bool use_prev_most_likely_path	= false;

	CVector impliedVol;
	CVector addonVariance;

	iterateMostLikelyPath(mostlikelyPath,mostlikelyVariance,impliedVol,
						  calibStrikes,idateCalib,idateCalibStart,previous_mostlikelyPath,
						  previous_mostlikelyVariance,
						  use_most_likely_variance,use_prev_most_likely_path,localVolInterp);


	targetVals = impliedVol;

	int max_iter = 15;
	bool suceed;
	for ( int iter = 0 ; iter < max_iter; iter++ )
	{

		CMatrix* temp = p_previous_mostlikelyPath;

		p_previous_mostlikelyPath			= p_mostlikelyPath;
		p_mostlikelyPath					= temp;
		temp								= p_previous_mostlikelyVariance;
		p_previous_mostlikelyVariance		= p_mostlikelyVariance;
		p_mostlikelyVariance				= temp;


		iterateMostLikelyPath(*p_mostlikelyPath,*p_mostlikelyVariance,impliedVol,
							  calibStrikes,idateCalib,idateCalibStart,*p_previous_mostlikelyPath,
							  *p_previous_mostlikelyVariance,
							  use_most_likely_variance,true,localVolInterp);

		suceed = true;
		for ( int istrike = 0 ; istrike < nstrikes; istrike++ )
		{
			if ( fabs(targetVals[istrike]-impliedVol[istrike]) > accuracy ){
				suceed = false;
				break;
			}
		}

		if ( suceed == true ){
			break;
		}
		else{
			targetVals = impliedVol;
		}

	}
				
	if ( suceed == false ){
		throw("localToImpliedVol iteration did not converge");
	}

	impliedVols = impliedVol;

	if ( returnMostlikelypath ){
		mostlikelyPaths = *p_mostlikelyPath;
	}
}		


/***************************************************************
**	Class  : LocalVolGrid 
**	Routine: localToImplied
**	Returns: impliedVols,mostlikelyPath
**	Action : 
****************************************************************/


void LocalVolGrid::createLocalVol(CVector& impliedVols,double tolerance,double accuracy)
{

	int max_iter = 50;
	int ncalib = m_calibStrikes.getsize();
	CMatrix mostlikelyPaths;
	CVector addonVariance;
	
	for ( int idateCalib = 0 ; idateCalib < ncalib; idateCalib++ )
	{	
		int nstrikes = m_calibStrikes[idateCalib].getsize();
		CVector targetVols(nstrikes);
		impliedVols.resize(nstrikes);
		CVector localVols(nstrikes);
		CVector strikes(nstrikes);
		CVector impliedsShort(nstrikes);

		

		if ( idateCalib > 0 )
		{
			for ( int istrike = 0 ; istrike < nstrikes; istrike++ )
			{
				localToImplied(impliedsShort,mostlikelyPaths,m_calibStrikes[idateCalib],
							   idateCalib,0,accuracy,false);
			}
		}
		
		for ( int istrike = 0 ; istrike < nstrikes; istrike++ )
		{
			targetVols[istrike] = m_mktVols->getVol(*m_calibStrikes[idateCalib][istrike],(DATE)m_calibDates[idateCalib]) ;	

			MlEqStrike stk;
			MlEqStrike::convertStrikes(stk,*m_calibStrikes[idateCalib][istrike]);
			strikes[istrike]   = stk.m_strike/m_spot;
		}

		
		int iter,i;
		double t,x;

		i = idateCalib;//MlEqMaths::Max(idateCalib-1,0);
		for (  iter = 0 ; iter < max_iter; iter++ )
		{	
			
			localToImplied(impliedVols,mostlikelyPaths, m_calibStrikes[idateCalib],idateCalib+1,idateCalib,accuracy,false);

			if ( idateCalib > 0 )
			{				
				for ( int istrike = 0 ; istrike < nstrikes; istrike++ )
				{

					if ( idateCalib == 0 ){
						t = 0.0;
					}
					else{
						t = m_calibTimes[idateCalib-1];
					}

					x = pow(impliedsShort[istrike],2.0)*t;
						
					impliedVols[istrike] = x + pow(impliedVols[istrike],2.0)*(m_calibTimes[idateCalib]-t);
					impliedVols[istrike] = sqrt(impliedVols[istrike]/m_calibTimes[idateCalib]);
				}
			}

			
			bool suceed = true;
			
			for ( int istrike = 0 ; istrike < nstrikes; istrike++ )
			{
				if ( fabs(targetVols[istrike]-impliedVols[istrike]) > tolerance )
				{
					suceed = false;
					break;
				}
			}

			if ( suceed == true ){
				break;
			}

//			update local vol surface here
			
			MlEqInterpolatorHandle pInterp = getObjectFromVector(m_localVolInterp,idateCalib);

			for ( int istrike = 0 ; istrike < nstrikes; istrike++ ){
				localVols[istrike] = pInterp->getValue(strikes[istrike]) ;
			}
			
			double fac;
			for ( int istrike = 0 ; istrike < nstrikes; istrike++ )
			{
				fac = targetVols[istrike]/impliedVols[istrike];

				if (idateCalib  > 0 ){
					fac = 1.0+(fac-1.0)*pow(targetVols[istrike]/localVols[istrike],2.0)*m_calibTimes[idateCalib]
							/(m_calibTimes[idateCalib]-m_calibTimes[idateCalib-1]);
				}

				localVols[istrike] *= fac;//pow(fac,1.5);
			}


			pInterp->reinitialize(strikes,localVols,idateCalib);

			if ( iter == max_iter ){
				throw("local vol method did not converge");
			}
		}
	}

}





/***************************************************************
**	Class  : LocalVolGrid 
**	Routine: calculateMostLikelyPath
**	Returns: mostlikelyPath,mostlikelyVariance
**	Action : 
****************************************************************/

/*
void LocalVolGrid::initialise(MlEqVolatilityStructureHandle locaVolGrid,MlEqAsset asset,
						 long horizonDate,int stepsPerDayRate,int numberSpacialPoints)
{

	const vector<DATE>& vcalibDates = locaVolGrid->getDates();
	int n = vcalibDates.size();
	
	int i;
	if ( vcalibDates[n-1] > horizonDate )
	{
		for ( i = 0 ; i < n; i++ ){
			if ( vcalibDates[i] >= horizonDate ){
				break;
			}
		}
		n = i;	
	}
	
	GVector< long > calibDates(n);
	for ( int i = 0 ; i < n; i++ ){
		calibDates[i] = vcalibDates[i];
	}
	
	vector <MlEqVolDataHandle >&	midVols	= locaVolGrid->getMidVolData();
	GVector< GVector < MlEqStrikeHandle > > calibStrikes;
	calibStrikes.resize(n);



	if ( midVols.size() > 1 )
	{
		if ( midVols.size() < n ){
			throw("incorrect set up in voldata initialization");
		}

		for ( int idate = 0 ; idate < n; idate++ )
		{	
			const vector< vector< MlEqStrikeHandle > > & strikes =  midVols[idate]->getStrikes();

			if ( strikes.size() != 1 ){
				throw("incorrect set up in voldata initialization");
			}
		
			calibStrikes[idate].resize(strikes[idate].size());
			
			for ( int istrike = 0 ; istrike < strikes[idate].size(); istrike++ ){
				calibStrikes[idate][istrike] = strikes[idate][istrike];
			}	
		}	
						
	}		
	else if ( midVols.size() == 1 )
	{

		const vector< vector< MlEqStrikeHandle > > & strikes =  midVols[0]->getStrikes();

		if ( strikes.size() < n ){
			throw("incorrect set up in voldata initialization");
		}


		for ( int idate = 0 ; idate < n; idate++ )
		{	
			calibStrikes[idate].resize(strikes[idate].size());		
			for ( int istrike = 0 ; istrike < strikes[idate].size(); istrike++ ){
				calibStrikes[idate][istrike] = strikes[idate][istrike];
			}	
		}	
	}
				
	initialise(locaVolGrid,asset,calibStrikes,calibDates,horizonDate,stepsPerDayRate,numberSpacialPoints);
}	
*/
	
/***************************************************************
**	Class  : LocalVolGrid 
**	Routine: calculateMostLikelyPath
**	Returns: mostlikelyPath,mostlikelyVariance
**	Action : 
****************************************************************/
	
/*	
void LocalVolGrid::initialise(GVector < MlEqInterpolatorHandle >& localVolGrid,GVector< long > & vcalibDates,MlEqAsset asset,
						 long horizonDate,int stepsPerDayRate,int numberSpacialPoints)
{	
	
	m_localVolInterp = localVolGrid ;

//	MlEqConstInterpolatorHandle pInterp = localVolGrid->getYInterpolator();
	
	int n = vcalibDates.getsize();
	
	int i;
	if ( vcalibDates[n-1] > horizonDate )
	{
		for ( i = 0 ; i < n; i++ ){
			if ( vcalibDates[i] >= horizonDate ){
				break;
			}
		}
		n = i;	
	}
	
	GVector< long > calibDates(n);
	for ( int i = 0 ; i < n; i++ ){
		calibDates[i] = vcalibDates[i];
	}
	
	GVector< GVector < MlEqStrikeHandle > > calibStrikes;

	calibStrikes.resize(n);

	double spot = asset.GetSpot(asset.GetCurrentDate());

	int m,k,whichData;
	int ninterp = localVolGrid.getsize();
	
	for ( int i = 0 ; i < n; i++ )
	{

		MlEqConstInterpolatorHandle xInterp = getObjectFromVector(localVolGrid,i);

		const CVector& xstrikes = xInterp->getXData(i);


do some work here

		m = xInterp->getNumberOfSlices();

		if( m )
		{
			k = i;
		}
		else{
			k = 0;
		}



		int nstrikes = xstrikes.getsize();
		calibStrikes[i].resize(nstrikes);
		for ( k = 0 ; k < nstrikes; k++ )
		{
			MlEqStrikeHandle stk= new MlEqStrike(xstrikes[k]*spot);
			calibStrikes[i][k] = stk;
		}
	}
			

	initialise(asset,calibStrikes,calibDates,horizonDate,stepsPerDayRate,numberSpacialPoints);
	
	
}	
*/



/**********************************************
***											***
***			Local Vol Monte Carlo			***
***											***
***********************************************/

// still needs debugging...


CLocalVolMC::CLocalVolMC(MlEqConstDateHandle hDate, MlEqAssetHandle hUdly)
:CSequentialBlackMC(hDate)
{
	m_nAssets = 1;
	m_assets.resize(1);
	m_assets[0] = hUdly;

	m_storedPaths = false;

	m_computeGreeks = false;
}


// here fixings should be [iasset][nPastFixings]

void CLocalVolMC::initialize( product& deriv, const CMatrix& Fixings, int npath, int rngFlag,int stepAYear)
{
	long nToday = m_hDate->GetDate();
	std::vector<long>	payoffDates;
	int n = (deriv.m_payoffDates).getsize() ;
	for(int i=0; i<n; i++){
		payoffDates.push_back( deriv.m_payoffDates[i] ) ;
	}

	initialize( payoffDates, Fixings, npath, rngFlag, stepAYear );

	m_mMCHelper.initialize( payoffDates, Fixings, nToday, payoffDates );	// for consistency
}


void CLocalVolMC::initialize( const std::vector<long>& payoffDates, int npath, int rngFlag,int stepAYear)
// this is the handle case...add today's date and use no past fixings
{
	long nToday = m_hDate->GetDate();

		// first insert today in the list if it is not there already
	std::set<long> horse;	
	horse.insert(nToday);	
	for(int i=0; i<payoffDates.size(); ++i)	{
		horse.insert(payoffDates[i]);		
	}

	std::vector<long> c_payoffDates;		// make a copy with today as an extra date
	for (std::set<long>::const_iterator it = horse.begin(); it != horse.end(); ++it){
		c_payoffDates.push_back(*it);
	}
	// these 2 vectors make sure that today's date is in the Monte Carlo 
	// and that the number of simulation dates doesn't include today...

	CMatrix fixings(m_nAssets, 0);
	initialize( c_payoffDates, fixings, npath, rngFlag,stepAYear);

	m_mMCHelper.initialize( c_payoffDates, CMatrix(), nToday, c_payoffDates );	// for consistency
}


//=====================================================
//
//		These 2 versions perform calibration...
//
//=====================================================
/*
void CLocalVolMC::initialize( product& deriv, const CMatrix& Fixings, int npath, int rngFlag )
{
	long nToday = m_hDate->GetDate();
	std::vector<long>	payoffDates;
	int n = (deriv.m_payoffDates).getsize() ;
	for(int i=0; i<n; i++){
		payoffDates.push_back( deriv.m_payoffDates[i] ) ;
	}

	initialize( payoffDates, Fixings, npath, rngFlag );

	m_mMCHelper.initialize( payoffDates, Fixings, nToday, payoffDates );	// for consistency
}


void CLocalVolMC::initialize( const std::vector<long>& payoffDates, int npath, int rngFlag )
// this is the handle case...add today's date and use no past fixings
{
	long nToday = m_hDate->GetDate();

		// first insert today in the list if it is not there already
	std::set<long> horse;	
	horse.insert(nToday);	
	for(int i=0; i<payoffDates.size(); ++i)	{
		horse.insert(payoffDates[i]);		
	}

	std::vector<long> c_payoffDates;		// make a copy with today as an extra date
	for (std::set<long>::const_iterator it = horse.begin(); it != horse.end(); ++it){
		c_payoffDates.push_back(*it);
	}
	// these 2 vectors make sure that today's date is in the Monte Carlo 
	// and that the number of simulation dates doesn't include today...

	CMatrix fixings(m_nAssets, 0);
	initialize( c_payoffDates, fixings, npath, rngFlag );

	m_mMCHelper.initialize( c_payoffDates, CMatrix(), nToday, c_payoffDates );	// for consistency
}
*/
// end comment

void CLocalVolMC::initialize( const std::vector<long>& payoffDates, const CMatrix& Fixings, int npath, int rngFlag , int stepAYear)
{
	m_nAssets = 1;

	initTimeGrids( payoffDates, stepAYear );
	initDiscounts( payoffDates, Fixings);
	initLocalVolGrid();
	initRandomGenerator( npath, rngFlag );
}

/*
//void CLocalVolMC_test::initialize( const std::vector<long>& payoffDates, const CMatrix& Fixings, int npath, int rngFlag,	int stepAYear )
void CLocalVolMC::initialize( const std::vector<long>& payoffDates, const CMatrix& Fixings, int npath, int rngFlag )
{
// calibration is made and the time grid therefore built in the routine
	m_nAssets = 1;

	initLocalVolGrid(payoffDates);

	initTimeGrids( payoffDates  );
	initDiscounts( payoffDates, Fixings);
	
	initRandomGenerator( npath, rngFlag );
}
*/

void CLocalVolMC::initRandomGenerator(int npath, int rngFlag)
{
	if( rngFlag == 6 ){
		throw "You musn't use a Mersenne Twister for this model !!!";
	}

	m_randomNumberFlag = rngFlag;
	m_nPaths = npath;	

	m_numberOfFactorsPerAsset.resize(m_nAssets);
	for(int i=0; i<m_nAssets; ++i){
		m_numberOfFactorsPerAsset[i] = 1;
	}

	m_randomsAreExternal = 0;
	

// initialize quasi with our flag

	initBBvar();

	int nDates = m_nDates;
	m_nDates = m_nBBDates;

	initializeRandomGenerator(m_nAssets);

	m_quasiGenerator = m_randomGenerator ;

// initialize the MT	

	m_nDates = nDates;
	m_randomNumberFlag = 6 ;	

	initializeRandomGenerator(m_nAssets);

	m_randomNumberFlag = rngFlag;

// resize storage stuff

	m_randoms.resize(m_nAssets);
	m_quasis.resize( m_nAssets );
	
	for(int i=0; i<m_nAssets; ++i){
		m_randoms[i].resize(m_nDates);
		m_quasis[i].resize( m_nBBDates );
	}	
}



void CLocalVolMC::initTimeGrids(  const std::vector<long>& payoffDates )
{

	m_nMcDates = payoffDates.size();
	m_mapToFixings.resize(m_nMcDates);
	m_termFwds.resize(m_nMcDates, 1);

	long nToday = m_hDate->GetDate();

	m_nDates = m_dates.getsize()-1;

	long nextDate = nToday; 
	int t=0;

	for(int i=0; i<m_nMcDates; ++i)
	{
		nextDate = payoffDates[i];

		if( nextDate > nToday )
		{
			while( ++t < m_nDates && m_dates[t] <= nextDate ){}		
		}

		m_mapToFixings[i] = t;	
	}


	m_dt.resize(m_nDates);
	m_sqrtdt.resize(m_nDates);

	for (int i = 0 ; i < m_nDates; i++ )
	{
		m_dt[i]			=  m_t[i+1]-m_t[i];	
		m_sqrtdt[i]		=  sqrt(m_dt[i]); 
	}

	m_nStart = 0;
	while( payoffDates[m_nStart] <= nToday ){
		m_nStart++;
	}
}


void CLocalVolMC::initTimeGrids(  const std::vector<long>& payoffDates, int stepAYear )
{

	m_nMcDates = payoffDates.size();
	m_mapToFixings.resize(m_nMcDates);
	m_termFwds.resize(m_nMcDates, 1);

	long nToday = m_hDate->GetDate();

	std::vector<long> sDates;
	sDates.push_back(nToday);

	m_nDates = 0;	// here m_nDates is actually the number of steps (or dates in the future to simulate...)

	long nextDate = nToday, prevDate = nToday; 
	double Tprev = 0., T = Tprev;

//  create simulation dates
	for(int i=0; i<m_nMcDates; ++i)
	{
		nextDate = payoffDates[i];

		if( nextDate > nToday )
		{

			T = m_hDate->GetYearFraction(nextDate) ;			
			int nSteps = int( stepAYear*(T-Tprev) ) + 1;
			double dt = (nextDate - prevDate) / double(nSteps) ;
			if(dt < 1.)	throw "You want too many steps";

			for (int j = 1 ; j < nSteps; j++ ) // start at 1 because we already included sim date prev at the previous i 
			{
				long date = long(prevDate + j*dt) ;
				sDates.push_back(date);
			}
			sDates.push_back(nextDate);

			m_nDates += nSteps ;			
		
			prevDate = nextDate;
			Tprev = T;
		}

		m_mapToFixings[i] = m_nDates;	
	}


	m_dt.resize(m_nDates);
	m_t.resize(m_nDates+1);
	m_sqrtdt.resize(m_nDates);
	m_dates.resize(m_nDates+1);

	m_t[0] = 0.0;
	m_dates[0] = sDates[0];
	for (int i = 0 ; i < m_nDates; i++ )
	{
		m_dates[i+1]	=  sDates[i+1];
		m_t[i+1]		=  m_hDate->GetYearFraction(m_dates[i+1]);
		m_dt[i]			=  m_t[i+1]-m_t[i];	
		m_sqrtdt[i]		=  sqrt(m_dt[i]); 
	}

	m_nStart = 0;
	while( payoffDates[m_nStart] <= nToday ){
		m_nStart++;
	}
}

void CLocalVolMC::initDiscounts(const std::vector<long>& payoffDates, const CMatrix& Fixings)
{
	if( Fixings.rows() != 1 ){
		throw "wrong size for fixings il local vol monte carlo";
	}


	MlEqAssetHandle hUnderlying = m_assets[0];

	m_forwards.resize(m_nMcDates);
	m_termFwds.resize(m_nMcDates, 1);
	m_fwds.resize(m_nMcDates, 1);
	m_discount.resize(m_nMcDates);
	MlEqZeroCurveHandle zc =	hUnderlying->GetPayZeroCurve(true);

	long nToday = m_hDate->GetDate();
	double spot = hUnderlying->GetSpot(nToday);


//  create simulation dates
	for(int i=0; i<m_nMcDates; ++i)
	{
		long nextDate = payoffDates[i];

		if( nextDate > nToday )
		{
			m_forwards[i]		=	hUnderlying->GetNaturalForward( nToday, nextDate, false);
			m_termFwds[i][0]	=	hUnderlying->GetQuantoForward( nToday, nextDate, false);
			if( i > 0){
				m_fwds[i][0]	=	m_termFwds[i][0] / m_termFwds[i-1][0] ;
			}
			
			m_discount[i] = zc->GetDiscountFactor(nToday, nextDate);		
		
		}
		else if( i < Fixings.cols() )		// these 2 cases may be exclusive...
		{
			m_forwards[i] = Fixings[0][i];
			m_discount[i] = 0.;///////////
		}
		else
		{
			m_forwards[i] = spot;
			m_discount[i] = 1.; 
		}
	}

	if( !m_computeGreeks )
	{
		m_currentPathValues.resize(m_nAssets, m_nMcDates);

		for(int i=0; i<m_nStart; ++i){
			m_currentPathValues[0][i] = m_forwards[i];
		}	
	}
	else
	{
		m_currentPathValues.resize(4, m_nMcDates);

		for(int i=0; i<m_nStart; ++i)
		{
			double fwd = m_forwards[i];

			for(int k=0; k<4; k++){
				m_currentPathValues[k][i] = fwd;
			}
		}
	}

}


void CLocalVolMC::initLocalVolGrid()
{
	MlEqAssetHandle hUnderlying = m_assets[0];

	long nToday = m_dates[0];
	double forward = hUnderlying->GetForward(nToday, m_dates[m_nDates-1] , false);
	double vol = hUnderlying->GetVolatility( forward, nToday, m_dates[m_nDates-1] );
	double T = m_hDate->GetYearFraction( m_dates[m_nDates-1] );

	double upBound =  4.* vol * sqrt(T) ;
	double dnBound = -6.* vol * sqrt(T) ;

	int nSpace = 500 * sqrt(T);	// need a lot..
	nSpace = (nSpace < 100)? 100:( (nSpace > 1000)? 1000:nSpace ) ;


	RCPtr<DupireLocalVol> lv_computer = new DupireLocalVol ;
	lv_computer->initialize( hUnderlying, dnBound, upBound, m_dates, nSpace );

	m_localvolgrid		= lv_computer->getLocalVolGrid();
	m_reducedSpotGrid	= lv_computer->getSpotGrid();

	m_min_yGrid = m_reducedSpotGrid[0];	// quite ugly but most efficient...
	m_dy		= m_reducedSpotGrid[1] - m_reducedSpotGrid[0];
	m_max_yIndex = m_reducedSpotGrid.getsize() - 1e-5;

	int row = m_localvolgrid.getsize() ;	// shoud be m_nDates+1
	m_space_size = m_reducedSpotGrid.getsize();

	int arraysize =  row * m_space_size ;	// we store the grid in an array because it is much faster
	m_arraylocalvolgrid.resize( arraysize );

	for(int i=0; i<row; ++i)
	{
		int tslice = i * m_space_size;
		for(int j=0; j<m_space_size; ++j){
				m_arraylocalvolgrid[ tslice + j ] = m_localvolgrid[i][j];
		}
	}

	if( m_computeGreeks )
	{
		lv_computer->parallelVegaBump();
		GVector<CVector> lvgrid	= lv_computer->getLocalVolGrid();

		m_arraybumpedlocalvolgrid.resize( arraysize );

		for(int i=0; i<row; ++i)
		{
			int tslice = i * m_space_size;
			for(int j=0; j<m_space_size; ++j){
					m_arraybumpedlocalvolgrid[ tslice + j ] = lvgrid[i][j];
			}
		}
	}

	m_quantoDrift = lv_computer->getQuantoDrift();
}



void CLocalVolMC::initLocalVolGrid( const std::vector<long>& payoffDates )
{
	throw "You don't need this function any longer";


	MlEqAssetHandle hUnderlying = m_assets[0];

	long nToday = m_hDate->GetDate();
	long maturityDate = payoffDates[payoffDates.size()-1];

	double forward = hUnderlying->GetForward(nToday, maturityDate , false);
	double vol = hUnderlying->GetVolatility( forward, nToday, maturityDate );
	double T = m_hDate->GetYearFraction( maturityDate );
	double upBound =  4.* vol * sqrt(T) ;
	double dnBound = -6.* vol * sqrt(T) ;

	int nSpace = 500 * sqrt(T);	// need a lot..
	nSpace = (nSpace < 100)? 100:( (nSpace > 1000)? 1000:nSpace ) ;

	RCPtr<FittedLocalVol> lv_computer = new FittedLocalVol ;
	lv_computer->FittedLocalVol::initialize( hUnderlying, dnBound, upBound, nSpace, payoffDates );

	m_localvolgrid		= lv_computer->getLocalVolGrid();
	m_reducedSpotGrid	= lv_computer->getSpotGrid();
	m_dates				= lv_computer->getDates();
	m_t					= lv_computer->getTimes();

	m_min_yGrid = m_reducedSpotGrid[0];	// quite ugly but most efficient...
	m_dy		= m_reducedSpotGrid[1] - m_reducedSpotGrid[0];
	m_max_yIndex = m_reducedSpotGrid.getsize() - 1e-5;

	int row = m_localvolgrid.getsize() ;	// shoud be m_nDates+1
	m_space_size = m_reducedSpotGrid.getsize();

	int arraysize =  row * m_space_size ;	// we store the grid in an array because it is much faster
	m_arraylocalvolgrid.resize( arraysize );

	for(int i=0; i<row; ++i)
	{
		int tslice = i * m_space_size;
		for(int j=0; j<m_space_size; ++j){
				m_arraylocalvolgrid[ tslice + j ] = m_localvolgrid[i][j];
		}
	}

	
	if( m_computeGreeks )
	{
		lv_computer->parallelVegaBump();
		GVector<CVector> lvgrid	= lv_computer->getLocalVolGrid();

		m_arraybumpedlocalvolgrid.resize( arraysize );

		for(int i=0; i<row; ++i)
		{
			int tslice = i * m_space_size;
			for(int j=0; j<m_space_size; ++j){
					m_arraybumpedlocalvolgrid[ tslice + j ] = lvgrid[i][j];
			}
		}
	}

	m_quantoDrift = lv_computer->getQuantoDrift();
}





long CLocalVolMC::GetNumberOfFutureDates()
{ 
	long n = m_nMcDates;	// useful for the Monte Carlo handle...
	return n-1;
}

const CVector& CLocalVolMC::GetDiscounts(int ipath)
{
	return m_discount;
}


void CLocalVolMC::createPayoutPath(CMatrix& payoutPathArray,CMatrix& simulatedPathArray)
{
	payoutPathArray = simulatedPathArray;
}

void CLocalVolMC::createPayoutPath(CMatrix& payoutPathArray,CMatrix& simulatedPathArray,CVector& payoutPathDiscounts,CVector& simulatedDiscounts)
{
	payoutPathArray = simulatedPathArray;
	payoutPathDiscounts = simulatedDiscounts; // = m_discount
}

	
double CLocalVolMC::GetDiscount(int ipath,int idate)	// for the handle
{
	if ( idate == 0 ){
		return 1.0;
	}

	if ( idate >= m_nMcDates ){
		throw("error out of bound in GetForwardDiscount");
	}

	return	m_discount[m_nStart+idate-1];
}



CMatrix& CLocalVolMC::GetPathArray(int ipath)
{
	if( m_storedPaths ){
		return (*m_pPath_array)[ipath];
	}

// otherwise make the path

	m_randomGenerator->generateRandoms(m_randoms[0], ipath);
	m_quasiGenerator->generateRandoms(m_quasis[0], ipath);
	BrownianBridgePath( m_randoms[0], m_quasis[0] );
	
	if( m_computeGreeks ){
		MakeGreeksPathArray( m_randoms[0] );
	}
	else{		
		MakePathArray( m_randoms[0] );
	}

	return m_currentPathValues;
}


void CLocalVolMC::initBBvar()
{
	std::vector<int> mapToBBfixings;

	long nToday = m_hDate->GetDate();
	long nextDate = nToday, prevDate = nToday; 

	double dt_min = 28.0; // number of days...4w

//  create simulation dates
	for(int i=0; i<m_nMcDates; ++i)
	{
		int idate = m_mapToFixings[i];
		nextDate = m_dates[idate];

		if( nextDate > nToday && double(nextDate - prevDate) >= dt_min)
		{	
			mapToBBfixings.push_back( idate );		
			prevDate = nextDate;
		}
	}

	if( mapToBBfixings.size() == 0 ){
		mapToBBfixings.push_back( m_nDates-1 );
	}

	m_nBBDates = mapToBBfixings.size();
	mapToBBfixings[m_nBBDates-1] = m_nDates-1;

	m_mapToBBfixings.resize(m_nBBDates);
	m_bbPtr.resize(3*m_nDates);
	m_bbVar.resize(m_nBBDates);

	int idate = 0;
	int j = 0;

	for(int bbdate=0; bbdate<m_nBBDates; ++bbdate)
	{
		int nextDate = m_mapToBBfixings[bbdate] = mapToBBfixings[bbdate];
		double T = m_t[nextDate];
		m_bbVar[bbdate] = sqrt(T - m_t[idate]);

		for ( ; idate < nextDate-1; idate++ )
		{
			double dt = m_dt[idate];
			double den = T - m_t[idate];
			double num = den - dt;

			m_bbPtr[j++] = dt / den;
			m_bbPtr[j++] = num / den;
			m_bbPtr[j++] = sqrt( dt*num /den );
		}

		m_bbPtr[j++] = 1.0;	
		m_bbPtr[j++] = 0.0;
		m_bbPtr[j++] = 0.0;

		idate++;
	}
}


void CLocalVolMC::BrownianBridgePath(CVector& randoms, const CVector& quasi)
{
	int idate = 0;
	double rand, prevRand = 0.0;

	double* dW = randoms.getPtr();
	const double* bbPtr = m_bbPtr.getConstPtr();

	for(int fdate=0; fdate<m_nBBDates; ++fdate)
	{
		int nextDate = m_mapToBBfixings[fdate];

		double qRand = quasi[fdate] * m_bbVar[fdate] ;

		for ( ; idate < nextDate; idate++ )
		{	
			rand =  (*bbPtr++) * qRand ;
			rand += (*bbPtr++) * prevRand ;
			rand += (*bbPtr++) * (*dW);
		
			*dW++ = rand - prevRand ;

			prevRand = rand;
		}
		prevRand = 0.0;
	}
}



CMatrix& CLocalVolMC::MakePathArray(CVector& randoms)
{
	double rspot = 	0.;
	int idate = 0;

	const double* dW = randoms.getConstPtr();
	const double* dt = m_dt.getConstPtr();
	const double* q = m_quantoDrift.getConstPtr();

	const double* plv = m_arraylocalvolgrid.getConstPtr() ;

	for(int fdate=m_nStart; fdate<m_nMcDates; ++fdate)
	{
		int nextDate = m_mapToFixings[fdate];

		for ( ; idate < nextDate; idate++ )
		{
			double lvol = computeLocalVol( idate, rspot, plv );	
			double drift = lvol * ( (*q++) -0.5*lvol );

			rspot  +=  drift * (*dt++) + lvol * (*dW++) ;
		}

		m_currentPathValues[0][fdate] = m_forwards[fdate] * exp( rspot ) ;
	}

	return m_currentPathValues;
}


CMatrix& CLocalVolMC::MakeGreeksPathArray(CVector& randoms)
{
// sticky strike greeks

	double rspot = 	0.;
	double rspot_up = 	0.01;
	double rspot_dn =  -0.01;
	double rspot_vup = 	0.;

	int idate = 0;

	const double* g = randoms.getPtr();
	const double* pdt = m_dt.getPtr();
	const double* q = m_quantoDrift.getPtr();

	const double* plv = m_arraylocalvolgrid.getConstPtr() ;
	const double* pblv = m_arraybumpedlocalvolgrid.getConstPtr() ;

	for(int fdate=m_nStart; fdate<m_nMcDates; ++fdate)
	{
		int nextDate = m_mapToFixings[fdate];

		for ( ; idate < nextDate; idate++ )
		{
			double lvol = computeLocalVol( idate, rspot, plv );	

			double lvol_up = computeLocalVol( idate, rspot_up, plv );	
			double lvol_dn = computeLocalVol( idate, rspot_dn, plv );

			double lvol_bp = computeLocalVol( idate, rspot_vup, pblv );

			double dW = *g++ ;
			double drift = *q++ ;
			double dt = *pdt++;

			rspot  +=  lvol * ( ( drift -0.5*lvol )* dt +  dW );

			rspot_up  +=  lvol_up * ( ( drift -0.5*lvol_up )* dt +  dW );
			rspot_dn  +=  lvol_dn * ( ( drift -0.5*lvol_dn )* dt +  dW );

			rspot_vup  +=  lvol_bp * ( ( drift -0.5*lvol_bp )* dt +  dW );
		}

		double fwd = m_forwards[fdate] ;

		m_currentPathValues[0][fdate] = fwd * exp( rspot ) ;

		m_currentPathValues[1][fdate] = fwd * exp( rspot_up ) ;
		m_currentPathValues[2][fdate] = fwd * exp( rspot_dn ) ;

		m_currentPathValues[3][fdate] = fwd * exp( rspot_vup ) ;
	}

	return m_currentPathValues;
}


/*
CMatrix& CLocalVolMC::MakeGreeksPathArray(CVector& randoms)
{
// sticky moneyness greeks

	double rspot = 	0.;
	double rspot_vup = 	0.;

	int idate = 0;

	const double* g = randoms.getPtr();
	const double* pdt = m_dt.getPtr();
	const double* q = m_qadj.getPtr();

	const double* plv = m_arraylocalvolgrid.getConstPtr() ;
	const double* pblv = m_arraybumpedlocalvolgrid.getConstPtr() ;


	for(int fdate=m_nStart; fdate<m_nMcDates; ++fdate)
	{
		int nextDate = m_mapToFixings[fdate];

		for ( ; idate < nextDate; idate++ )
		{
			double lvol = computeLocalVol( idate, rspot, plv );	
	
			double lvol_bp = computeLocalVol( idate, rspot_vup, pblv );

			double dW = *g++ ;
			double drift = *q++ ;
			double dt = *pdt++;

			rspot  +=  lvol * ( ( drift -0.5*lvol )* dt +  dW );
			rspot_vup  +=  lvol_bp * ( ( drift -0.5*lvol_bp )* dt +  dW );
		}

		double fwd = m_forwards[fdate] ;
		m_currentPathValues[3][fdate] = fwd * exp( rspot_vup ) ;

		double sim_spot = fwd * exp( rspot ) ;
		m_currentPathValues[0][fdate] = sim_spot;

		m_currentPathValues[1][fdate] = 1.01 * sim_spot ;
		m_currentPathValues[2][fdate] = 0.99 * sim_spot ;
	}

	return m_currentPathValues;
}
*/




inline
double CLocalVolMC::computeLocalVol(int time, double rspot, const double* lv_0 )
{
	double igrid = (rspot - m_min_yGrid) / m_dy;

	if( igrid > m_max_yIndex )	return *( lv_0 + time * m_space_size + int(m_max_yIndex) );	
	if( igrid < 1e-5)			return *( lv_0 + time * m_space_size );	

	int iigrid = int(igrid);
	double wght = igrid - iigrid;
	const double* lv = lv_0 + time * m_space_size + iigrid;

	double localvol =  (1.-wght) * (*lv++) + wght * (*lv);
	return localvol ;
}



void CLocalVolMC::generatePaths()	// should be used only for handle
{
	if ( m_storedPaths || m_nDates == 0 ){
		return;
	}
	m_path_array.resize(m_nPaths,m_nAssets,m_nMcDates); //
	m_pPath_array	= &m_path_array;
	m_iasset		= 0; // ??

	for(int ipath=0; ipath<m_nPaths; ipath++)
	{
		(*m_pPath_array)[ipath] = GetPathArray( ipath );
	}

	m_storedPaths = true;
}







CMultiLocalVolMC::CMultiLocalVolMC(MlEqConstDateHandle hDate, MlEqAssetHandle hUdly)
:CLocalVolMC(hDate, hUdly)
{
	m_assets = hUdly->GetAssets();
	m_nAssets = m_assets.size();

	m_lvmcComponent.resize(m_nAssets);
	std::string	payCurr = m_assets[0]->GetPayZeroCurve(true)->GetName() ;
	
	for(int iasset=0; iasset<m_nAssets; ++iasset)
	{
		m_lvmcComponent[iasset] = new CLocalVolMC(m_hDate, m_assets[iasset]);
	
		if ( iasset>0 && m_assets[iasset]->GetPayZeroCurve(true)->GetName() != payCurr ){
			throw("pay currencies must be identical for all assets");
		}
	}
}



CMultiLocalVolMC::CMultiLocalVolMC(MlEqConstDateHandle hDate, const std::vector<MlEqAssetHandle>& vhUdly)
:CLocalVolMC(hDate, vhUdly[0])
{
	m_assets = vhUdly;
	m_nAssets = m_assets.size();

	m_lvmcComponent.resize(m_nAssets);
	std::string	payCurr = m_assets[0]->GetPayZeroCurve(true)->GetName() ;

	for(int iasset=0; iasset<m_nAssets; ++iasset)
	{
		m_lvmcComponent[iasset] = new CLocalVolMC(m_hDate, m_assets[iasset]);
		
		if ( iasset>0 && m_assets[iasset]->GetPayZeroCurve(true)->GetName() != payCurr ){
			throw("pay currencies must be identical for all assets");
		}
	}
}


void CMultiLocalVolMC::initialize( const std::vector<long>& payoffDates, const CMatrix& Fixings, int npath, int rngFlag, int stepAYear )
{
	initTimeGrids( payoffDates, stepAYear );	// we need that for the rng...

	m_currentPathValues.resize(m_nAssets, m_nMcDates);
	initRandomGenerator(npath, rngFlag);
	
	CMatrix correl;
	getCorrelationMatrix(m_assets, correl);	
	m_cholesky = Cholesky(correl, m_nDates);
	m_qCholesky = Cholesky(correl, m_nBBDates);

	if( Fixings.rows() != m_nAssets ){
		throw "wrong fixings size in local vol mc";
	}

	CMatrix fixings( 1, Fixings.cols() );
	for(int iasset=0; iasset<m_nAssets; ++iasset)
	{		
		fixings[0] = Fixings[iasset]; 
		m_lvmcComponent[iasset]->initialize(payoffDates, fixings, 1/*m_nPaths*/, 0/*rng*/, stepAYear); 
	}
	
	m_discount = m_lvmcComponent[0]->GetDiscounts(0);	// do better	
}



CMatrix& CMultiLocalVolMC::GetPathArray(int ipath)
{
	if( m_storedPaths ){
		return (*m_pPath_array)[ipath];
	}

	m_randomGenerator->generateRandoms(m_randoms,m_numberOfFactorsPerAsset,m_cholesky,ipath,m_rndTemp);
	m_quasiGenerator-> generateRandoms(m_quasis, m_numberOfFactorsPerAsset,m_qCholesky,ipath,m_qRndTemp);


	for(int iasset=0; iasset<m_nAssets; ++iasset)
	{
		BrownianBridgePath( m_randoms[iasset], m_quasis[iasset] );
		m_currentPathValues[iasset] = ( m_lvmcComponent[iasset]->MakePathArray( m_randoms[iasset] ) ) [0] ; 
	}

	return m_currentPathValues;
}




double	CMultiLocalVolMC::GetForwardValue(int idate,int iasset)
{
	double res = m_lvmcComponent[iasset]->GetForwardValue(idate,0);
	return res;
}




