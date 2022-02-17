


#include "stdafx.h"
#include "mleqobjects.h"

#include "utility.h"
#include "mleqshortmodels.h"

#include "cMatrix.h"
#include <assert.h>
#include <math.h>
#include "utility.h"
//#include "montecarlo.h"
#include "SciSobol.h"
#include "MlEqNormal.h"
#include "pdeproducts.h"




/***************************************************************
**	Class  : calibrateHermiteOptions 
**	Routine: ObjectiveFcn
**	Returns: nothing
**	Action : 
****************************************************************/


void calibrateHermiteOptions::ObjectiveFcn(double* gVals,double* xVals)
{
	CVector result(m_fixedStrikes.getsize());

	m_currentHermitCoeff[0] = xVals[0];
	m_currentHermitCoeff[1] = xVals[1];
	m_currentHermitCoeff[2] = -xVals[1]+xVals[2];

	for ( int i = 3 ; i < m_currentHermitCoeff.getsize(); i++ )
	{
		m_currentHermitCoeff[i] = xVals[i];
	}

	int returnVolFlag = 0;
	HermiteOptionsNew(result);

//	HermiteOptionsNew(result,m_fixedStrikes,m_currentHermitCoeff,m_maturity,m_forward,returnVolFlag,m_ntimesteps,m_npaths,m_nspacialPoints,m_nstdev,m_randomGenerator);

	double res = 0.0;
	for ( int i = 0 ; i < result.getsize();i++ ){
		res += m_calibWeights[i]*pow((result[i]-m_calibPrices[i]),2.0);
//		res += fabs((result[i]-m_calibPrices[i]));
	}

	res /= result.getsize();
	res /= m_forward;

	*gVals = res;
}


/***************************************************************
**	Class  : calibrateHermiteOptions 
**	Routine: ConstraintFcn
**	Returns: nothing
**	Action : 
****************************************************************/

void  calibrateHermiteOptions::ConstraintFcn(double* contraintVals, double* xVals)
{

	if ( m_varSwapVol.getsize() == 0 ){
		LSGRGSolver::ConstraintFcn(contraintVals, xVals);
		return;
	}

	double cval = 0.0;
	double fac = 1.0;
	for ( int i = 0 ; i < m_currentHermitCoeff.getsize(); i++ )
	{
		cval += pow(m_currentHermitCoeff[i],2.0)/fac;
		fac  *= (double)(i+2);
	}

	cval = sqrt(cval);
//	cval = cval-m_varSwapVol[0];

	contraintVals[0] = cval;
	contraintVals[1] = -cval;

}


/***************************************************************
**	Class  : calibrateHermiteOptions 
**	Routine: calibrate
**	Returns: nothing
**	Action : 
****************************************************************/

void  calibrateHermiteOptions::calibrate(CVector& x)
{

	int maximizeFlag = 0;
	int returnCode;

	x = m_hermitCoeffGuess;
	solve( x, maximizeFlag, returnCode);

	x = 	m_currentHermitCoeff;
}


/***************************************************************
**	Class  : calibrateHermiteOptions
**	Routine: ~calibrateHermiteOptions
**	Returns: nothing
**	Action : 
****************************************************************/

calibrateHermiteOptions::~calibrateHermiteOptions()
{

	int nherm = m_hermitCoeffGuess.getsize();

	for ( int i = 0 ; i < m_nspacialPoints; i++ ){
		for ( int j = 0 ; j < m_npaths; j++ ){
			for ( int k = 0 ; k < nherm; k++ ){
					delete [] m_cashedVar[i][j][k] ;

			}
			delete [] m_cashedVar[i][j];
		}
		delete [] m_cashedVar[i];
	}




}


/***************************************************************
**	Class  : calibrateHermiteOptions 
**	Routine: calibrate
**	Returns: nothing
**	Action : 
****************************************************************/

void  calibrateHermiteOptions::initialize(GVector< MlEqStrikeHandle >& calibStrikes,
										  CVector& calibVols,double	forward,
										  MlEqDateHandle	Today,
										  long				maturityDate,
										  CVector&			varSwapVol,
										  int				ntimesteps,
										  int				npaths,
										  int				nspacialPoints,
										  int				nstdev,
										  CVector&			hermitCoeffGuess
										  )


{

	m_calibStrikes				=	calibStrikes;

	MlEqStrike fstrike;
	m_fixedStrikes.resize(m_calibStrikes.getsize());

	for ( int i = 0 ; i < m_calibStrikes.getsize(); i++ )
	{	
		MlEqStrike::convertStrikes(fstrike,*m_calibStrikes[i]);
		m_fixedStrikes[i] = fstrike.m_strike;
	}

	m_calibVols		=	calibVols;
	m_forward		=	forward;
	m_maturityDate	=	maturityDate;
	m_Today			=	Today;
	m_maturity		=	m_Today->GetYearFraction(maturityDate);

	m_calibPrices.resize(m_fixedStrikes.getsize());
	for ( int i = 0 ; i < m_calibPrices.getsize(); i++ ){
		m_calibPrices[i] = Bs(m_forward,calibVols[i],m_maturity,m_fixedStrikes[i],1.0,1.0);
	}

	m_calibWeights.resize(m_calibPrices.getsize());
	double volshift = 0.01;
	for ( int i = 0 ; i < m_calibPrices.getsize(); i++ ){
		m_calibWeights[i] = (Bs(m_forward,calibVols[i]+volshift,m_maturity,m_fixedStrikes[i],1.0,1.0)-m_calibPrices[i])/(m_forward*volshift*0.4*sqrt(m_maturity));
	}


	m_varSwapVol			= 	varSwapVol;
	m_ntimesteps			=	ntimesteps;
	m_npaths				=	npaths;
	m_nspacialPoints		=	nspacialPoints;
	m_nstdev				=   nstdev;
	m_hermitCoeffGuess		=	hermitCoeffGuess;
	m_currentHermitCoeff	=	m_hermitCoeffGuess;



	int ndim				= ntimesteps;	
	int	randomNumberFlag	= 2;
	int nsteps = m_ntimesteps;

	if ( randomNumberFlag == 0 )
	{
			long idum = -123;
			m_randomGenerator = new	randomGenerator(idum,1,nsteps);
		
	}
	else if ( randomNumberFlag == 1 )
	{
			int dimension	= nsteps;
			int skip		= 1000;
			int leap		= 1;

			m_randomGenerator = new Sobol(dimension,skip,leap,0,0);

	}
	else
	{
			CUnitcube* D = new CUnitcube(ndim);
					//	select alpha
					//	CNiederreiteralpha *alpha = new CNiederreiteralpha;

			CAlpha* alpha;
			if ( randomNumberFlag == 2 ){
						alpha = new CPrimerootalpha;
			}
			else if ( randomNumberFlag == 3 ){
						alpha = new CNiederreiteralpha;
			}
			else if ( randomNumberFlag == 4 ){
						alpha = new CBakeralpha;
			}
			else if ( randomNumberFlag == 5 ){
						alpha = new CPrimerootalpha;
			}
			else{
						throw("randomNumberFlag must be between 1 and 5");
			}

			CArithmeticmeanweights* weights = new CArithmeticmeanweights;
			int d = 1;// d >= r 
			CBakerperiodization* periodization=new CBakerperiodization;
			CNTparameters* par = new CNTparameters(alpha,weights,periodization);
			CNTintegrator* NTintegrator= new CNTintegrator;

			m_randomGenerator = new Diophantine(ndim,npaths,alpha,weights,D,periodization,NTintegrator,par);
	}



	m_gaussPoint.resize(nsteps);
	m_gaussWeight.resize(nsteps);

	MlEqMaths::dGauleg(0.0,1.0,m_gaussPoint,m_gaussWeight,nsteps,true);			

	m_spacialGaussPoint.resize(m_nspacialPoints);
	m_spacialGaussWeight.resize(m_nspacialPoints);
	MlEqMaths::dGauleg(-nstdev,+nstdev,m_spacialGaussPoint,m_spacialGaussWeight,m_nspacialPoints,true);			

// start cashing hermits

	double var,I,omega,val;
	CVector randoms(nsteps);
	int nherm = m_currentHermitCoeff.getsize();
	CVector rands(randoms.getsize());
	CVector res(m_fixedStrikes.getsize());

	double lastOmega,t,delta_t;

/*	m_cashedVar.resize(m_nspacialPoints,npaths);
	for ( int i = 0 ; i < m_nspacialPoints; i++ ){
		for ( int j = 0 ; j < npaths; j++ ){
			m_cashedVar[i][j].resize(nherm,nherm);}
	}
*/

//	double****					m_cashedVar;// don't do this at home

	m_cashedVar = new double***[m_nspacialPoints];
	for ( int i = 0 ; i < m_nspacialPoints; i++ ){
		m_cashedVar[i] = new double**[npaths];
		for ( int j = 0 ; j < npaths; j++ ){
			m_cashedVar[i][j] = new double*[nherm];
			for ( int k = 0 ; k < nherm; k++ ){
				m_cashedVar[i][j][k] = new double[nherm];
			}
		}
	}


	for ( int ipoint = 0 ; ipoint < m_nspacialPoints; ipoint++ )
	{
		lastOmega = m_spacialGaussPoint[ipoint];

		val = 0.0;
		for ( int istrike = 0 ; istrike < res.getsize(); istrike++ ){
			res[istrike] = 0.0;
		}

		for ( int ipath = 0 ; ipath < npaths; ipath++ )
		{
			m_randomGenerator->generateRandoms(randoms,ipath);

			for ( int i = 0 ; i < rands.getsize(); i++ )
			{
				if ( i ){
					delta_t = (m_gaussPoint[i]-m_gaussPoint[i-1]);
				}
				else{
					delta_t = m_gaussPoint[i];
				}

				rands[i] = randoms[i]*sqrt(delta_t);
			}

			for ( int i = 1 ; i < rands.getsize(); i++ ){
				rands[i] += rands[i-1];
			}


			var = 0.0;
			for ( int i = 0 ; i < nherm; i++ )
			{
				for ( int j = 0 ; j < nherm; j++ )
				{
					I = 0.0;
					for ( int istep = 0 ; istep < nsteps; istep++ )
					{
						t = m_gaussPoint[istep];
						omega = lastOmega*t+rands[istep]-t*rands[nsteps-1];

						I += m_gaussWeight[istep]*
							 Hermit(i,omega,t)*Hermit(j,omega,t);
					}
					double z = (double)(i+1)*(double)(j+1)*I*pow(m_maturity,(double)(i+j)/2.0+1.0);
					m_cashedVar[ipoint][ipath][i][j] = z;
				}
			}


		}
	}




	CVector initialGuess;

	initialGuess = m_hermitCoeffGuess;

		
	iVector linearVariableIndex;
	CMatrix xValBounds;
	CMatrix ObjectiveBounds;
	CMatrix NonZeroPartialDerivativesSpecification;


	double InitialTolerance		= 1e-3;
	double FinalTolerance		= 1e-3;
	double StoppingTolerance	= 1e-3;

	int outputFlag = 1;

	xValBounds.resize(initialGuess.getsize(),2);
	xValBounds[0][0] = -0.0;
	xValBounds[0][1] =  0.6;

	for ( int i = 1 ; i < xValBounds.rows(); i++ )
	{
		xValBounds[i][0] = -0.2;
		xValBounds[i][1] = 0.2;
	}

	if ( m_varSwapVol.getsize() > 0 )
	{
		ObjectiveBounds.resize(1,2);
		ObjectiveBounds[0][0] =  m_varSwapVol[0]-0.003;
		ObjectiveBounds[0][1] =  m_varSwapVol[0]+0.003;
	}

	LSGRGSolver::initialize(   initialGuess,linearVariableIndex,
							   xValBounds,ObjectiveBounds,
							   InitialTolerance,FinalTolerance,
							   StoppingTolerance,
							   NonZeroPartialDerivativesSpecification,outputFlag);


}


/***************************************************************
**	Class  : calibrateHermiteOptions 
**	Routine: calibrate
**	Returns: nothing
**	Action : 
****************************************************************/



										
void calibrateHermiteOptions::HermiteOptionsNew(CVector& result)
{

	if ( m_maturity < 1e-3 ){
		throw("maturity zero entered");
	}

	for ( int istrike = 0 ; istrike < m_fixedStrikes.getsize(); istrike++ ){
				result[istrike] = 0.0;
	}			

	int nsteps = m_ntimesteps;
	double deltaT = 1.0/(double)nsteps;

	CVector coeff(m_currentHermitCoeff.getsize());

	double n = 1.0;
	for (int  i = 0 ; i < coeff.getsize(); i++ )
	{
		coeff[i] = m_currentHermitCoeff[i]/(n*pow(m_maturity,(double)i/2.0));
		n*= (double)(i+2);
	}


//  estimate at-the money vol



	double stateVar,var;
	CVector randoms(nsteps);
	int nherm = coeff.getsize();



	double sqrt_2pi =  sqrt(2.0*3.141592654);

	double value = 0.0,val;

	CVector rands(randoms.getsize());
	CVector res(result.getsize());

	double lastOmega,density,spot;

	for ( int ipoint = 0 ; ipoint < m_nspacialPoints; ipoint++ )
	{
		lastOmega = m_spacialGaussPoint[ipoint];

		stateVar = 0.0;
		for ( int i = 0 ; i < nherm; i++ ){
			stateVar  += coeff[i]*Hermit(i+1,lastOmega,1.0)*pow(m_maturity,(double)(i+1)/2.0);
		}

		val = 0.0;
		for ( int istrike = 0 ; istrike < res.getsize(); istrike++ ){
			res[istrike] = 0.0;
		}

		for ( int ipath = 0 ; ipath < m_npaths; ipath++ )
		{
			var = 0.0;
			for ( int i = 0 ; i < nherm; i++ )
			{
				for ( int j = 0 ; j < i; j++ )
				{
					var += coeff[i]*coeff[j]*m_cashedVar[ipoint][ipath][i][j];
				}
			}
			var *= 2;// these were the offdiagonal terms

			// add diagonal terms
			for ( int i = 0 ; i < nherm; i++ )
			{
				int j = i;
				var += coeff[i]*coeff[j]*m_cashedVar[ipoint][ipath][i][j];
			}

			spot = exp(stateVar-0.5*var);
			val  += spot;
			spot *= m_forward;

			for ( int istrike = 0 ; istrike < m_fixedStrikes.getsize(); istrike++ ){
				res[istrike] += MlEqMaths::Max(spot-m_fixedStrikes[istrike],0.0);
			}			

		}


		density = m_spacialGaussWeight[ipoint]*exp(-0.5*lastOmega*lastOmega)/sqrt_2pi; 

		for ( int istrike = 0 ; istrike < m_fixedStrikes.getsize(); istrike++ ){
				res[istrike] *= density;
			}			


		val *= density;
		value += val;

		for ( int istrike = 0 ; istrike < m_fixedStrikes.getsize(); istrike++ ){
				result[istrike] += res[istrike];
			}			

	}


	value /= m_npaths;

	for ( int istrike = 0 ; istrike < m_fixedStrikes.getsize(); istrike++ ){
				result[istrike] /= m_npaths;
		}			


/*	if ( returnVolFlag )
	{
		double accuracy   = 1e-10;
		double lower_x    = 0.03; 
		double upper_x    = 1.5; 

		int rootfind_flag = eROOTFIND_BRENT_NOGROW;

		CVector impliedVols(strikes.getsize());
		for (int i = 0; i < strikes.getsize(); i++ )
		{

			impliedVols[i] =	MlEqBSImpliedVol(
								result[i],
								forward,
								maturity,
								strikes[i],
								1,1,rootfind_flag ,	 
								accuracy   ,		 
								lower_x    ,		 
								upper_x);
			
		}

		result = impliedVols;

	}




*/

}


/***************************************************************
**	Class  : none 
**	Routine: impliedLogContractVol
**	Returns: double
**	Action : 
****************************************************************/

	
double impliedLogContractVol(MlEqAsset& asset,long nMaturity,int ngauss,double upper,double lower)
{	
	
	CVector gaussPoint(ngauss),gaussWeight(ngauss);

	MlEqMaths::dGauleg(lower,upper,gaussPoint,gaussWeight,ngauss,true);			

	double fwd = asset.GetForward(nMaturity,false);
	MlEqStrike stk(fwd);
	double vol,st;

	double volatm = asset.GetVolatility(stk,nMaturity);

	MlEqConstDateHandle	date = asset.GetDateHandle();
	double mat = date->GetYearFraction(nMaturity);

	double res = 0.0;
	for ( int i = 0 ; i < ngauss; i++ )
	{
		MlEqNormalizedStrike normStrike;
		normStrike.initialize(asset,gaussPoint[i], nMaturity,true,date);
		vol = asset.GetVolatility(normStrike,nMaturity);

		MlEqStrike strike;
		MlEqStrike::convertStrikes(strike,normStrike);

		st = strike.m_strike/fwd;

		if ( normStrike.m_strike < 0 ){
			res += gaussWeight[i]/pow(st,1.0)*Bs(1,vol,mat,st,1,-1);
		}
		else
		{
			res += gaussWeight[i]/pow(st,1.0)*Bs(1,vol,mat,st,1,1);
		}
	}

	res *= volatm*sqrt(mat);
	res *= 2.0/mat;
	res = sqrt(res);

	return res;
}




class expectedVolInfo
{
	public:

	double mat;
	double fwd;
	double eps;
	long nMaturity;

	MlEqAssetHandle pAsset;

};

/***************************************************************
**	Class  : none 
**	Routine: expectedVol
**	Returns: nothing
**	Action : 
****************************************************************/


void expectedVol(double x, double z[], void *ptr)
{	

	expectedVolInfo* p = (expectedVolInfo * ) ptr;

	double	mat			= p->mat;
	double	eps			= p->eps;
	double	fwd			= p->fwd;
	long nMaturity      = p->nMaturity;

	MlEqAssetHandle pAsset = p->pAsset;

	double vol_up,vol_do,vol;

	double fac = 1.0;//0.85;
	x*= fac;

	MlEqStrike strike(x);
	vol = pAsset->GetVolatility(strike,nMaturity);

	MlEqStrike strike_up(x*(1.0+eps));
	vol_up = pAsset->GetVolatility(strike_up,nMaturity);

	MlEqStrike strike_do(x*(1.0-eps));
	vol_do = pAsset->GetVolatility(strike_do,nMaturity);



	double density = 

	( Bs(fwd,vol_up,mat,x*(1.0+eps),1.0,1) +
	  Bs(fwd,vol_do,mat,x*(1.0-eps),1.0,1) -
	  2.0*Bs(fwd,vol,mat,x,1.0,1)
	  );



	density /= pow(eps*x,2.0);

	z[1] = density*vol*vol;//*MlEqMaths::Max(x-fwd,0.0);//vol*vol;


}	

static int counter;

/***************************************************************
**	Class  : none 
**	Routine: expectedVol
**	Returns: double
**	Action : 
****************************************************************/


double expectedVol(MlEqAssetHandle& asset,long nMaturity,int ngauss,double upper,double lower)
{	
	
	double fwd = asset->GetForward(nMaturity,false);
	MlEqStrike stk(fwd);

	double volatm = asset->GetVolatility(stk,nMaturity);

	MlEqConstDateHandle	date = asset->GetDateHandle();
	double mat = date->GetYearFraction(nMaturity);


	double eps=1e-4,hInit=0.1,hMin=1e-2,bump=0.0008,scale = 10;
	int nfunctions = 1;

	CVector integrationLimits(2);
	integrationLimits[0] = fwd*exp(lower*volatm*sqrt(mat));
	integrationLimits[1] = fwd*exp( upper*volatm*sqrt(mat));

	CVector results(2);
	counter = 0;

	expectedVolInfo info;
	info.fwd = fwd;
	info.mat = mat;
	info.pAsset = asset;
	info.eps = 0.001;
	info.nMaturity = nMaturity;

	double integral;
	Runge_Kutta_Integrate_new(integral,expectedVol,integrationLimits,hMin,eps,bump,&info, hInit ,scale);

	integral = sqrt(integral);

	return integral;
}


double PatHaagenImpliedVol(double strike,double mat,double fwd,double correl,double beta,double volvol,double shortvol)
{
	if( fabs( volvol ) < 1e-6 )
	{
		double f_b = pow(0.5*(fwd+strike), 1.-beta) ;

		double detvol = pow( (1.-beta)*shortvol/f_b, 2 ) / 24. * mat;
		detvol += (1.-beta)*(2.+beta) * pow( 2.*(fwd-strike)/(fwd+strike) , 2 ) / 24. ;
		detvol = shortvol/f_b * ( 1. + detvol ) ;
		return detvol;
	}


	double alpha = shortvol;


	if( fabs( strike-fwd ) < 1e-6 )		// atm case
	{
		double f_b = pow( fwd, 1.-beta) ;
		double volatm = pow(alpha*(1.-beta)/f_b, 2) / 24. 
					  + 0.25 * correl*beta*alpha*volvol / f_b
					  + ( 2. - 3.*correl*correl) * volvol*volvol / 24. ;

		volatm = alpha/f_b * ( 1. + volatm * mat );

		return volatm;
	}

	double s = fwd*strike;
	double ls = log(fwd/strike);

	double z = volvol/alpha*pow(s,(1.0-beta)/2.0)*ls;

	double x_z = ( sqrt(1.0-2.0*correl*z+z*z)+z-correl)/(1.0-correl);
		   x_z = log(x_z);


	double vol;

	vol = pow(s,(1.0-beta)/2.0)*(1.0+pow((1.0-beta)*ls,2.0)/24.0+pow((1.0-beta)*ls,4.0)/1920.0);
	vol = 1.0/vol;
	vol *= alpha*z/x_z;

	vol *= (1.0 + 

			( pow((1.0-beta)*alpha,2.0)/(24.0*pow(s,1.0-beta))+0.25*correl*beta*volvol*alpha/pow(s,(1.0-beta)/2.0)+
			  (2.0-3.0*correl*correl)/24.0*pow(volvol,2.0) 
		   )*mat
		   );

	return vol;

}

// simplified formula if beta = 1
double PatHaagenImpliedVol(double strike,double mat,double fwd,double correl,double volvol,double shortvol)
{
	if( fabs(volvol) < 1e-6 )	return shortvol;

	if( fabs( strike-fwd ) < 1e-6 )		// atm case
	{
		double volatm =  0.25 * correl*shortvol*volvol 
					  + ( 2.-3.*correl*correl) * volvol*volvol / 24. ;

		volatm = shortvol * ( 1. + volatm * mat );
		return volatm;
	}

	double ls = log(fwd/strike);
	double z = volvol/shortvol*ls;

	double x_z = ( sqrt(1.-2.*correl*z+z*z) + z-correl ) / (1.-correl);
	x_z = log(x_z);

	double vol = shortvol* z/x_z;

	vol *=  1. + ( 0.25*correl*volvol*shortvol + (2.-3.*correl*correl)/24.*volvol*volvol ) * mat ;
	return vol;
}


double BlackScholesSwaptionPricer(long valuationDate, GVector<long>& SwapFixingDates,double swaptionStrike,double vol,DayCountConventionEnum & dayCount ,MlEqZeroCurve&	Z, int callPut=1)
{
//  callPut = 1 : payer Swaption

	MlEqDate DateToDouble(valuationDate,dayCount);

	double FSR,res;
	double DV01 = 0.0;
	double delta;
	for ( int i = 1 ; i < SwapFixingDates.getsize(); i++ )
	{
		delta  = DateToDouble.GetYearFraction(SwapFixingDates[i])-DateToDouble.GetYearFraction(SwapFixingDates[i-1]);
		DV01 += delta*Z.GetDiscountFactor(SwapFixingDates[i]);
	}

	FSR = (Z.GetDiscountFactor(SwapFixingDates[0])-Z.GetDiscountFactor(SwapFixingDates[SwapFixingDates.getsize()-1]))/DV01;
	double mat = DateToDouble.GetYearFraction(SwapFixingDates[0]);
	res =  DV01*Bs(FSR,vol,mat,swaptionStrike,1.0,callPut);
	return res;
}	
	
	
void HermiteOptionsNew(CVector& result, CVector& strikes,const  CVector& hermitCoeff,double maturity,double forward,int returnVolFlag,int ntimesteps,int npaths,int nspacialPoints,int nstdev)
{
	randomGeneratorHandle RandomGenerator;
	HermiteOptionsNew(result, strikes,hermitCoeff,maturity,forward,returnVolFlag,ntimesteps,npaths,nspacialPoints,nstdev,RandomGenerator);
}



void HermiteOptionsNew(CVector& result, CVector& strikes,const  CVector& hermitCoeff,double maturity,double forward,int returnVolFlag,int ntimesteps,int npaths,int nspacialPoints,int nstdev,	randomGeneratorHandle& RandomGenerator)

{

	if ( maturity < 1e-3 ){
		throw("maturity zero entered");
	}

	for ( int istrike = 0 ; istrike < strikes.getsize(); istrike++ ){
				result[istrike] = 0.0;
	}			

	int nsteps = ntimesteps;
	double deltaT = 1.0/(double)nsteps;

	CVector coeff(hermitCoeff.getsize());

	double n = 1.0;
	for (int  i = 0 ; i < coeff.getsize(); i++ )
	{
		coeff[i] = hermitCoeff[i]/(n*pow(maturity,(double)i/2.0));
		n*= (double)(i+2);
	}


//  estimate at-the money vol

	int ndim				= nsteps;	
	int	randomNumberFlag	= 2;

;

	if ( !RandomGenerator)
	{
		if ( randomNumberFlag == 0 )
		{
			long idum = -123;
			RandomGenerator = new	randomGenerator(idum,1,nsteps);
		
		}
		else if ( randomNumberFlag == 1 )
		{
			int dimension	= nsteps;
			int skip		= 1000;
			int leap		= 1;

			RandomGenerator = new Sobol(dimension,skip,leap,0,0);

		}
		else
		{
			CUnitcube* D = new CUnitcube(ndim);
					//	select alpha
					//	CNiederreiteralpha *alpha = new CNiederreiteralpha;

			CAlpha* alpha;
			if ( randomNumberFlag == 2 ){
						alpha = new CPrimerootalpha;
			}
			else if ( randomNumberFlag == 3 ){
						alpha = new CNiederreiteralpha;
			}
			else if ( randomNumberFlag == 4 ){
						alpha = new CBakeralpha;
			}
			else if ( randomNumberFlag == 5 ){
						alpha = new CPrimerootalpha;
			}
			else{
						throw("randomNumberFlag must be between 1 and 5");
			}

			CArithmeticmeanweights* weights = new CArithmeticmeanweights;
			int d = 1;// d >= r 
			CBakerperiodization* periodization=new CBakerperiodization;
			CNTparameters* par = new CNTparameters(alpha,weights,periodization);
			CNTintegrator* NTintegrator= new CNTintegrator;

			RandomGenerator = new Diophantine(ndim,npaths,alpha,weights,D,periodization,NTintegrator,par);

		}

	}

	double stateVar,var,I,omega;
	CVector randoms(nsteps);
	int nherm = coeff.getsize();

	CVector gaussPoint(nsteps),gaussWeight(nsteps);
	MlEqMaths::dGauleg(0.0,1.0,gaussPoint,gaussWeight,nsteps,true);			

	CVector spacialGaussPoint(nspacialPoints),spacialGaussWeight(nspacialPoints);
	MlEqMaths::dGauleg(-nstdev,+nstdev,spacialGaussPoint,spacialGaussWeight,nspacialPoints,true);			

	double sqrt_2pi =  sqrt(2.0*3.141592654);

	double value = 0.0,val;

	CVector rands(randoms.getsize());
	CVector res(result.getsize());

	double lastOmega,t,delta_t,spot,density;

	for ( int ipoint = 0 ; ipoint < nspacialPoints; ipoint++ )
	{
		lastOmega = spacialGaussPoint[ipoint];

		stateVar = 0.0;
		for ( int i = 0 ; i < nherm; i++ ){
			stateVar  += coeff[i]*Hermit(i+1,lastOmega,1.0)*pow(maturity,(double)(i+1)/2.0);
		}

		val = 0.0;
		for ( int istrike = 0 ; istrike < res.getsize(); istrike++ ){
			res[istrike] = 0.0;
		}

		for ( int ipath = 0 ; ipath < npaths; ipath++ )
		{
			RandomGenerator->generateRandoms(randoms,ipath);

			for ( int i = 0 ; i < rands.getsize(); i++ )
			{
				if ( i ){
					delta_t = (gaussPoint[i]-gaussPoint[i-1]);
				}
				else{
					delta_t = gaussPoint[i];
				}

				rands[i] = randoms[i]*sqrt(delta_t);
			}

			for ( int i = 1 ; i < rands.getsize(); i++ ){
				rands[i] += rands[i-1];
			}


			var = 0.0;
			for ( int i = 0 ; i < nherm; i++ )
			{
				for ( int j = 0 ; j < nherm; j++ )
				{
					I = 0.0;
					for ( int istep = 0 ; istep < nsteps; istep++ )
					{
						t = gaussPoint[istep];
						omega = lastOmega*t+rands[istep]-t*rands[nsteps-1];

						I += gaussWeight[istep]*
							 Hermit(i,omega,t)*Hermit(j,omega,t);
					}
					var += coeff[i]*coeff[j]*(double)(i+1)*(double)(j+1)*I*pow(maturity,(double)(i+j)/2.0+1.0);
				}
			}

			spot = exp(stateVar-0.5*var);
			val  += spot;
			spot *= forward;

			for ( int istrike = 0 ; istrike < strikes.getsize(); istrike++ ){
				res[istrike] += MlEqMaths::Max(spot-strikes[istrike],0.0);
			}			

		}


		density = spacialGaussWeight[ipoint]*exp(-0.5*lastOmega*lastOmega)/sqrt_2pi; 

		for ( int istrike = 0 ; istrike < strikes.getsize(); istrike++ ){
				res[istrike] *= density;
			}			


		val *= density;
		value += val;

		for ( int istrike = 0 ; istrike < strikes.getsize(); istrike++ ){
				result[istrike] += res[istrike];
			}			

	}


	value /= npaths;

	for ( int istrike = 0 ; istrike < strikes.getsize(); istrike++ ){
				result[istrike] /= npaths;
		}			


	if ( returnVolFlag )
	{
		double accuracy   = 1e-10;
		double lower_x    = 0.03; 
		double upper_x    = 1.5; 

		int rootfind_flag = eROOTFIND_BRENT_NOGROW;

		CVector impliedVols(strikes.getsize());
		for (int i = 0; i < strikes.getsize(); i++ )
		{

			impliedVols[i] =	MlEqBSImpliedVol(
								result[i],
								forward,
								maturity,
								strikes[i],
								1,1,rootfind_flag ,	 
								accuracy   ,		 
								lower_x    ,		 
								upper_x);
			
		}

		result = impliedVols;

	}

}




void VolOptions(CMatrix& result, CVector& strikes,const  CVector& hermitCoeff,double maturity,double forward,int returnVolFlag,int ntimesteps,int npaths,int nspacialPoints,int nstdev)
{

	if ( maturity < 1e-3 ){
		throw("maturity zero entered");
	}

	result.resize(strikes.getsize()+1,2);

	int nsteps = ntimesteps;
	double deltaT = 1.0/(double)nsteps;

	CVector coeff(hermitCoeff.getsize());

	double n = 1.0;
	for (int  i = 0 ; i < coeff.getsize(); i++ )
	{
		coeff[i] = hermitCoeff[i]/(n*pow(maturity,(double)i/2.0));
		n*= (double)(i+2);
	}


//  estimate at-the money vol

	int ndim				= nsteps;	
	int	randomNumberFlag	= 2;
	randomGeneratorHandle RandomGenerator;


	if ( randomNumberFlag == 0 )
	{
		long idum = -123;
		RandomGenerator = new	randomGenerator(idum,1,nsteps);
	
	}
	else if ( randomNumberFlag == 1 )
	{
		int dimension	= nsteps;
		int skip		= 1000;
		int leap		= 1;

		RandomGenerator = new Sobol(dimension,skip,leap,0,0);

	}
	else
	{
		CUnitcube* D = new CUnitcube(ndim);
				//	select alpha
				//	CNiederreiteralpha *alpha = new CNiederreiteralpha;

		CAlpha* alpha;
		if ( randomNumberFlag == 2 ){
					alpha = new CPrimerootalpha;
		}
		else if ( randomNumberFlag == 3 ){
					alpha = new CNiederreiteralpha;
		}
		else if ( randomNumberFlag == 4 ){
					alpha = new CBakeralpha;
		}
		else if ( randomNumberFlag == 5 ){
					alpha = new CPrimerootalpha;
		}
		else{
					throw("randomNumberFlag must be between 1 and 5");
		}

		CArithmeticmeanweights* weights = new CArithmeticmeanweights;
		int d = 1;// d >= r 
		CBakerperiodization* periodization=new CBakerperiodization;
		CNTparameters* par = new CNTparameters(alpha,weights,periodization);
		CNTintegrator* NTintegrator= new CNTintegrator;

		RandomGenerator = new Diophantine(ndim,npaths,alpha,weights,D,periodization,NTintegrator,par);

	}



	double stateVar,var,I,omega;
	CVector randoms(nsteps);
	int nherm = coeff.getsize();

	CVector gaussPoint(nsteps),gaussWeight(nsteps);
	MlEqMaths::dGauleg(0.0,1.0,gaussPoint,gaussWeight,nsteps,true);			

	CVector spacialGaussPoint(nspacialPoints),spacialGaussWeight(nspacialPoints);
	MlEqMaths::dGauleg(-nstdev,+nstdev,spacialGaussPoint,spacialGaussWeight,nspacialPoints,true);			

	double sqrt_2pi =  sqrt(2.0*3.141592654);

	double value = 0.0,val;

	CVector rands(randoms.getsize());
	CVector res(result.rows());

	double lastOmega,t,delta_t,density;

	double fwd=0,xfwd,fwdValue=0.0;

	for ( int ipoint = 0 ; ipoint < nspacialPoints; ipoint++ )
	{
		lastOmega = spacialGaussPoint[ipoint];

		stateVar = 0.0;
		for ( int i = 0 ; i < nherm; i++ ){
			stateVar  += coeff[i]*Hermit(i+1,lastOmega,1.0)*pow(maturity,(double)(i+1)/2.0);
		}

		val = 0.0;
		xfwd = 0.0;
	
		for ( int istrike = 0 ; istrike < res.getsize(); istrike++ ){
			res[istrike] = 0.0;
		}

		for ( int ipath = 0 ; ipath < npaths; ipath++ )
		{
			RandomGenerator->generateRandoms(randoms,ipath);

			for ( int i = 0 ; i < rands.getsize(); i++ )
			{
				if ( i ){
					delta_t = (gaussPoint[i]-gaussPoint[i-1]);
				}
				else{
					delta_t = gaussPoint[i];
				}

				rands[i] = randoms[i]*sqrt(delta_t);
			}

			for ( int i = 1 ; i < rands.getsize(); i++ ){
				rands[i] += rands[i-1];
			}

			var = 0.0;
			for ( int i = 0 ; i < nherm; i++ )
			{
				for ( int j = 0 ; j < nherm; j++ )
				{
					I = 0.0;
					for ( int istep = 0 ; istep < nsteps; istep++ )
					{
						t = gaussPoint[istep];
						omega = lastOmega*t+rands[istep]-t*rands[nsteps-1];

						I += gaussWeight[istep]*
							 Hermit(i,omega,t)*Hermit(j,omega,t);
					}
					var += coeff[i]*coeff[j]*(double)(i+1)*(double)(j+1)*I*pow(maturity,(double)(i+j)/2.0+1.0);
				}
			}
			
						
//			val		+= spot;
			xfwd	+= var;
			
			for ( int istrike = 0 ; istrike < strikes.getsize(); istrike++ ){
				res[istrike] += MlEqMaths::Max(var-strikes[istrike],0.0);
			}			
		
			if ( var < 0.0 ){
				int iii = 0.0;
			}
		}	
			
			
		density = spacialGaussWeight[ipoint]*exp(-0.5*lastOmega*lastOmega)/sqrt_2pi; 
		
		for ( int istrike = 0 ; istrike < strikes.getsize(); istrike++ ){
				res[istrike] *= density;
			}			
		
		
		val  *= density;
		xfwd *= density;

		fwd += xfwd;
		value += val;

		for ( int istrike = 0 ; istrike < strikes.getsize(); istrike++ ){
				result[istrike][0] += res[istrike];
			}			
	}


	value	/= npaths;
	fwd		/= npaths;

	result[strikes.getsize()][1] = fwd;

	for ( int istrike = 0 ; istrike < strikes.getsize(); istrike++ ){
				result[istrike][0] /= npaths;
		}			


	double accuracy   = 1e-16;
	double lower_x    = 0.01; 
	double upper_x    = 1.5; 

	int rootfind_flag = eROOTFIND_BRENT_NOGROW;

	CVector impliedVols(strikes.getsize());
	for (int i = 0; i < strikes.getsize(); i++ )
	{
		double x;
		try {		

			x =	MlEqBSImpliedVol(
								result[i][0],
								fwd,
								maturity,
								strikes[i],
								1,1,rootfind_flag ,	 
								accuracy   ,		 
								lower_x    ,		 
								upper_x);			

		} catch(...) {

			x = -10000.0;
		}

		result[i][1] = x;
	}	
}








void AmericanDownInPut(double& price,CMatrix& AmPrices,int method,double DownBarrier,
					   double spot,double mat,double r,double q,
					   double SpotBarrierVols,double SpotBarrierVolSlope,
					   double FwdBarrierVols,double FwdBarrierVolsSlope,
					   bool PayAtTheEnd,
					   int npoints,int nplotPoints,bool plotCummulaticeHitProb)
{


	

	if ( method == 3 )
	{

		if ( nplotPoints ){
			AmPrices.resize(nplotPoints,2);
		}

		BootstrapStoppingTimes(  price,AmPrices,					 
								  DownBarrier, spot, mat,
								  SpotBarrierVols,r-q, r,
								  SpotBarrierVolSlope,
								  FwdBarrierVols,
								  FwdBarrierVolsSlope,
								  PayAtTheEnd,
								  npoints,plotCummulaticeHitProb);
	}
	else{
		throw("only method 3 is implemented");
	}

}

	   

void BootstrapStoppingTimes(double& price,CMatrix& HitProb,					 
					 double Barrier,double spot,double mat,
					 double SpotBarrierVols,double drift,double r,
					 double SpotBarrierVolSlope,
					 double FwdBarrierVols,
					 double FwdBarrierVolsSlope,bool payAtTheEnd,
					 int npoints,bool plotCummulativeHitProb)
{


	int n;

	double ud = (spot < Barrier)?1.:-1. ;

	CVector timeToMat(npoints);
	double delta_t = mat/(double)npoints;
	
	timeToMat[0] = mat;
	for ( n = 1 ; n < npoints; n++ ){
		timeToMat[n] = timeToMat[n-1] - delta_t;
	}

//  determine stopping densities

	CVector stoppingDensity(npoints);
	double xDigital, integral,Digital,fwd,t,deltaVol,eps=1e-5;

	for ( n = 1 ; n < npoints; n++ )
	{

		double tt1 = timeToMat[n];
		integral=0.0;

		for ( int m = 1 ; m < n; m++ )
		{
			double tt2 = timeToMat[m];

			t = tt2-tt1;
			deltaVol = FwdBarrierVolsSlope*eps;

			fwd = Barrier*exp(drift*t);

			xDigital = 

			( Bs(fwd,FwdBarrierVols,t,Barrier,1.0,ud) - 
			  Bs(fwd,FwdBarrierVols+ud*deltaVol,t, Barrier*(1.0+ud*eps),1.0,ud)
			)/(eps*Barrier);

			if ( xDigital < 0.0 )
			{
				if( xDigital < -1e-10 ){
					throw("negative value for digital encountered");
				}
				else{
					xDigital = 0.0;
				}
			}

			integral += stoppingDensity[m]*xDigital ;
		}

		deltaVol = SpotBarrierVolSlope*eps*Barrier;


		t = mat-tt1;
		fwd = spot*exp(drift*t);


		Digital = 
			( Bs(fwd,SpotBarrierVols,t, Barrier,1.0,ud) - 
			  Bs(fwd,SpotBarrierVols+ud*deltaVol,t, Barrier*(1.0+ud*eps),1.0,ud)
			)/(eps*Barrier);


		if ( Digital < 0.0 )
		{
			if( Digital < -1e-10 ){
				throw("negative value for digital encountered");
			}
			else{
				Digital = 0.0;
			}
		}


		t = timeToMat[n-1]-timeToMat[n];
		fwd = Barrier*exp(drift*t);

		deltaVol = FwdBarrierVolsSlope*eps;

		xDigital = 
					( Bs(fwd,FwdBarrierVols,t, Barrier,1.0,ud) - 
					  Bs(fwd,FwdBarrierVols+ud*deltaVol,t, Barrier*(1.0+ud*eps),1.0,ud) )
					/(eps*Barrier);

		if ( xDigital < 0.0 )
		{
			if( xDigital < -1e-10 ){
				throw("negative value for digital encountered");
			}
			else{
				xDigital = 0.0;
			}
		}

		double windowProb = (Digital-integral)/xDigital;

		if ( windowProb < 0.0 )
		{
			if ( windowProb > -1e-5 ){
				windowProb = 0.0;
			}
			else{
				throw("negative hitting probability encountered");
			}
		}


		stoppingDensity[n] = windowProb;
	}

// calculate American digital now



	double df;
	price = 0.0;
	for ( n = 1 ; n < npoints; n++ )
	{
		t = mat-0.5*(timeToMat[n]+timeToMat[n-1]);
		df = exp(-r*t);

		if ( !payAtTheEnd ){
			price += stoppingDensity[n]*df;
		}else{
			price += stoppingDensity[n];
		}

	}

	if ( payAtTheEnd ){
		price *= exp(r*mat);
	}


	if ( HitProb.rows() > 0 )
	{

		double deltaT = mat/(double)HitProb.rows();

		int l = 0,k;
		for ( n = 0 ; n < HitProb.rows(); n++ )
		{
			t = (n+1)*deltaT;
			HitProb[n][0] = t;
			for ( k = l ; k < npoints; k++ )
			{
				if ( double kt = (k+1)*delta_t > t ){
					l = k;
					break;
				}

				HitProb[n][1] += stoppingDensity[k];
			}

		}

		if ( plotCummulativeHitProb )
		{
			for ( n = 0 ; n < HitProb.rows(); n++ )
			{
				if (n) {
					HitProb[n][1] += HitProb[n-1][1];
				}
			}
		}
	}


}

/////

void BootstrapStoppingTimes2(double& price,CMatrix& HitProb,					 
							 double Barrier,double spot,double mat,
							 double SpotBarrierVols,double drift,double r,
							 double SpotBarrierVolSlope,
							 double FwdBarrierVols,
							 double FwdBarrierVolsSlope,
							 double FwdBarrierStrikeVols,
							 double strike,
							 double cp,
							 double rebate,
							 bool payAtTheEnd,
							 int npoints,bool plotCummulativeHitProb)
{


	int n;

	double ud = (spot < Barrier)?1.:-1. ;

	CVector timeToMat(npoints);
	double delta_t = mat/(double)npoints;
	
	timeToMat[0] = mat;
	for ( n = 1 ; n < npoints; n++ ){
		timeToMat[n] = timeToMat[n-1] - delta_t;
	}

//  determine stopping densities

	CVector stoppingDensity(npoints);
	double xDigital, integral,Digital,fwd,t,deltaVol,eps=1e-5;

	for ( n = 1 ; n < npoints; n++ )
	{

		double tt1 = timeToMat[n];
		integral=0.0;

		for ( int m = 1 ; m < n; m++ )
		{
			double tt2 = timeToMat[m];

			t = tt2-tt1;
			deltaVol = FwdBarrierVolsSlope*eps;

			fwd = Barrier*exp(drift*t);

			xDigital = 

			( Bs(fwd,FwdBarrierVols,t,Barrier,1.0,ud) - 
			  Bs(fwd,FwdBarrierVols+ud*deltaVol,t, Barrier*(1.0+ud*eps),1.0,ud)
			)/(eps*Barrier);

			if ( xDigital < 0.0 )
			{
				if( xDigital < -1e-10 ){
					throw("negative value for digital encountered");
				}
				else{
					xDigital = 0.0;
				}
			}

			integral += stoppingDensity[m]*xDigital ;
		}

		deltaVol = SpotBarrierVolSlope*eps*Barrier;


		t = mat-tt1;
		fwd = spot*exp(drift*t);


		Digital = 
			( Bs(fwd,SpotBarrierVols,t, Barrier,1.0,ud) - 
			  Bs(fwd,SpotBarrierVols+ud*deltaVol,t, Barrier*(1.0+ud*eps),1.0,ud)
			)/(eps*Barrier);


		if ( Digital < 0.0 )
		{
			if( Digital < -1e-10 ){
				throw("negative value for digital encountered");
			}
			else{
				Digital = 0.0;
			}
		}


		t = timeToMat[n-1]-timeToMat[n];
		fwd = Barrier*exp(drift*t);

		deltaVol = FwdBarrierVolsSlope*eps;

		xDigital = 
					( Bs(fwd,FwdBarrierVols,t, Barrier,1.0,ud) - 
					  Bs(fwd,FwdBarrierVols+ud*deltaVol,t, Barrier*(1.0+ud*eps),1.0,ud) )
					/(eps*Barrier);

		if ( xDigital < 0.0 )
		{
			if( xDigital < -1e-10 ){
				throw("negative value for digital encountered");
			}
			else{
				xDigital = 0.0;
			}
		}

		double windowProb = (Digital-integral)/xDigital;

		if ( windowProb < 0.0 )
		{
			if ( windowProb > -1e-5 ){
				windowProb = 0.0;
			}
			else{
				throw("negative hitting probability encountered");
			}
		}


		stoppingDensity[n] = windowProb;
	}

// calculate call now



	double df, dt;
	price = 0.0;
	for ( n = 1 ; n < npoints; n++ )
	{
		t = mat-0.5*(timeToMat[n]+timeToMat[n-1]);
		dt = (timeToMat[n-1]-timeToMat[n]) / mat ;
		df = exp(-r*t);
		double vprice = ::Bs(Barrier*exp(drift*t), FwdBarrierStrikeVols, t, strike, 1., cp, 0.);

		if ( !payAtTheEnd )
		{
			price += stoppingDensity[n] * vprice + rebate * (dt-stoppingDensity[n]) / df;
		}
		else
		{		
			price += stoppingDensity[n] * (vprice - rebate) + dt * rebate ;
		}

	}

//	if ( payAtTheEnd ){
//		price *= exp(r*mat);
//	}

/*
	if ( HitProb.rows() > 0 )
	{

		double deltaT = mat/(double)HitProb.rows();

		int l = 0,k;
		for ( n = 0 ; n < HitProb.rows(); n++ )
		{
			t = (n+1)*deltaT;
			HitProb[n][0] = t;
			for ( k = l ; k < npoints; k++ )
			{
				if ( double kt = (k+1)*delta_t > t ){
					l = k;
					break;
				}

				HitProb[n][1] += stoppingDensity[k];
			}

		}

		if ( plotCummulativeHitProb )
		{
			for ( n = 0 ; n < HitProb.rows(); n++ )
			{
				if (n) {
					HitProb[n][1] += HitProb[n-1][1];
				}
			}
		}
	}

*/
}



void AMDIP(double& price,double DownBarrier,double spot,double mat,
		   double SpotBarrierVols,double r,double q,
		   double SpotBarrierVolSlope,
		   double FwdBarrierVols,
		   double FwdBarrierVolsSlope,int npoints,int niter)
{


//	set up zeroth order BlackScholes approximation 

	double mu = r-q;
	double deltaT = mat/((double)npoints+1.0);
	double t;
	double volvol = 0.0,deltaVol;
	int nsig=4,ngaussPoints=1;


	CVector hitprob(npoints);
	for (int  n = 0 ; n < npoints; n++ )
	{
		t = (n+1)*deltaT;
		double fwd = spot*exp(mu*t);

		double_knock_out_rebate(hitprob[n],t,spot,fwd,1.0,
							  SpotBarrierVols,0.0,spot*50,DownBarrier,1.0,volvol,
							  nsig,ngaussPoints);
		if ( n ){
			hitprob[n] = hitprob[n]-hitprob[n-1];
		}
	}


	double fwd = spot*exp(mu*mat);

	double eps = 1e-3;
	deltaVol = SpotBarrierVolSlope*eps;
	double sq = 1.0/sqrt(2.0*3.141592654);	


	CVector AmDip(npoints);
	double integral,d1,d2,xmat,y,val;

	for ( int iter = 0 ; iter < niter; iter++ )
	{

		for (int  n = 0 ; n < npoints; n++ )
		{
			t = (n+1)*deltaT;

			double fwd	= spot*exp(mu*t);
			deltaVol	= SpotBarrierVolSlope*eps;

			double
			EurDip = (   Bs(fwd,SpotBarrierVols,t,DownBarrier,1.0,-1) - 
				Bs(fwd,SpotBarrierVols+deltaVol,t,DownBarrier*(1.0-eps),1.0,-1)
			 )/(eps*DownBarrier);

			integral = 0.0;
			for ( int m = 0 ; m < n; m++ )
			{
				xmat = t-(m+1)*deltaT;
				y = FwdBarrierVols*sqrt(xmat);
				d1 = (mu*xmat)/y+0.5*y;
				d2 = d1-y;

				val = normal(-d2)-0.5+
 DownBarrier*exp(-q*xmat)*sqrt(xmat)*sq*exp(-0.5*d1*d1)*FwdBarrierVolsSlope;

				integral += hitprob[m]*val;
			}

			AmDip[n] = 2.0*(EurDip-integral);

			hitprob[n] = AmDip[n];
			if ( n ) {
				hitprob[n] = hitprob[n]-hitprob[n-1];
			}


		}

		price = AmDip[npoints-1];
	}

}



double BlackScholesDIP(double T,double vol,double DownBarrier,double Spot,double drift)
{

	if ( T  < 1e-3 &&  DownBarrier < Spot ){
		return 0.0;
	}

	double lambda	= drift/vol-vol/2.0;
	double l		= 1.0/(vol*sqrt(T))*log(DownBarrier/Spot);
	double alpha    = exp(2.0*lambda*l*sqrt(T));
	double d4		= lambda*sqrt(T)-l;
	double d5       = -lambda*sqrt(T)-l;

//	double price = 1.0-normal(d4)+alpha*(1.0-normal(d5));
	double price = normal(-d4)+alpha*normal(-d5);

	return price;
}


class DataContainer
{
public:

	double stock_forward	;
	double mat				;
	double discount_factor	;
	double cp				;
	double strike			;

public:

	double Tshort;
	double Tlong;
	double DownBarrier;
	double Spot;
	double drift;


};



double BlackScholesCapletPricer(long valuationDate, long RateSettingDate,long RateMaturityDate,double strike,double vol,DayCountConventionEnum & dayCount,MlEqZeroCurve&	discountCurve,int callput=1)
{
	GVector<long> SwapFixingDates(2);

	SwapFixingDates[0] = RateSettingDate;
	SwapFixingDates[1] = RateMaturityDate;

	double res = BlackScholesSwaptionPricer(valuationDate,SwapFixingDates,strike,vol,dayCount,discountCurve,callput);
	return res;
}




/***************************************************************
**	Class  :  
**	Routine: BlackScholesBarrierOption
**	Returns: 
**	Action : 
****************************************************************/




void blackScholesBarrier::initialize(
										double strike,
										double callPut,
										double upperBarrier,
										double lowerBarrier,
										double maturity,
										double vol,
										double forward,
										double spot,
										double rate
						)
{

	m_strike		= strike;
	m_callPut		= callPut;
	m_upperBarrier	= upperBarrier;
	m_lowerBarrier	= lowerBarrier;
	m_maturity		= maturity;
	m_vol			= vol;
	m_forward		= forward;
	m_spot			= spot;
	m_rate			= rate;


	m_lp = 1.0/(m_vol*sqrt(m_maturity))*log(m_upperBarrier/m_spot);
	m_lm = 1.0/(m_vol*sqrt(m_maturity))*log(m_lowerBarrier/m_spot);

	double mu	= log(forward/spot)/m_maturity;
	m_la		= mu/m_vol-m_vol/2.0;
	m_lap		= m_la+m_vol;	
	m_d			= m_lp-m_lm;
	m_k			= log(m_strike/m_spot)/(m_vol*sqrt(m_maturity));


	m_mu		= log(m_forward/m_spot)/m_maturity;

}




double blackScholesBarrier::I(double n,double x )
{

	double res;
	if ( m_callPut == 1 )
	{

		if ( m_strike >= m_lowerBarrier )
		{
			res = normal(m_lp+2.0*n*m_d-x)-normal(m_k+2.0*n*m_d-x);
			res *= exp(-2.0*n*m_d*x);
		}
		else
		{
			res = normal(2.0*m_lp-m_k+2.0*n*m_d+x)-normal(m_lp+2.0*n*m_d+x);
			res *= exp(2.0*x*(n*m_d+m_lp));
		}

		return res;
	}

	if ( m_strike <= m_lowerBarrier )
	{
		res = normal(m_lp+2.0*n*m_d-x)-normal(m_lp+(2.0*n-1.0)*m_d-x);
		res *= exp(-2.0*n*m_d*x);
		return res;
	}
	else 
	{
		res = normal(m_k+2.0*n*m_d-x)-normal(m_lp+(2.0*n-1.0)*m_d-x);
		res *= exp(-2.0*n*m_d*x);
		return res;

	}


}


double blackScholesBarrier::J(double n,double x )
{

	double res;
	if ( m_callPut == 1 )
	{

		if ( m_strike >= m_lowerBarrier )
		{
			res = normal(2.0*m_lp-m_k+2.0*n*m_d+x)-normal(m_lp+2.0*n*m_d+x);
			res *= exp(2.0*x*(n*m_d+m_lp));
		}
		else
		{
			res = normal(m_lp+(2.0*n+1.0)*m_d+x)-normal(m_lp+2.0*n*m_d+x);
			res *= exp(2.0*x*n*m_d+m_lp);
		}

		return res;
	}


	if ( m_strike <= m_lowerBarrier )
	{
		res = normal(m_lp+(2.0*n+1.0)*m_d+x)-normal(m_lp+2.0*n*m_d+x);
		res *= exp(2.0*x*n*m_d+m_lp);
		return res;
	}
	else
	{
		res = normal(m_lp+(2.0*n+1.0)*m_d+x)-normal(2.0*m_lp-m_k+2.0*n*m_d+x);
		res *= exp(2.0*x*(n*m_d+m_lp));
		return res;
	}
}




double blackScholesBarrier::price()
{

	double eps = 1e-10;
	double res;
	double x,xx,nn;
	int nmax = 20;

	int n;
	double sqT = sqrt(m_maturity);


		xx = I(0,m_lap*sqT)-J(0,m_lap*sqT);
		for (  n = 1; n < nmax; n++ )
		{

			nn = (double)n;

			x =		I(nn,m_lap*sqT)-J(nn,m_lap*sqT);
			x +=	I(-nn,m_lap*sqT)-J(-nn,m_lap*sqT);

			xx += x;

			if ( fabs(x)/fabs(xx)< eps ){
				break;
			}

		}
		if ( n == nmax ){
			throw("analytic barrier pricer did not converge");
		}
		
		res = m_callPut*xx*m_spot*exp((m_mu-m_rate)*m_maturity);

		xx = I(0,m_la*sqT)-J(0,m_la*sqT);
		for (  n = 1; n < nmax; n++ )
		{
			nn = (double)n;

			x =		I(nn,m_la*sqT)-J(nn,m_la*sqT);
			x +=	I(-nn,m_la*sqT)-J(-nn,m_la*sqT);

			xx += x;

			if ( fabs(x)/fabs(xx)< eps ){
				break;
			}

		}
		if ( n == nmax ){
			throw("analytic barrier pricer did not converge");
		}
		
		res -= m_callPut*xx*m_strike*exp(-m_rate*m_maturity);
		return res;
}






	
extern void rootfind_underdetermined_solve(const CVector& initial_x, const CVector& bump_x, const CVector& tolerances, 
								   p_dev_underdet_func fn, void* vp, int max_tries, 
								   int max_restarts, const CMatrix& weights, CVector& found_x);


	
	
	





/***************************************************************
**	Class   : calibrateEffectiveLocalVol
**	Function: calibrateEffectiveLocalVol
**	Returns : int
**	Comment : 
****************************************************************/






void calibrateEffectiveLocalVol::calibrateLV(CMatrix& res,CVector& initialGuess)
{

	double FinalTolerance		= 1e-5;
	double StoppingTolerance	= 1e-5;
	double InitialTolerance		= 0.5*FinalTolerance;
	int	   outputFlag			= 0;

    CMatrix NonZeroPartialDerivativesSpecification;
	CMatrix data(m_targetvols.getsize(),2);
	for ( int i = 0 ; i < data.rows(); i++ )
	{
		data[i][0] = m_fixedCalibStrikes[i]/m_shortFwd;
		data[i][1] = m_targetvols[i];
	}


	CMatrix ObjectiveBounds;
	CMatrix xValBounds(initialGuess.getsize(),2);
	
	xValBounds[0][0] = 0.01;
	xValBounds[0][1] = 1.0;
	
	xValBounds[1][0] = -1.5;
	xValBounds[1][1] =  0.0;
	
	xValBounds[2][0] = -0.5;
	xValBounds[2][1] = 1.0;
	
	m_fitLV.initialize(*this,initialGuess,data,
					   xValBounds,ObjectiveBounds,
					   InitialTolerance,FinalTolerance,
					   StoppingTolerance,
					   NonZeroPartialDerivativesSpecification,outputFlag);
	

	CVector out;
	out = initialGuess;

	int maximizeFlag = 0;
	int returnCode;

	m_fitLV.solve(out,maximizeFlag,returnCode);

	CVector prices(m_fixedCalibStrikes.getsize());
	calcOptions(prices,out,m_fixedCalibStrikes);

	m_quadraticLV->reinitialize(out);
	for ( int i = 0 ; i < m_results.rows(); i++ ){
		m_results[i][2] = getLocalVol(m_calibSpots[i]) ;
	}

	for ( int i = 0 ; i < out.getsize(); i++ ){
		m_results[i][3] = out[i] ;
	}


	m_quadraticLV->setCoeff(out);


	res = m_results;


}		
	
	
	
void fitLocalVol::ObjectiveFcn(double* gVals,double* xVals)
{	
	
	for ( int i = 0 ; i < m_currentVols.getsize(); i++ ){
		m_currentVols[i] = xVals[i];
	}

	m_cLV->m_quadraticLV->reinitialize(m_currentVols);
	m_cLV->calcOptions(m_currentVals,m_currentVols,m_cLV->m_fixedCalibStrikes);

	double z;
	double objVal = 0.0;
	for ( int i = 0 ; i < m_currentVals.getsize(); i++ )
	{
		z = m_currentVals[i]-m_cLV->m_targetvols[i];
		objVal += z*z*m_weight[i];
	}

	*gVals = objVal;
}


void fitLocalVol::initialize(calibrateEffectiveLocalVol& calibLV,CVector& initialGuess,CMatrix& data,
					   CMatrix& xValBounds,CMatrix& ObjectiveBounds,
					   double InitialTolerance,double FinalTolerance,
					   double StoppingTolerance,
					   CMatrix& NonZeroPartialDerivativesSpecification,int outputFlag)
{

	m_cLV	= &calibLV;
	m_currentVals.resize(data.rows());
	m_currentVols.resize(initialGuess.getsize());
	m_weight.resize(data.rows());

	for ( int i = 0 ; i < m_weight.getsize(); i++ ){
		m_weight[i] = 1.0;
	}



	if ( ObjectiveBounds.rows() == 0 )
	{
		ObjectiveBounds.resize(1,2);
		ObjectiveBounds[0][0] = -1e99;
		ObjectiveBounds[0][1] = 1e99;
	}


	iVector linearVariableIndex;
	LSGRGSolver::initialize(  initialGuess,linearVariableIndex,
							  xValBounds,ObjectiveBounds,
							  InitialTolerance,FinalTolerance,
							  StoppingTolerance,
							  NonZeroPartialDerivativesSpecification,outputFlag);
}


/***************************************************************
**	Class   : calibrateEffectiveLocalVol
**	Function: calibrate
**	Returns : void
**	Comment : 
****************************************************************/


void calibrateEffectiveLocalVol::init(MlEqAsset& asset,MlEqDateHandle startDate,long endDate,
							   GVector< MlEqStrikeHandle >& calibStrikes,GVector< MlEqStrikeHandle >& calibSpots,CVector& targetvols,
							   CVector& guess,int ngauss,CVector& tanhWingInfo,double finalTol,double stoppingTol,bool initSkewMultsOpt)
{
	m_targetvols	= targetvols;
	m_ngauss		= ngauss;
	m_calibStrikes	= calibStrikes;
	m_fixedCalibStrikes.resize(m_calibStrikes.getsize());
	m_asset			= &asset;

	for ( int i = 0 ; i < m_fixedCalibStrikes.getsize(); i++ )
	{
		MlEqStrike stk;
		MlEqStrike::convertStrikes(stk,*(m_calibStrikes[i]));
		m_fixedCalibStrikes[i] = stk.m_strike;
	}

	m_startDate = startDate;
	m_endDate	= endDate;
	m_shortFwd	= asset.GetForward(asset.GetCurrentDate(),m_startDate->GetDate(),false);
	m_longFwd	= asset.GetForward(asset.GetCurrentDate(),m_endDate,false);
	m_deltaT	= m_startDate->GetYearFraction(m_endDate);
	m_sqdT		= sqrt(m_deltaT);
	m_longMat	= asset.GetDateHandle()->GetYearFraction(m_endDate);

	if ( calibStrikes.getsize() != targetvols.getsize() ){
		throw(" number of target prices supplied does not match number of strikes");
	}

	m_targetPrices.resize(targetvols.getsize());
	for ( int i = 0 ; i < targetvols.getsize(); i++ ){
		m_targetPrices[i] = Bs(m_shortFwd,targetvols[i],m_longMat,m_fixedCalibStrikes[i],1,1);
	}


	double eps = 1e-3;


	double	volatm		= asset.GetVolatility(MlEqStrike(m_shortFwd),m_startDate->GetDate());
	m_shortMat			= asset.GetDateHandle()->GetYearFraction(m_startDate->GetDate());
			
	double  sqT			= sqrt(m_shortMat);
	m_sqShortMat		= sqT;
	m_volatm			= volatm;
	m_fwdvol			= asset.GetVolatility(MlEqStrike(m_longFwd),m_startDate->GetDate(),m_endDate);

	
//  start setting up calibration points in local vol
	
	int ndim = calibSpots.getsize();
	m_calibSpots.resize(ndim);

	m_results.resize(ndim,4);
	for ( int i = 0 ; i < ndim ; i++ ){
		m_results[i][0]		= calibSpots[i]->m_strike;
	}

	for ( int i = 0 ; i < ndim; i++ )
	{
		MlEqStrike stk;
		MlEqStrike::convertStrikes(stk,*calibSpots[i]);
		m_calibSpots[i] = stk.m_strike;
		m_results[i][1]  = m_calibSpots[i];
	}
	

//  start initiializing local vol surface

	if ( guess.getsize() == 0 )
	{
		m_guess.resize(ndim);
		for ( int i = 0 ; i < ndim ; i++ ){
			m_guess[i] = asset.GetVolatility(MlEqStrike(m_calibSpots[i]),m_startDate->GetDate(),m_endDate);
		}
	}
	else{
		m_guess	= guess;
	}

	
	int n	= m_guess.getsize();
	CVector xData(n);
	CVector yData(n);
		

	CMatrix xValBounds(3,2);

	xValBounds[0][0] = 0.0;
	xValBounds[0][1] = 1.0;
	xValBounds[1][0] = -1.50;
	xValBounds[1][1] = +0.0;
	xValBounds[2][0] = 0.0;
	xValBounds[2][1] = 1.0;

	double cL ;
	double cR ;
	double yPower ;
	int	   addTanhWings;

	if ( tanhWingInfo.getsize() )
	{
		addTanhWings	= tanhWingInfo[0];
		cL				= tanhWingInfo[1];
		cR				= tanhWingInfo[2];
		yPower			= tanhWingInfo[3];
	}
	else
	{
		addTanhWings	= 1.0;
		cL				= 1.5;
		cR				= 5.0;
		yPower			= 2.0;
	}

// second argument 

	for ( int i = 0 ; i < n; i++ ){
		xData[i] = log(m_calibSpots[i]/m_shortFwd)/(m_volatm*m_sqdT);
	}


	m_quadraticLV = new MlEqAnalyticCurveWithTanhWingEdge();
	m_quadraticLV->initialize(xData,m_guess,xValBounds,addTanhWings,cL,cR,yPower,finalTol,stoppingTol);


	double lowerEdge = log(m_calibSpots[0]/m_shortFwd)/(m_volatm*m_sqShortMat);
	double upperEdge = log(m_calibSpots[ndim-1]/m_shortFwd)/(m_volatm*m_sqShortMat);;

	CVector coeff;

	int npoints;
	if ( m_calibSpots.getsize() < 5 ){
		npoints = 5;
	}
	else{
		npoints = m_calibSpots.getsize();
	}

	coeff = m_guess; 
	m_quadraticLV->initialize(lowerEdge,upperEdge,addTanhWings, cL, cR, yPower,npoints,coeff);

	CVector test(m_calibSpots.getsize());

//  set up future skews

	int xdim = 7;

	GVector<CVector> xDat(1);
	xDat[0].resize(xdim);

	xDat[0][0] =  0.6;
	xDat[0][1] =  0.7;
	xDat[0][2] =  0.8;
	xDat[0][3] =  0.9;
	xDat[0][4] =  1.0;
	xDat[0][5] =  1.2;
	xDat[0][6] =  1.4;

	GVector< CVector > yDat(1);
	yDat[0].resize(xdim);
	for ( int i = 0 ; i < xdim ; i++ )
	{
		MlEqStrikeHandle stk = new 	MlEqForwardBasedStrike(*m_asset,xDat[0][i],m_endDate,m_startDate,false);
		yDat[0][i] = asset.GetVolatility(*stk,m_startDate->GetDate(),m_endDate)/m_fwdvol;
	}


	if ( initSkewMultsOpt )
	{
		cL = 5.0;
		cR = 5.0;
		yPower = 1.0;
		addTanhWings = 1;

		m_skewMultipliers = new MlEqMonotonicSplineInterpolator(xDat,yDat,addTanhWings,cL,cR,yPower);

	}

	
//  start to set up numerical integrators
	

	m_ngauss = ngauss;
	m_gaussPoints.resize(ngauss);
	m_gaussWeights.resize(ngauss);

	sqT			= sqrt(m_shortMat);
	double nstdev = 3.0;
	double lower = exp(volatm*nstdev*sqT*(-nstdev));
	double upper = exp(volatm*nstdev*sqT*(nstdev));
	
	if ( lower > m_calibSpots[0] ){
		throw("incorrect dimensioning of effctive local vol calibration");
	}

	if ( upper > m_calibSpots[m_calibSpots.getsize()-1] ){
		throw("incorrect dimensioning of effctive local vol calibration");
	}
	
	MlEqMaths::dGauleg(lower,upper,m_gaussPoints,m_gaussWeights,ngauss,true);				

	m_vols.resize(m_ngauss);
	for ( int i = 0 ; i < m_ngauss; i++ )
	{
		m_vols[i] = m_asset->GetVolatility(MlEqStrike(m_gaussPoints[i]*m_shortFwd),m_asset->GetCurrentDate(),m_startDate->GetDate());
	}

}	
/***************************************************************
**	Class   : calibrateEffectiveLocalVol
**	Function: calibrate
**	Returns : void
**	Comment : 
****************************************************************/



void calibrateEffectiveLocalVol::init(MlEqAsset& asset,MlEqAnalyticCurveWithTanhWingEdgeHandle& quadraticInterp,MlEqDateHandle startDate,long endDate,
							   GVector< MlEqStrikeHandle >& calibStrikes,GVector< MlEqStrikeHandle >& calibSpots,CVector& targetvols,
							   CVector& guess,int ngauss,bool initSkewMultsOpt)
{
	m_targetvols	= targetvols;

	m_calibStrikes	= calibStrikes;
	m_fixedCalibStrikes.resize(m_calibStrikes.getsize());
	m_asset			= &asset;

	for ( int i = 0 ; i < m_fixedCalibStrikes.getsize(); i++ )
	{
		MlEqStrike stk;
		MlEqStrike::convertStrikes(stk,*(m_calibStrikes[i]));
		m_fixedCalibStrikes[i] = stk.m_strike;
	}

	m_startDate = startDate;
	m_endDate	= endDate;
	m_shortFwd	= asset.GetForward(asset.GetCurrentDate(),m_startDate->GetDate(),false);
	m_longFwd	= asset.GetForward(asset.GetCurrentDate(),m_endDate,false);
	m_deltaT	= m_startDate->GetYearFraction(m_endDate);
	m_sqdT		= sqrt(m_deltaT);
	m_longMat	= asset.GetDateHandle()->GetYearFraction(m_endDate);

	if ( calibStrikes.getsize() != targetvols.getsize() ){
		throw(" number of target prices supplied does not match number of strikes");
	}

	m_targetPrices.resize(targetvols.getsize());
	for ( int i = 0 ; i < targetvols.getsize(); i++ ){
		m_targetPrices[i] = Bs(m_shortFwd,targetvols[i],m_longMat,m_fixedCalibStrikes[i],1,1);
	}


	double eps = 1e-3;

	double	volatm		= asset.GetVolatility(MlEqStrike(m_shortFwd),m_startDate->GetDate());
	m_shortMat			= asset.GetDateHandle()->GetYearFraction(m_startDate->GetDate());
			
	int xtest = asset.GetCurrentDate();
	xtest = m_startDate->GetDate();

	double  sqT			= sqrt(m_shortMat);
	m_sqShortMat		= sqT;
	m_volatm			= volatm;
	m_fwdvol			= asset.GetVolatility(MlEqStrike(m_longFwd),m_startDate->GetDate(),m_endDate);

	
//  start setting up calibration points in local vol
	
	int ndim = calibSpots.getsize();
	m_calibSpots.resize(ndim);

	m_results.resize(ndim,4);
	for ( int i = 0 ; i < ndim ; i++ ){
		m_results[i][0]		= calibSpots[i]->m_strike;
	}

	for ( int i = 0 ; i < ndim; i++ )
	{
		MlEqStrike stk;
		MlEqStrike::convertStrikes(stk,*calibSpots[i]);
		m_calibSpots[i] = stk.m_strike;
		m_results[i][1]  = m_calibSpots[i];
	}
	

//  start initiializing local vol surface

	m_quadraticLV = quadraticInterp;

	if ( guess.getsize() == 0 )
	{
		const CVector& guess = m_quadraticLV->getInitialGuess();
		m_guess = guess;
	}
	else{
		m_guess	= guess;
	}

	CVector test(m_calibSpots.getsize());

//  set up future skews

	int xdim = 7;

	GVector<CVector> xDat(1);
	xDat[0].resize(xdim);

	xDat[0][0] =  0.6;
	xDat[0][1] =  0.7;
	xDat[0][2] =  0.8;
	xDat[0][3] =  0.9;
	xDat[0][4] =  1.0;
	xDat[0][5] =  1.2;
	xDat[0][6] =  1.4;

	GVector< CVector > yDat(1);
	yDat[0].resize(xdim);
	for ( int i = 0 ; i < xdim ; i++ )
	{
		MlEqStrikeHandle stk = new 	MlEqForwardBasedStrike(*m_asset,xDat[0][i],m_endDate,m_startDate,false);
		yDat[0][i] = asset.GetVolatility(*stk,m_startDate->GetDate(),m_endDate)/m_fwdvol;
	}


	if ( initSkewMultsOpt )
	{
		double cL = 5.0;
		double cR = 5.0;
		double yPower = 1.0;
		double addTanhWings = 1;

		m_skewMultipliers = new MlEqMonotonicSplineInterpolator(xDat,yDat,addTanhWings,cL,cR,yPower);

	}
	
//  start to set up numerical integrators
	

	m_ngauss = ngauss;
	m_gaussPoints.resize(ngauss);
	m_gaussWeights.resize(ngauss);



	double nstdev = 3.0;
	double lower = exp(volatm*sqT*(-nstdev));
	double upper = exp(volatm*sqT*(nstdev));
	
	if ( lower > m_calibSpots[0] ){
		throw("incorrect dimensioning of effctive local vol calibration");
	}

	if ( upper > m_calibSpots[m_calibSpots.getsize()-1] ){
		throw("incorrect dimensioning of effctive local vol calibration");
	}
	
	MlEqMaths::dGauleg(lower,upper,m_gaussPoints,m_gaussWeights,ngauss,true);				

	m_vols.resize(m_ngauss);
	for ( int i = 0 ; i < m_ngauss; i++ )
	{
		m_vols[i] = m_asset->GetVolatility(MlEqStrike(m_gaussPoints[i]*m_shortFwd),m_asset->GetCurrentDate(),m_startDate->GetDate());
	}

}	

void calibrateEffectiveLocalVol::calcOptions(CVector& prices, const CVector& vols,GVector< MlEqStrikeHandle >& strikes)
{
	int n = strikes.getsize();
	CVector fixedStrikes(n);
	MlEqStrike stk;

	for ( int i = 0 ; i < n; i++ )
	{	
		MlEqStrike::convertStrikes(stk,*strikes[i]);
		fixedStrikes[i] = stk.m_strike;
	}

	calcOptions(prices,vols,fixedStrikes);

}


double calibrateEffectiveLocalVol::getLocalVol(double xspot)
{
	double vol;

	double x = log(xspot/m_shortFwd)/(m_volatm*m_sqdT);
	vol = m_quadraticLV->getValue(x);
	return vol;
}


void calibrateEffectiveLocalVol::calcOptions(CVector& prices, const CVector& vols,CVector& fixedStrikes)
{	
	
	CVector xVol;
	xVol = vols;
	m_quadraticLV->reinitialize(xVol);
	
	int n = prices.getsize();
	if ( n != fixedStrikes.getsize() )
	{
		n = fixedStrikes.getsize();
		prices.resize(n);
	}
	
	for ( int i = 0 ; i < n ; i++ ){
		prices[i] = 0.0;

	}
		
	double vol,z,spot,fwd,skewMult,x,density,eps,xspot;
	double fwdCall;

	CVector xcall(m_ngauss);
	CVector xput(m_ngauss);

	eps = 1e-3;
	for ( int i = 0 ; i < n ; i++ )
	{

		spot = m_shortFwd;
		fwd  = spot*m_longFwd/m_shortFwd;

		vol = getLocalVol(spot);
//		x = log(fixedStrikes[i]/fwd)/(m_fwdvol*m_sqdT);
		x = fixedStrikes[i]/fwd;

		skewMult = m_skewMultipliers->getValue(x);
		vol *= skewMult;

		fwdCall = Bs(fwd,vol,m_longMat-m_shortMat,fixedStrikes[i],1,1);
		prices[i] = fwdCall;

		double test = 0.0;
		for ( int k = 0 ; k < m_ngauss; k++ )
		{

			spot = m_gaussPoints[k]*m_shortFwd;
			fwd  = spot*m_longFwd/m_shortFwd;


			if ( i == 0 )
			{
				xcall[k] = Bs(m_shortFwd,m_vols[k],m_shortMat,spot,1,1);
				xput[k]  = Bs(m_shortFwd,m_vols[k],m_shortMat,spot,1,-1);
			}

//			x = log(fixedStrikes[i]/fwd)/(m_fwdvol*m_sqdT);
			x = fixedStrikes[i]/fwd;

			vol = getLocalVol(spot);
			skewMult = m_skewMultipliers->getValue(x);
			vol *= skewMult;

			density =	-2.0*Bs(fwd,vol,m_longMat-m_shortMat,fixedStrikes[i],1,1);
			xspot = spot*(1.0+eps);
			fwd   = xspot*m_longFwd/m_shortFwd;
			vol = getLocalVol(xspot);
//			x = log(fixedStrikes[i]/fwd)/(m_volatm*m_sqdT);
			x = fixedStrikes[i]/fwd;

			skewMult = m_skewMultipliers->getValue(x);
			vol *= skewMult;

			density +=	Bs(fwd,vol,m_longMat-m_shortMat,fixedStrikes[i],1,1);
				

			xspot = spot*(1.0-eps);
			fwd   = xspot*m_longFwd/m_shortFwd;
			vol = getLocalVol(xspot);
//			x = log(fixedStrikes[i]/fwd)/(m_volatm*m_sqdT);
			x = fixedStrikes[i]/fwd;

			skewMult = m_skewMultipliers->getValue(x);

			vol *= skewMult;
			density +=	Bs(fwd,vol,m_longMat-m_shortMat,fixedStrikes[i],1,1);
			
			density /= pow(spot*eps,2.0);

			if ( spot >= m_shortFwd ){
				z =	 density*xcall[k];
			}
			else{
				z =	 density*xput[k];
			}

			z *= m_gaussWeights[k]*m_shortFwd;

			test += density*m_gaussWeights[k]*m_shortFwd;
			prices[i] += z;

		}
	}

	double accuracy		= 1e-5;
	double lower_x		= 0.0;
	double upper_x		= 1.5;


	for ( int istrike = 0.0 ; istrike < n; istrike++ )
	{


		prices[istrike] = MlEqBSImpliedVol(
								prices[istrike],		
								m_longFwd,			
								m_longMat,	
								fixedStrikes[istrike],
								1.0,
								1,		
								eROOTFIND_BRENT_GROWRANGE,	
								accuracy,
								lower_x,
								upper_x	
								);
	}

}	


void calibrateEffectiveLocalVolFull(CMatrix& result,MlEqAsset& asset, GVector <MlEqAnalyticCurveWithTanhWingEdgeHandle >& quadraticInterpolator,
									MlEqDateHandle startDate,GVector< long > endDates,CMatrix& guess,int ngauss,bool usePrevRes,
									GVector< GVector< MlEqStrikeHandle > >& calibStrikes,GVector< GVector< MlEqStrikeHandle > >& calibSpots)
{

	int ndates = quadraticInterpolator.getsize();

	if ( endDates.getsize()-1 != ndates ){
		throw("inconsistent input encountered: number of calibration dates must coincide with number of interpolators provided");
	}

	if ( guess.rows() ){
		if ( guess.rows() != ndates ){
			throw("incorrect inputs encountered");
		}
	}

	if ( calibStrikes.getsize() != 1 ){
		if ( calibStrikes.getsize() != ndates ){
			throw("inconsistent set up in calibration strikes");
		}
	}

	if ( calibSpots.getsize() != 1 ){
		if ( calibSpots.getsize() != ndates ){
			throw("inconsistent set up in calibration strikes");
		}
	}

	result.resize(ndates,3);

	CMatrix res;

	for ( int idate = 0 ; idate < ndates; idate++ )
	{
			calibrateEffectiveLocalVol cEffLv;
			startDate->PutDate(endDates[idate]);

			CVector currentGuess;

			if ( guess.rows() )
			{
				if ( guess.rows() == ndates ){
					currentGuess = guess[idate];
				}
				else{
					currentGuess = guess[0];
				}
			}
			else
			{
				int jdate = idate;
				if ( usePrevRes && idate ){
					jdate = idate-1;
				}

				const CVector& xguess =  quadraticInterpolator[idate]->getCoeff();
				currentGuess = xguess;
			}

	
			GVector< MlEqStrikeHandle > currentCalibStrikes;
			if ( calibStrikes.getsize() == ndates ){
				currentCalibStrikes = calibStrikes[idate];
			}
			else{
				currentCalibStrikes = calibStrikes[0];
			}

			GVector< MlEqStrikeHandle > currentCalibSpots;
			if ( calibSpots.getsize() == ndates ){
				currentCalibSpots = calibSpots[idate];
			}
			else{
				currentCalibSpots = calibSpots[0];
			}

			CVector targetvols(currentCalibStrikes.getsize());
			for ( int istrike = 0 ; istrike < targetvols.getsize(); istrike++ ){
				targetvols[istrike] = asset.GetVolatility(*currentCalibStrikes[istrike], endDates[idate+1]);
			}

int test = startDate->GetDate();

			cEffLv.init(asset,quadraticInterpolator[idate],startDate,endDates[idate+1],
						currentCalibStrikes,currentCalibSpots,targetvols,currentGuess,ngauss);

			cEffLv.calibrateLV(res,currentGuess);
			CVector coeff(3);

			for ( int i = 0 ; i < coeff.getsize(); i++ ){
				coeff[i] = res[i][3];
			}

			quadraticInterpolator[idate]->reinitialize(coeff);

			result[idate] = coeff;
	}
}




void calibrateEffectiveLocalVolGridFull(CMatrix& result,MlEqAsset& asset, GVector <MlEqAnalyticCurveWithTanhWingEdgeHandle >& quadraticInterpolator,
									MlEqDateHandle startDate,GVector< long > endDates,CMatrix& guess,bool usePrevRes,
									GVector< GVector< MlEqStrikeHandle > >& calibStrikes,GVector< GVector< MlEqStrikeHandle > >& calibSpots,
									double lowSpot,double highSpot,int nx,int nt)
{


//	create pde
	pdeLocalVolEffectiveHandle pdeDriver = new pdeLocalVolEffective();

	double spot = asset.GetSpot(asset.GetCurrentDate());
	double lx = log(lowSpot/spot);
	double ux = log(highSpot/spot);

	const int keep = PdeEnums::PDE_KEEP_LAST;
	pdeDriver->init(quadraticInterpolator,asset,startDate,
					endDates[endDates.getsize()-1],endDates,
					NULL,1, nx,&nx,lx,nt,ux,keep);

	int ndates = quadraticInterpolator.getsize();

	if ( endDates.getsize() != ndates ){
		throw("inconsistent input encountered: number of calibration dates must coincide with number of interpolators provided");
	}

	if ( guess.rows() ){
		if ( guess.rows() != ndates ){
			throw("incorrect inputs encountered");
		}
	}

	if ( calibStrikes.getsize() != 1 ){
		if ( calibStrikes.getsize() != ndates ){
			throw("inconsistent set up in calibration strikes");
		}
	}

	if ( calibSpots.getsize() != 1 ){
		if ( calibSpots.getsize() != ndates ){
			throw("inconsistent set up in calibration strikes");
		}
	}


// start setting up pde


	CVector CallPut(calibStrikes.getsize());
	for ( int i = 0 ; i < CallPut.getsize(); i++ ){
		CallPut[i] = 1.0;
	}



	result.resize(ndates,3);
	CMatrix res;	

	//77777777777777777777
// test
			plainVanillaOptionHelper plainVanilla;
			plainVanilla.initialize(CallPut,calibStrikes[ndates-1],spot);
//			pdeDriver->reinitialize_nd(calibStrikes[ndates-1].getsize(),&plainVanilla);

			pdeDriver->init(quadraticInterpolator,asset,startDate,
							endDates[endDates.getsize()-1],endDates,
							&plainVanilla,calibStrikes[ndates-1].getsize(), nx,&nx,lx,nt,ux,keep);



			CVector xresult(calibStrikes[ndates-1].getsize());
			pdeDriver->pde_integrate(xresult);


	//777777777777777

	
	for ( int idate = 0 ; idate < ndates; idate++ )
	{				
				
			
			CVector currentGuess;
			if ( guess.rows() )
			{
				if ( guess.rows() == ndates ){
					currentGuess = guess[idate];
				}
				else{
					currentGuess = guess[0];
				}
			}
			else
			{
				int jdate = idate;
				if ( usePrevRes && idate ){
					jdate = idate-1;
				}

				const CVector& xguess =  quadraticInterpolator[idate]->getCoeff();
				currentGuess = xguess;
			}


			GVector< MlEqStrikeHandle > currentCalibStrikes;
			if ( calibStrikes.getsize() == ndates ){
				currentCalibStrikes = calibStrikes[idate];
			}
			else{
				currentCalibStrikes = calibStrikes[0];
			}
			
			GVector< MlEqStrikeHandle > currentCalibSpots;
			if ( calibSpots.getsize() == ndates ){
				currentCalibSpots = calibSpots[idate];
			}
			else{
				currentCalibSpots = calibSpots[0];
			}
			
			CVector targetvols(currentCalibStrikes.getsize());
			for ( int istrike = 0 ; istrike < targetvols.getsize(); istrike++ ){
				targetvols[istrike] = asset.GetVolatility(*currentCalibStrikes[istrike], endDates[idate+1]);
			}

			
			plainVanillaOptionHelper plainVanilla;
			plainVanilla.initialize(CallPut,calibStrikes[idate],spot);
			pdeDriver->reinitialize_nd(calibStrikes[idate].getsize(),&plainVanilla);
			
			calibrateEffectiveLocalVolGrid cEffLv;
			
			cEffLv.init(idate,pdeDriver,asset,quadraticInterpolator[idate],startDate,endDates[idate],
						currentCalibStrikes,currentCalibSpots,targetvols);

			cEffLv.calibrateLV(res,currentGuess);
			
			CVector coeff(3);
			for ( int i = 0 ; i < coeff.getsize(); i++ ){
				coeff[i] = res[i][3];
			}
			
			quadraticInterpolator[idate]->reinitialize(coeff);			
			result[idate] = coeff;
			
	}		
		
}		
		

/***************************************************************
**	Class   : calibrateEffectiveLocalVolGrid
**	Function: calcOptions
**	Returns : void
**	Comment : 
****************************************************************/



void calibrateEffectiveLocalVolGrid::calcOptions(CVector& prices, const CVector& vols,CVector& fixedStrikes)
{

	m_pdeDriver->reinitialize(m_idate,vols);
	const GVector<long>& mapinv = m_pdeDriver->getMapInv();

	m_pdeDriver->pde_integrate(prices,mapinv[m_idate]);

	double accuracy		= 1e-5;
	double lower_x		= 0.0;
	double upper_x		= 1.5;

	for ( int istrike = 0.0 ; istrike < fixedStrikes.getsize(); istrike++ )
	{

		prices[istrike] = MlEqBSImpliedVol(
								prices[istrike],		
								m_longFwd,			
								m_longMat,	
								fixedStrikes[istrike],
								1.0,
								1,		
								eROOTFIND_BRENT_GROWRANGE,	
								accuracy,
								lower_x,
								upper_x	
								);
	}

}



/***************************************************************
**	Class   : calibrateEffectiveLocalVolGrid
**	Function: calcOptions
**	Returns : void
**	Comment : 
****************************************************************/


void calibrateEffectiveLocalVolGrid::init(int idate,pdeLocalVolEffectiveHandle pdeDriver,MlEqAsset& asset,MlEqAnalyticCurveWithTanhWingEdgeHandle& quadraticInterp,MlEqDateHandle startDate,long endDate,
										  GVector< MlEqStrikeHandle >& calibStrikes,GVector< MlEqStrikeHandle >& calibSpots,CVector& targetvols)
{


	CVector guess;

	calibrateEffectiveLocalVol::init(asset,quadraticInterp,startDate,endDate,
							   calibStrikes,calibSpots,targetvols,
							   guess,0,false);

	m_pdeDriver = pdeDriver;
	m_idate = idate;
  	m_pdeDriver->reinitialize(m_idate,m_quadraticLV);

}

/***************************************************************
**	Class   : calibrateEffectiveLocalVolGrid
**	Function: calcOptions
**	Returns : void
**	Comment : 
****************************************************************/

void calibrateEffectiveLocalVolGrid::calcOptions(CVector& prices, const CVector& vols,GVector< MlEqStrikeHandle >& strikes)
{

	int n = strikes.getsize();
	CVector fixedStrikes(n);
	MlEqStrike stk;

	for ( int i = 0 ; i < n; i++ )
	{	
		MlEqStrike::convertStrikes(stk,*strikes[i]);
		fixedStrikes[i] = stk.m_strike;
	}

	calcOptions(prices,vols,fixedStrikes);





}
