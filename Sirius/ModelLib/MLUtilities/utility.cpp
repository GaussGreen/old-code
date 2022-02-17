
#include "stdafx.h"
#include "mleqobjects.h"

#include "utility.h"
#include "cMatrix.h"
#include <assert.h>
#include "montecarlo.h"
#include <math.h>
#include "SciSobol.h"
#include "MlEqPde.h"

const static double EPS = 1e-10;

void ConvertConstantNotionalRamStrike(NotionalTypeEnum* pType, long nStartDate, long nValuationDate, double* pfNotional, double* pfStrike)
{
	if (*pType != RAMConstantNotional) return;	

	if (nStartDate > nValuationDate){
		*pfStrike /= 100.0;
	} else {
		*pType = ConstantShare;
		*pfNotional /= 100.0;		

	}

}


void trim(std::string* psz)
{
	size_t nFirst = psz->find_first_not_of(" \n\r\t\f\v");
	size_t nLast = psz->find_last_not_of(" \n\r\t\f\v");	
	if (nFirst == psz->npos){
		psz->erase();
	} else {
		psz->assign(psz->substr(nFirst, nLast + 1 - nFirst));
	}	
}






double bridge(double delta_t,double U,double L,double X_i,double X_f,double vol)
{	

//	this is the probability of not hitting the barrier
	
	if ( MlEqMaths::Min(X_i,X_f) <= L+EPS )
		return 0.0;
	
	if ( MlEqMaths::Max(X_i,X_f) >= U-EPS )
		return 0.0;
	
	if ( L < EPS )
	{
		L = 1e-10;
	}

	int n;
	int Nmax = 40;
	double eps = 1e-20;
	double dn,zz,y;
	double z = 0.0;
	
	double a = 1.0/vol*log(U/X_i);
	double c = 1.0/vol*log(U/L);
	double x = 1.0/vol*log(X_f/X_i);
	
	for ( n = 0 ; n < Nmax; n++ )
	{
		dn = (double) n;
		zz = exp(2.0*dn*c*(x-dn*c)/delta_t)-exp(2.0*(x+dn*c-a)*(a-dn*c)/delta_t);
		y = zz;
	
		if ( n != 0 )
		{
			dn = -dn;
			zz += exp(2.0*dn*c*(x-dn*c)/delta_t)-exp(2.0*(x+dn*c-a)*(a-dn*c)/delta_t);
			y += zz;
		}
		z += zz;
		if ( fabs(y)/z < eps || fabs(y) < 1e-8 )
		{
			if ( z < 0.0 )
				z = 0.0;
			else if ( z > 1.0 )
			{
				z = 1.0;
			}
			break;
		}
	}
	
	
	if ( n == Nmax )
		throw("bridge function did not converge");
	
	return z;
}	



static double BS(double forward,double vol,double maturity,double strike,double discount_factor,int cp)
{	
	double res,d1,d2;
	
	if ( vol < EPS  || maturity < 1e-6)
	{
		res = discount_factor*MlEqMaths::Max((double)cp*(forward-strike),0.0);
		return res;
	}

	if ( strike < 0.0 && forward > 0.0)
	{
		res = discount_factor*MlEqMaths::Max((double)cp*(forward-strike),0.0);
		return res;
	}

	double volsqrt= vol*sqrt(maturity);
	
	if ( forward/strike < 0.0 )
	{
		res = 0.0;
		return res;
	}
	
	d1 = (log(forward/strike)+0.5*volsqrt*volsqrt)/volsqrt;
	d2 = d1 - volsqrt;
	
	if ( cp == 1)
		res = discount_factor*( forward*normal(d1)-strike*normal(d2) );
	else if ( cp == -1 )
		res = discount_factor*( -forward*normal(-d1)+strike*normal(-d2) );
	else
		assert(0);
		
	return res;
}	

	
double Bs(double forward,double vol,double maturity,double strike,double discount_factor,int cp,double blend)
{

	if ( ( 1.0-blend ) < 0.001 )
	{
		blend = 0.999;
	}

	if ( blend < 0.001 )
	{
		blend = 0.001;
	}

	double std	 = vol*sqrt(maturity);
	double cb	 = normal(std/2.0)-normal((1.0-blend)*std/2.0);
		   cb	/= normal((1.0-blend)*std/2.0)-0.5;
		   cb	*= forward;

	if ( maturity < 1e-6 ){
		cb = 0.0;
	}

	return BS(forward+cb,vol*(1.0-blend),maturity,strike+cb,discount_factor,cp);
}


double Bs(double forward,double vol,double maturity,double strike,double discount_factor,int cp,double blend,double volvol,int nsig /*=4*/,int npoints /*=12*/)
{

	int N = npoints;
	
	if ( N <= 1 ){
		return Bs(forward,vol,maturity,strike,discount_factor,cp,blend );
	}
	if ( nsig < 1 ){
		throw(" number of standart deviations entered cannot be smaller than one");
	}


	CVector gaussWeight(N);
	CVector gaussPoint(N);

	double mat = maturity;
	double sqrt_mat = sqrt(maturity);
	double stochFac,prob,vols;
	double sq = 1.0/sqrt(2.0*3.141592654);	
			
	MlEqMaths::dGauleg(-nsig, nsig,gaussPoint,gaussWeight,N,true);

	double res = 0.0;
	for (int i = 0 ; i < N; i++ )
	{
		stochFac	=	exp(-0.5*volvol*volvol*mat+volvol*sqrt_mat*gaussPoint[i]);
		prob		=	gaussWeight[i]*sq*exp(-0.5*gaussPoint[i]*gaussPoint[i]);	
		vols		=	vol*stochFac;
		res		    +=	prob*Bs(forward,vols,maturity,strike,discount_factor,cp,blend);
	}

	return res;
}

void double_knock_out_rebate(double & price,double mat,double spot,double forward,
						  double discount_factor,double vol,double blend,double barrier_up,
						  double barrier_do,double rebate,double volvol,int nsig /*=4*/,int npoints /*=12*/)
{

	int N = npoints;
	
	if ( N < 1 || fabs(volvol) < 1e-3 ){
		double_knock_out_rebate(price,mat,spot,forward,
						  discount_factor,vol,blend,barrier_up,
						  barrier_do,rebate);

		return;
	}
	if ( nsig < 1 ){
		throw(" number of standart deviations entered cannot be smaller than one");
	}


	CVector gaussWeight(N);
	CVector gaussPoint(N);


	double sqrt_mat = sqrt(mat);
	double stochFac,prob,vols;
	double sq = 1.0/sqrt(2.0*3.141592654);	
			
	MlEqMaths::dGauleg(-nsig, nsig,gaussPoint,gaussWeight,N,true);

	price = 0.0;
	double res ;
	for (int i = 0 ; i < N; i++ )
	{
		stochFac	=	exp(-0.5*volvol*volvol*mat+volvol*sqrt_mat*gaussPoint[i]);
		prob		=	gaussWeight[i]*sq*exp(-0.5*gaussPoint[i]*gaussPoint[i]);	
		vols		=	vol*stochFac;
		double_knock_out_rebate(res,mat,spot,forward,
						  discount_factor,vols,blend,barrier_up,
						  barrier_do,rebate);

		price += prob*res;
	}


}

void double_knock_out_rebate(double & price,double mat,double spot,double forward,
						  double discount_factor,double vol,double blend,double barrier_up,
						  double barrier_do,double rebate)
{


	if ( ( 1.0-blend ) < 0.001 ){
		blend = 0.999;
	}

	if ( blend < 0.001 ){
		blend = 0.001;
	}

	double std	 = vol*sqrt(mat);
	double cb	 = normal(std/2.0)-normal((1.0-blend)*std/2.0);
		   cb	/= normal((1.0-blend)*std/2.0)-0.5;
		   cb	*= forward;



	double_knock_out_rebate(price,mat,spot+cb,forward+cb,
						  discount_factor,vol*(1.0-blend),barrier_up+cb,
						  barrier_do+cb,rebate);

	
}

void double_knock_out_rebate(double & price,double mat,double spot,double forward,
						  double discount_factor,double vol,double barrier_up,
						  double barrier_do,double rebate)


{

	int n;
	double dn,la,xi,xip,lp,lm,alphap,d4p,d4m,mu,x,xx,d5p,d5m;

	if ( ( spot > barrier_up )|| ( spot < barrier_do ) )
	{
		price = rebate; //AB: before it was {price = 0;} , as we are already knock in we pay rebate 
		return;
	}


	// sometimes there is overflow with too big blend (near to 1) cuz 
	// vol is about zero, and forward is near to spot 

	double eps = 1e-10;

	mu = log(forward/spot)/mat;

	if(barrier_do<0.0000001) barrier_do=0.00001;//AB: this line added
	if(barrier_up/spot>100) 
		barrier_up=100*spot;//AB: this line added cuz overflow

	lp = 1.0/(vol*sqrt(mat))*log(barrier_up/spot);



	lm = 1.0/(vol*sqrt(mat))*log(barrier_do/spot);
	la = mu/vol-vol/2.0;
	d4p = la*sqrt(mat)-lp;
	d4m = la*sqrt(mat)-lm;
	alphap = exp(2.0*la*lp*sqrt(mat));
	d5p = -la*sqrt(mat)-lp;
	d5m = -la*sqrt(mat)-lm;
	double y;

	int Nmax = 40;
	
	x = 0.0;
	for ( n = 0 ; n < Nmax; n++ )
	{

		dn = (double)n;
		xi = 2.0*dn*(lp-lm);
		xip = xi+2.0*(lp-lm);

		xx = exp(-la*sqrt(mat)*xi)*(normal(xi-d4p)-normal(xi-d4m))
			- alphap*exp(la*sqrt(mat)*xi)*(normal(xip-d5m)-normal(xi-d5p));

		y = xx;
		x += xx;

		if ( n != 0)
		{
			dn = -dn;
			xi = 2.0*dn*(lp-lm);
			xip = xi+2.0*(lp-lm);

			xx = exp(-la*sqrt(mat)*xi)*(normal(xi-d4p)-normal(xi-d4m))
				- alphap*exp(la*sqrt(mat)*xi)*(normal(xip-d5m)-normal(xi-d5p));

			y += xx;
			x += xx;
		}

		if ( fabs(y/x) < eps )
			break;
	}

	if ( n > Nmax ){
		throw("maximum number of iterations exceeded in barrier calculation");
	}

	price = (1.0-x)*rebate*discount_factor;	
}	




curve::curve(CVector& xArray,CVector& yArray)
{
	initialize(xArray,yArray);
}

int curve::initialize(CVector& xArray,CVector& yArray)
{
	m_xVals = xArray;
	m_yVals = yArray;

	m_method = 0;

	return 0;
}

		
double curve::getYValue(double x)
{
	double z;
	z = MlEqMaths::linearInterp(m_xVals,m_yVals,x,m_method);
	return z;
}

double curve::getXValue(int i)
{
	return m_xVals[i];
}
	
int curve::getsize()
{
	int n = m_xVals.getsize();
	return n;
}


curve& curve::operator=(const curve& rhs)
{
	
	this->m_method = rhs.m_method;
	this->m_xVals = rhs.m_xVals;
	this->m_yVals = rhs.m_yVals;

	return *this;
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


bool implied_vol_fcn(double x,void* vp,double* f)
{			
	DataContainer* p = (DataContainer*)vp;

    double res = Bs(p->stock_forward,x,p->mat,p->strike,p->discount_factor,p->cp);
	*f =  res;
	return true;
}	


/****************************************************************
**	class:		EqSpBSImpliedVol
**	Summary:		return the implied vol of an option
**	References:		
**	Description:	
**	Example:		
****************************************************************/

double MlEqBSImpliedVol(
						double option_val,		// option value
						double forward,			// forward price
						double maturity,		// maturity
						double strike,			// strike
						double discount_factor,	// discount factor
						int cp,					// 1 : call, -1 : put, 0 : forward
						int rootfind_flag ,		// rootfinder flag 
						double accuracy   ,		// accuracy 
						double lower_x    ,			// lower bound to search 
						double upper_x    			// upper bound to search 
)
{	
	double result;
	int max_tries = 50;
	
	if ( cp != 1 && cp != -1 && cp != 2 )
		return 0;

	DataContainer impl_data;
	
	impl_data.stock_forward = forward;
	impl_data.mat = maturity;
	impl_data.discount_factor = discount_factor;
	impl_data.cp = cp;
	impl_data.strike = strike;
	
	
	rootfind_solve( rootfind_flag, implied_vol_fcn,  lower_x,  upper_x,
				 accuracy,  max_tries,  option_val, &impl_data, &result,NULL);
	
	return result;
}	



inline double Hermit(int n, double x, double y )
{
	double res;
	if ( n < 0 ){
		throw( "index must be greater or equal to zero in Hermit polynomial calculation");
	}
	if ( n == 0 ){
		return 1.0;
	}
	if ( n == 1 ){
		res = x;
		return x;
	}
	if ( n == 2 ){
		res = x*x-y;
		return res;
	}
	if ( n == 3 ){
		res = pow(x,3.0)-3.0*x*y;
		return res;
	}
	if ( n == 4 ){
		res = pow(x,4.0)-6.0*x*x*y+3.0*y*y;
		return res;
	}
	res = x*Hermit(n-1,x,y)-(n-1.0)*y*Hermit(n-2,x,y);
	return res;
}



class  hermitContainer
{
	public:

	CVector coeff;
	double mat;
	double normalization;
	CVector strikes;
	double fwd;
	double sqrtMat;
	double sqrt_2pi;
	double elasticity;
	int ngauss;
	double volvol;
	CVector gaussWeight;
	CVector gaussPoint;
	int nstates;
};




/****************************************************************
**	class:			hermiteSkew
**	Summary:		stochHermitNormalization
**	References:		
**	Description:	
**	Example:		
****************************************************************/


void stochHermitNormalization(double x, double z[], void *vp)
{

	hermitContainer* hc = (hermitContainer*) vp;

	CVector& coeff = hc->coeff;
	double mat = hc->mat;

	int i;

	double res = 0.0;
	double volvol = hc->volvol;

	double stateVar = 0.0;
	double stochFac=1,stochProb=1,scale=1.0,xscale;
	for ( int istate = 0 ; istate < hc->nstates; istate++ )
	{
		if ( hc->ngauss > 1 )
		{
			stochFac	=	exp(-0.5*volvol*volvol*mat+volvol*hc->sqrtMat*hc->gaussPoint[istate]);
			stochProb	=	hc->gaussWeight[istate]/hc->sqrt_2pi*exp(-0.5*hc->gaussPoint[istate]*hc->gaussPoint[istate]);	
			scale		=   pow(stochFac,hc->elasticity);
		}

		xscale = 1.0;

		stateVar = 0.0;
		for ( i = 0 ; i < coeff.getsize(); i++ )
		{
			if( i == 0 ){
				xscale = stochFac;
			}else{
				xscale = scale;
			}

			stateVar += coeff[i]*xscale*Hermit(i+1, x*hc->sqrtMat, mat );
		}

		res += stochProb*exp(stateVar);

	}

	z[1] = res*exp(-0.5*x*x)/hc->sqrt_2pi;

}

/****************************************************************
**	class:			none
**	Summary:		stochHermitOptions
**	References:		
**	Description:	
**	Example:		
****************************************************************/


void stochHermitOptions(double x, double z[], void *vp)
{


	int i;

	hermitContainer* hc = (hermitContainer*) vp;

	CVector& coeff			= hc->coeff;
	double mat				= hc->mat;
	CVector& strikes		= hc->strikes;
	double normalization	= hc->normalization;
	double fwd				= hc->fwd;
	double volvol			= hc->volvol;
	double vol				= coeff[0];


	for (i = 0 ; i < strikes.getsize(); i++ ){
			z[i+1] = 0.0;
	}

	double stateVar = 0.0;

	double res = 0.0;
	double stochFac=1,stochProb=1,scale=1.0,xscale;
	for ( int istate = 0 ; istate < hc->nstates; istate++ )
	{
		if ( hc->ngauss > 1 )
		{
			stochFac	=	exp(-0.5*volvol*volvol*mat+volvol*hc->sqrtMat*hc->gaussPoint[istate]);
			stochProb	=	hc->gaussWeight[istate]/hc->sqrt_2pi*exp(-0.5*hc->gaussPoint[istate]*hc->gaussPoint[istate]);	
			scale		=   pow(stochFac,hc->elasticity);
		}

		xscale = 1.0;
		stateVar = 0.0;
		for ( i = 0 ; i < coeff.getsize(); i++ )
		{
			if( i == 0 ){
				xscale = stochFac;
			}else{
				xscale = scale;
			}
			stateVar += coeff[i]*xscale*Hermit(i+1, x*hc->sqrtMat, mat );
		}

		double spot = fwd*normalization*exp(stateVar);
		for (i = 0 ; i < strikes.getsize(); i++ ){
			z[i+1] += stochProb*MlEqMaths::Max(spot-strikes[i],0.0);
		}

	}

	double prob = exp(-0.5*x*x)/hc->sqrt_2pi;
	for (i = 0 ; i < strikes.getsize(); i++ ){
			z[i+1] *= prob;
	}

}


/****************************************************************
**	class:			none
**	Summary:		hermitNormalization
**	References:		
**	Description:	
**	Example:		
****************************************************************/

void hermitNormalization(double x, double z[], void *vp)
{

	hermiteSkew* hc = (hermiteSkew*) vp;

	CVector& coeff = hc->m_coeff;
	double mat = hc->m_maturity;

	double stateVar		 = 0.0;
	double stateVariance = 0.0;
	double n = mat;
	for (int i = 0 ; i < hc->m_coeff.getsize(); i++ ){
		stateVar  += coeff[i]*Hermit(i+1,x*hc->m_sqrtMat, mat );

	}

	z[1] = exp(-0.5*x*x)/hc->m_sqrt_2pi*exp(stateVar);

}



/****************************************************************
**	class:			none
**	Summary:		hermitOptions
**	References:		
**	Description:	
**	Example:		
****************************************************************/

void hermitOptions(double x, double z[], void *vp)
{
	int i;

//	hermitContainer* hc = (hermitContainer*) vp;
	hermiteSkew* hc = (hermiteSkew*) vp;

	CVector& coeff			= hc->m_coeff;
	double mat				= hc->m_maturity;
	CVector& strikes		= hc->m_strikes;
	double normalization	= hc->m_normalization;
	double fwd				= hc->m_forward;

	double stateVar		 = 0.0;
	double stateVariance = 0.0;
	double n = mat;
	for ( i = 0 ; i < coeff.getsize(); i++ ){
		stateVar  += coeff[i]*Hermit(i+1,x*hc->m_sqrtMat, mat );
	}

	double spot = fwd*normalization*exp(stateVar);
	double prob = exp(-0.5*x*x)/hc->m_sqrt_2pi;

	for (i = 0 ; i < strikes.getsize(); i++ ){
		z[i+1] = prob*MlEqMaths::Max(spot-strikes[i],0.0);
	}

}


/****************************************************************
**	class:			hermiteSkew
**	Summary:		initialize
**	References:		
**	Description:	
**	Example:		
****************************************************************/


void hermiteSkew::initialize(double maturity,double forward,const CVector& hermitCoeff,bool adjustATM)
{

	m_maturity		= maturity;
	m_forward		= forward;
	m_hermitCoeff	= hermitCoeff;
	m_sqrtMat		= sqrt(maturity);
	m_sqrt_2pi		= sqrt(2.0*3.141592654);

	if ( maturity < 1e-3 ){
		throw("maturity zero entered");
	}

	m_coeff.resize(hermitCoeff.getsize());
	double n = 1.0;
	for (int  i = 0 ; i < m_coeff.getsize(); i++ )
	{
		m_coeff[i] = hermitCoeff[i]/(n*pow(maturity,(double)i/2.0));
		n*= (double)(i+2);
	}
	m_coeff[0] *= m_factor;

	double eps=1e-8,hInit=0.1,hMin=1e-2,bump=0.0008,scale = 10;

	CVector integrationLimits(2);
	integrationLimits[0] = -5.0;
	integrationLimits[1] = 5.0;

//  calculate martingale compensator first



/*	int N = 15;

	CVector gaussWeight(N);
	CVector gaussPoint(N);

	double mat = maturity;
			
	MlEqMaths::dGauleg(integrationLimits[0],integrationLimits[1],gaussPoint,gaussWeight,N,true);

	m_normalization = 0.0;
	for ( int n = 0 ; n < N ; n++ )
	{
		double stateVar	= 0.0;
		double x = gaussPoint[n];

		for (int i = 0 ; i < m_coeff.getsize(); i++ ){
			stateVar  += m_coeff[i]*Hermit(i+1,x*m_sqrtMat, mat );
		}

		m_normalization += gaussWeight[n]*exp(-0.5*x*x)/m_sqrt_2pi*exp(stateVar);
	}
*/


	Runge_Kutta_Integrate_new(m_normalization,hermitNormalization,integrationLimits,
							  hMin,eps,bump,this,hInit,scale);



	m_normalization = 1.0/m_normalization;
	m_logNormalization = log(m_normalization);


	int maxiter = 4;

	if ( !adjustATM || m_iter > maxiter){
		return;
	}

	m_iter++;

	CVector stk(1),res(1);
	stk[0] = m_forward;
	calculateOptions(res, stk,true);

	if ( fabs(res[0] - hermitCoeff[0]) < 0.0025 ){
			maxiter = maxiter;
			return;
	}
	else{
			m_factor *= m_hermitCoeff[0]/res[0];
			initialize(maturity,forward,hermitCoeff,adjustATM);
	}

}

/****************************************************************
**	class:			hermiteSkew
**	Summary:		adjustSlope
**	References:		
**	Description:	
**	Example:		
****************************************************************/

double hermiteSkew::newtonUpdate(double hermit1,double forward)
{
	CVector res(2),strikes(2);
	double bump = 0.05;

	double ns = 0.01;

	strikes[0] = forward*exp(ns*m_hermitCoeff[0]*sqrt(m_maturity));
	strikes[1] = forward*exp(-ns*m_hermitCoeff[0]*sqrt(m_maturity));

	m_hermitCoeff[1] = hermit1;
	m_hermitCoeff[2] = -m_hermitCoeff[1];
	m_hermitCoeff[3] = -m_hermitCoeff[2]/2.0;	

	initialize(m_maturity,forward,m_hermitCoeff);
	calculateOptions(res, strikes);
	double xslope = (res[0]-res[1])/(2.0*ns*m_hermitCoeff[0]);

	return xslope;
}

/****************************************************************
**	class:			hermiteSkew
**	Summary:		adjustSlope
**	References:		
**	Description:	
**	Example:		
****************************************************************/

double hermiteSkew::getStateVar(double x )
{
	double stateVar		 = 0.0;
	for (int i = 0 ; i < m_coeff.getsize(); i++ ){
		stateVar  += m_coeff[i]*Hermit(i+1,x*m_sqrtMat, m_maturity );
	}

	stateVar += m_logNormalization;
	return stateVar;
}


/****************************************************************
**	class:			hermiteSkew
**	Summary:		adjustSlope
**	References:		
**	Description:	
**	Example:		
****************************************************************/

double hermiteSkew::adjustSlope(double targetSlope,double forward,double eps)
{


	double new_x,deriv,res;
	
	double val	= newtonUpdate(m_hermitCoeff[1],forward);
	deriv		= val/m_hermitCoeff[1];

	double old_x = m_hermitCoeff[1];
	double old_val = val;

	int max_iter = 7;

	double absoluteEps = eps*0.02;

	int n;
	for ( n = 0 ; n < max_iter; n++ )
	{	
		new_x = -(val-targetSlope)/deriv+old_x;
		
		old_val = val;
		val	= newtonUpdate(new_x,forward);

		if ( fabs(val-targetSlope) < absoluteEps )
		{
		  res = new_x;
		  return res;
		}	
		if ( fabs(val/targetSlope-1.0) < eps )
		{		
		  res = new_x;
		  return res;
		}

		deriv = (val-old_val)/(new_x-old_x);
		old_x = new_x;		

	}
	
	if ( n == max_iter ){
		throw("rootfind did not find solution" );
	}

	return -1e99;
}


/****************************************************************
**	class:			hermiteSkew
**	Summary:		calculateOptions
**	References:		
**	Description:	
**	Example:		
****************************************************************/

void hermiteSkew::calculateOptions(CVector& result, CVector& strikes,bool returnVolFlag)
{
	result.resize(strikes.getsize());

	double eps=1e-8,hInit=0.1,hMin=1e-2,bump=0.0008,scale = 10;

	CVector integrationLimits(2);
	integrationLimits[0] = -5.0;
	integrationLimits[1] = 5.0;
	
	m_strikes = strikes;

	Runge_Kutta_Integrate_new_vector(result, hermitOptions, integrationLimits,hMin,eps,bump,this, hInit,scale);

//  calculate implied vols

	if ( !returnVolFlag ){
		return;
	}

	double accuracy   = 1e-10;
	double lower_x    = 0.03; 
	double upper_x    = 1.5; 

	int rootfind_flag = eROOTFIND_BRENT_NOGROW;

	CVector impliedVols(strikes.getsize());
	for (int i = 0; i < strikes.getsize(); i++ )
	{


		impliedVols[i] =	MlEqBSImpliedVol(
							result[i],
							m_forward,
							m_maturity,
							strikes[i],
							1,1,rootfind_flag ,	 
							accuracy   ,		 
							lower_x    ,		 
							upper_x);
		
	}

	result = impliedVols;
}


void HermiteOptions(CVector& result, CVector& strikes,const  CVector& hermitCoeff,double maturity,double forward,int returnVolFlag,double elasticity,double volvol,int ngauss)
{

	if ( volvol < 1e-2 || ngauss == 0 )
	{
		bool adjustATM = true;
		hermiteSkew hSk;
		hSk.initialize( maturity, forward, hermitCoeff,adjustATM);
		hSk.calculateOptions(result, strikes, returnVolFlag ? true : false);

		return;
	}


	CVector coeff(hermitCoeff.getsize());

	if ( maturity < 1e-3 ){
		throw("maturity zero entered");
	}

	double n = 1.0;
	for (int  i = 0 ; i < coeff.getsize(); i++ )
	{
		coeff[i] = hermitCoeff[i]/(n*pow(maturity,(double)i/2.0));
		n*= (double)(i+2);
	}

//  estimate at-the money vol

	hermitContainer hC;

	hC.coeff = coeff;
	hC.mat   = maturity;
	hC.fwd = forward;
	hC.sqrtMat = sqrt(maturity);
	hC.sqrt_2pi = sqrt(2.0*3.141592654);	
	hC.elasticity = elasticity;
	hC.ngauss  = ngauss;
	hC.volvol = volvol;

	CVector gaussWeight,gaussPoint;

	if ( volvol < 1e-2 ){
		ngauss = 0;
	}

	if ( ngauss )
	{
		gaussPoint.resize(ngauss);
		gaussWeight.resize(ngauss);

		double a = -4;
		double b =  4;
		MlEqMaths::dGauleg(a, b,gaussPoint,gaussWeight,ngauss,true);
	}

	hC.gaussWeight = gaussWeight;
	hC.gaussPoint  = gaussPoint;

	if ( ngauss ){
		hC.nstates = ngauss;
	}else{
		hC.nstates = 1;
	}


	double eps=1e-8,hInit=0.1,hMin=1e-2,bump=0.0008,scale = 10;

	CVector integrationLimits(2);
	integrationLimits[0] = -5.0;
	integrationLimits[1] = 5.0;

//  calculate martingale compensator first

	double normalization;

	if ( ngauss == 0 )
	{
		Runge_Kutta_Integrate_new(normalization,hermitNormalization,integrationLimits,
								  hMin,eps,bump,&hC,hInit,scale);
	}
	else
	{
		Runge_Kutta_Integrate_new(normalization,stochHermitNormalization,integrationLimits,
								  hMin,eps,bump,&hC,hInit,scale);
	}


	hC.normalization = 1.0/normalization;
	hC.strikes = strikes;

//  start option pricing

	result.resize(strikes.getsize());

	if ( ngauss == 0 )
	{
		Runge_Kutta_Integrate_new_vector(result, hermitOptions, integrationLimits,hMin,eps,bump,&hC, hInit,scale);
	}
	else
	{
		Runge_Kutta_Integrate_new_vector(result, stochHermitOptions, integrationLimits,hMin,eps,bump,&hC, hInit,scale);
	}


//  calculate implied vols

	if ( !returnVolFlag ){
		return;
	}


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










void STLVectorFromCVector(vector < double > & stlVec,const  CVector & in)
{
	stlVec.resize(in.getsize());
	for ( int i = 0 ; i < in.getsize(); i++ )
	{
		stlVec[i] = in[i];
	}
}

void CVectorFromSTLVector(CVector& vec,const vector < double > & in)
{
	vec.resize(in.size());
	for ( int i = 0 ; i < in.size(); i++ )
	{
		vec[i] = in[i];
	}

}

void STLVectorVectorFromCMatrix(vector < vector < double > > & stlVec,const CMatrix& in )
{
	stlVec.resize(in.rows());
	int i,j;


	for ( i = 0 ; i < in.rows(); i++ )
	{
		stlVec[i].resize(in.cols());
		for ( j = 0 ; j < in.cols(); j++ )
		{
			stlVec[i][j] = in[i][j];
		}
	}
}


void STLVectorVectorTransposeFromCMatrix(vector < vector < double > > & stlVec,const CMatrix& in )
{
	stlVec.resize(in.cols());
	int i,j;

	for ( i = 0 ; i < in.cols(); i++ )
	{
		stlVec[i].resize(in.rows());
		for ( j = 0 ; j < in.rows(); j++ )
		{
			stlVec[i][j] = in[j][i];
		}
	}
}


void CMatrixFromSTLVectorVector(CMatrix& mat, const vector < vector < double > > & in )
{
	if ( in.size() == 0 )
		return;

	mat.resize(in.size(),in[0].size());
	int i,j;
	for ( i = 0 ; i < mat.rows(); i++ )
	{
		for ( j = 0 ; j < mat.cols(); j++ )
		{
			mat[i][j] = in[i][j];
		}
	}
}


/****************************************************************
**	Class  : product 
**	Routine: accruePayout
**	Returns: nothing
**	Action : accrues payoffs over all simulation paths
**  
****************************************************************/


void product::accruePayout(CMatrix& result,CMatrix& value)
{
 	    int i,j;
		for ( i = 0 ; i < result.rows(); i++ )
		  for ( j = 0 ; j < result.cols(); j++ )
				result[i][j] += value[i][j];	
}



/****************************************************************
**	Class  : product 
**	Routine: createResult
**	Returns: nothing
**	Action : averages over paths and calculates errors
**  
****************************************************************/


void product::createResult(CMatrix& result,MlEqMonteCarlo& mc)
{
	int i;
	int npaths = mc.m_nPaths;

	for ( i = 0 ; i < result.rows(); i++ )
	{
//			result[i][0] /= (double) npaths;
			mc.averageMonteCarloResults(result[i]);

/*
//			result[i][1] = result[i][1]/(double) npaths-pow(result[i][0],2.0);
			result[i][1] -= pow(result[i][0],2.0);

			if ( result[i][1] <  0.0 ){
				if ( fabs(result[i][1]/result[i][0]) < 1e-5 ){
					result[i][1] = 0.0;}};

			result[i][1] /= (double) npaths;		
			result[i][1] = 2.0*sqrt(result[i][1]);		
*/
	}		
}

/****************************************************************
**	Class  : product 
**	Routine: fill
**	Returns: nothing
**	Action : calculates second moment of price (for error calculation)
**  
****************************************************************/



void product::fillValues(CMatrix& value,CMatrix& path_array,MlEqMonteCarlo& m)
{
	for (int i = 0 ; i < value.rows(); i++ )
	{
		value[i][1] = value[i][0]*value[i][0];
	}
	return;
}

/****************************************************************
**	Class  : product 
**	Routine: setUp
**	Returns: nothing
**	Action : set up
**  
****************************************************************/

void product::setUp(CMatrix& value,MlEqMonteCarlo& mc)
{
//	mc.initializeMcHelper(*this);
}

/****************************************************************
**	Class  : MlEqDayCounter 
**	Routine: MlEqDayCounter
**	Returns: nothing
**	Action : set up
**  
****************************************************************/

/*
MlEqDayCounter::MlEqDayCounter(DayCountConvention_ dcType )
:DateToDouble(dcType)
{}
*/


double getObjectFromCMatrix(CMatrix& objects,int iRow,int iCol,int numberOfSlices)
{

	int nsize = objects.size();
	if ( nsize == 1 )
	{
		if ( objects.cols() <= iCol ){
			throw("error indexing Column array");
		}

		return objects[0][iCol];
	}
	else if ( nsize == numberOfSlices )
	{
		if ( objects.cols() <= iCol ){
			throw("error indexing strikearray");
		}

		return objects[iRow][iCol];
	}
	else{
		throw("indexing error ");
	}
}


double  getObjectFromCVector(CVector&  objects,int index)
{
	int nsize = objects.getsize();
	if ( nsize == 1 )
	{
		if ( objects.getsize() <= index ){
			throw("error indexing Column array");
		}

		return objects[0];
	}
	else 
	{
		if ( objects.getsize() <= index ){
			throw("error indexing strikearray");
		}

		return objects[index];
	}

};









