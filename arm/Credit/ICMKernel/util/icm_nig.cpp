#include "ICMKernel\util\icm_nig.h"
#include "gpclosedforms\bessel.h"

#include "ICMKernel/glob/icm_maths.h"
#include <nagd01.h>

using namespace ARM ;

ICM_Nig::ICM_Nig (double alpha, double BetaNig, double Mu, double Delta, double Rho)
{
	setParameters(alpha, BetaNig, Mu, Delta, Rho );
}

ICM_Nig::ICM_Nig (double	alpha, double BetaNig, double Rho)
{
	setParameters(alpha,
				  BetaNig,
				  - alpha * BetaNig/sqrt(alpha*alpha-BetaNig*BetaNig),
				  alpha,
				  Rho);
}

void ICM_Nig::setParameters(double Alpha, double BetaNig, double Mu, double Delta, double Rho)
{
	Set_NIG_Alpha(Alpha) ;
	Set_NIG_Beta(BetaNig) ;
	Set_NIG_Mu(Mu) ;
	Set_NIG_Delta(Delta) ;
	SetRho(Rho) ;
}

double ICM_Nig::generateRandom()
{
	// random NIG variable
	double itsAlpha, itsBetaNig, itsMu, itsDelta ;
	
	itsAlpha = Get_NIG_Alpha() ;
	itsBetaNig = Get_NIG_Beta() ;
	itsMu = Get_NIG_Mu();
	itsDelta = Get_NIG_Delta();


	double normal_1 = NAG_random_normal(0.,1.);
	double V = normal_1*normal_1 ;
	double ksi = itsDelta*itsDelta ;
	double psi = itsAlpha*itsAlpha - itsBetaNig*itsBetaNig ;
	double vega = sqrt(ksi/psi) ;
	
	double y1 = vega + vega *(vega*V-sqrt(4 * vega * ksi * V+vega*vega*V*V))/(2*ksi) ;
	double y2 =vega*vega/y1 ;

	double proba = vega/(vega+y1) ;

	//Bernouilli test with proba
	double Y ;

	if(NAG_random_continuous_uniform()<proba)
		Y = y1 ;
	else
		Y = y2 ;

	double normal_2 = NAG_random_normal(0.,1.) ;

	double X = itsMu + itsBetaNig * Y + sqrt(Y)*normal_2 ;
	return X ;
}

double ICM_Nig::GetDensity(double x)
{
	double itsAlpha, itsBetaNig, itsMu, itsDelta ;
	
	itsAlpha = Get_NIG_Alpha() ;
	itsBetaNig = Get_NIG_Beta() ;
	itsMu = Get_NIG_Mu();
	itsDelta = Get_NIG_Delta();

	double gamma = sqrt(itsAlpha*itsAlpha-itsBetaNig*itsBetaNig);

	double aux = bessel_fractional_K(1., itsAlpha*sqrt(itsDelta*itsDelta+(x-itsMu)*(x-itsMu)));
	
	double res = aux * itsDelta*itsAlpha*
					exp(itsDelta*gamma+itsBetaNig*(x-itsMu))/
					(PI*sqrt(itsDelta*itsDelta+(x-itsMu)*(x-itsMu))) ;
		
	return res ;	
}

static double __stdcall FunctionToIntegrate (double y,Nag_User *comm)
{
	double itsAlpha, itsBetaNig, itsMu, itsDelta ;
	
	ICM_Nig* ptr = (ICM_Nig*) comm->p;
	
	itsAlpha = ptr->Get_NIG_Alpha();
	itsBetaNig = ptr->Get_NIG_Beta();
	itsMu = ptr->Get_NIG_Mu();
	itsDelta = ptr->Get_NIG_Delta();

	double gamma = sqrt(itsAlpha*itsAlpha-itsBetaNig*itsBetaNig);

	double x = ptr->GetitsX(); 

	double aux = (x-(itsMu+itsBetaNig*y))/sqrt(y) ;
	
	double aux1 = cdfNormal(aux);
	
	return aux1 * (gamma*itsDelta) * 
		   exp(-(itsDelta*gamma-gamma*gamma*y)*(itsDelta*gamma-gamma*gamma*y)/(2*gamma*gamma*y))/
		   (sqrt(2*PI)*gamma*pow(y,1.5)) ;
}

//----------------------------------------------------------------------------
// Computing NIG CDF
//----------------------------------------------------------------------------

double ICM_Nig::GetDistribution(double x) 
{
	double resultat = 0. ;
		
	SetX(x) ;

	resultat = IntegrateByGaussLaguerre(0.0, 0.5, FunctionToIntegrate);

	return resultat;
}

//----------------------------------------------------------------------------
// Computing CDF Inverse by dichotomy (to be Optimized)
//----------------------------------------------------------------------------

double ICM_Nig::NIG_Inverse(double y)
{
	// By dichotomy
	double x_inf, x_sup, y_inf, y_sup ;
	x_inf = -10000 ;
    x_sup = 10000 ;

	y_inf = GetDistribution(x_inf);
	y_sup = GetDistribution(x_sup);

	int count = 0 ;
	double x, res ;
	double aux = fabs((y_sup+y_inf)/2. - y );

	while (aux >=0.000000000001 && count <= 100) 
	{
		x = (x_sup+x_inf)/2 ;
		res = GetDistribution(x) ;

		if (res < y)
		{
			x_inf = x ;
			y_inf = GetDistribution(x_inf) ;
		}
		else
		{
			x_sup = x ;
			y_sup = GetDistribution(x_sup) ;
		}

		count++;
	}

	x = (x_sup+x_inf)/2 ;
		
	return x;
}

//----------------------------------------------------------------------------
// Computing Implied Default Barrier
//----------------------------------------------------------------------------

double ICM_Nig::GetDefaultBarrier (double MktDefProb)
{
	double factor = 1.0/GetitsRho();

	ICM_Nig* nigcopula = new ICM_Nig(factor*Get_NIG_Alpha(),
									 factor*Get_NIG_Beta(),
									 factor*Get_NIG_Mu(),
									 factor*Get_NIG_Delta(),
									 factor);

	double Barrier = nigcopula->NIG_Inverse(MktDefProb);

	delete nigcopula ;

	return Barrier;
}

//----------------------------------------------------------------------------
// Computing Conditional Default Probability
//-----------------------------------------------------------------------------

double ICM_Nig::getConditionalDefaultProbability(double MktLevel, double MktDefProb)
{
	double DefaultBarrier = GetDefaultBarrier(MktDefProb);

	double rho = GetitsRho();

	double factor = sqrt((1-rho*rho)/rho);
	
	ICM_Nig* nigcopula = new ICM_Nig(factor*Get_NIG_Alpha(),factor*Get_NIG_Beta(),factor*Get_NIG_Mu(),factor*Get_NIG_Delta(), factor);
	
	double Proba = nigcopula->GetDistribution((DefaultBarrier-rho*MktLevel)/sqrt(1-rho*rho)) ;

	delete nigcopula ;

	return Proba;
}

//----------------------------------------------------------------------------
// Gauss Laguerre Integration by NAG
//-----------------------------------------------------------------------------

double ICM_Nig::IntegrateByGaussLaguerre(double a, double b, double(__stdcall fun)(double x,Nag_User *comm))
{
	/*
	This fucntion computes an estimate of the definite integral of a 
	function of known analytical form, using a Gaussian quadrature formula
	with a specified number of abscissae. The Formula is provided for a 
	semi-infinite interval (Gauss-Laguerre)
	*/
	
	static Integer nstor[1] = {64};
	
	NagError fail;
	INIT_FAIL(fail); 
	double ans;
	
	Nag_GaussFormulae gaussformula;
	Boolean success = TRUE;
	// fail.print = TRUE;
	
	gaussformula = Nag_Laguerre;
	
	Nag_User comm;
	comm.p = this ;
	ans = d01tac(gaussformula, fun, a, b, nstor[0],&comm, &fail);

	if (fail.code == NE_NOERROR || fail.code == NE_QUAD_GAUSS_NPTS_RULE)
		success= TRUE;
	else
	{
		success = FALSE;
	}
	
	if (success)
		return ans;
	else
		ICMTHROW(ERR_INVALID_ARGUMENT,"exit(EXIT_FAILURE);"); 
}

