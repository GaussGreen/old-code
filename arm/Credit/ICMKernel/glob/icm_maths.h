/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_MATHS.H
	PROJECT:	GLOB
	
	DESCRIPTION:	some useful math tools

  -----------------------------------------------------------------

 	CREATION:	November 5, 2004

	LAST MODIF:	November 5, 2004
  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

# ifndef _MATHS
# define _MATHS

#include <nag.h>		
#include <nags.h>		
#include <nagg01.h>		
#include <nagg05.h>		
#include "ICMKernel\util\icm_integrator.h"
#include <vector>

typedef ARM_Vector Array;
typedef unsigned int Size;

// 17783 using namespace std;

# define FMIN(a,b) (((a)<(b))?(a):(b)) 
# define FMAX(a,b) (((a)>(b))?(a):(b)) 
# define SQR(a) ((a)*(a))

// replicates Sgn as in visual basic, the signum of a real number
# define VB_SGN(a) ((a) > (0) ? (1) : (((a) < 0) ? (-1) : (0)))


// functions for dates
# define MATHTIME(a)		((a) * ONEOVER365)
# define DAYSTIME(a)		((a) * 365)

// some constants
# define MY_PI				3.14159265358979
# define ONEOVER365			2.73972602739726027397260273972e-3
# define SQRTTWOPI			2.506628274631000502415765284811
# define ONEOVERSQRTTWOPI	0.39894228040143267793994605993438
# define ONEOVERSQRTTWOPI	0.39894228040143267793994605993438
# define SQRT_TWO			1.4142135623730950488016887242097
# define TWO_PI				6.283185307179586476925286766559

# define _PLUS_INFINITY_	1e25
# define _MINUS_INFINITY_	-1e25

//# define _BAD_RET_VALUE_	1e99

// -------------------------------------------------------------
// ICM INTEGRATOR

# define	ZEPS	1e-6
# define	EPS_GH	3.0e-14				// Relative precision.
# define	PIM4_GH 0.7511255444649425	//	1/Pi^1/4.
# define	MAXIT_GH 10					//	Maximum iterations.

// -------------------------------------------------------------

// BLACK SHOLES -- TOOL BOX

# define CalBs(g, s, k, v, t) ((log((s) / (k)) + (g) * (v) * (v) * (t) * 0.5) / (v * sqrt(t))) 

extern double BSFormula(
					bool	IsCall,			// true = cap, false = floor
					double	Sigma,
					double	Spot,
					double	Strike,
					double	NbD,
					bool	Hedges,			// true if one wants all greeks to be computed
					double* Price,
					double* Delta,
					double* Gamma,
					double* Vega,
					double* Teta);

double	 BSNnValues(
					bool	Call,			// true = cap, false = floor
					double	Sigma,
					double	Spot,
					double	Strike,
					double	NbD,
					double* Nd1,
					double* Nd2,
					double* nd1,
					double* nd2);

extern double	NormalFormula(
					bool	IsCall,			// true = cap, false = floor
					double	Sigma,
					double	Spot,
					double	Strike,
					double	NbD,
					bool	Hedges,			// true if one wants all greeks to be computed
					double* Price,
					double* Delta,
					double* Gamma,
					double* Vega,
					double* Teta);

double StandardGaussianDensity(double x);


// --------------------------------------------------------------------------------------
// BIVARIATE
// http://www.mathfinance.de/FF/cpplib.html
// --------------------------------------------------------------------------------------

double	ndf(double t);
double	nc(double x);
double	fxy(double x, double y, double a, double b, double rho);
double	Ntwo(double a, double b, double rho);
double	ND2(double a, double b, double rho);


// density, Cumulative Prob and it's Inverse of the Student Distribution
double StudentDensity(int NStudent, double x);

double StudentCum(int NStudent, double x);

double InvStudentCum(int NStudent, double p);

double Chi2Density (int N, double x); 

double gammln(double xx);
// Returns the value ln|gamma(xx)| for xx > 0.

double beta(double z, double w);
// Returns the value of the beta function B(z; w).

double betai(double a, double b, double x);
// Returns the incomplete beta function I x (a; b).

double betacf(double a, double b, double x);
// Used by betai: Evaluates continued fraction for incomplete beta function by modied Lentz's method ( x 5.2).

// NAG Functions


//	----------------------------------------------------------------------------------------
//	g01fac
static inline double NAG_deviates_normal(double p)
{
	NagError fail; 
	INIT_FAIL(fail) ;
	double ret = g01fac(Nag_LowerTail, p,&fail);
	// if (fail.code==NE_NOERROR) 
	// 	ICMLOG("NAG: "<<fail.message<<"["<<fail.errnum<<"]"); 
	return ret ;
}
//	----------------------------------------------------------------------------------------
// g01ecc
static inline double NAG_deviates_normal_dist(double p)
{
	NagError fail; 
	INIT_FAIL(fail) ;
	double ret=nag_deviates_normal_dist(p,&fail);
	// if (fail.code==NE_NOERROR) 
	// 	ICMLOG("NAG: "<<fail.message<<"["<<fail.errnum<<"]"); 
	return ret ;
}
//	----------------------------------------------------------------------------------------
//	g01ebc
static inline double NAG_prob_students_t(double t,double df)
{
	NagError fail; 
	INIT_FAIL(fail) ;
	double ret= g01ebc(Nag_LowerTail, t,df,&fail);
	// if (fail.code==NE_NOERROR) 
	// 	ICMLOG("NAG: "<<fail.message<<"["<<fail.errnum<<"]"); 
	return ret ;

}
//	----------------------------------------------------------------------------------------
// g01fbc
static inline double NAG_deviates_students_t(double t,double df)
{
	NagError fail; 
	INIT_FAIL(fail) ;
	double ret= g01fbc(Nag_LowerTail, t,df,&fail);
	// if (fail.code==NE_NOERROR) 
	// 	ICMLOG("NAG: "<<fail.message<<"["<<fail.errnum<<"]"); 
	return ret ;
}
//	----------------------------------------------------------------------------------------
// g01hac
static inline double NAG_bivariate_normal_dist(double x,double y,double rho)
{
	NagError fail; 
	INIT_FAIL(fail) ;
	double ret=  g01hac(x,y,rho,&fail);
	// if (fail.code==NE_NOERROR) 
	// 	ICMLOG("NAG: "<<fail.message<<"["<<fail.errnum<<"]"); 
	return ret ;

}

static inline void NAG_poisson_dist(double rlambda, int k, double& probainfequalk, 
									  double& probasupk, double& probaequalk)
{
	NagError fail; 
	INIT_FAIL(fail) ;
	g01bkc(rlambda,k,&probainfequalk,&probasupk,&probaequalk,&fail);
	// if (fail.code==NE_NOERROR) 
	// 	ICMLOG("NAG: "<<fail.message<<"["<<fail.errnum<<"]"); 
	//return ret ;
}

//	----------------------------------------------------------------------------------------
// s15abc
static inline double NAG_cumul_normal(double x)
{
	double ret=  s15abc(x);
	return ret ;
}
//	----------------------------------------------------------------------------------------
// g05eac
// static inline void NAG_ref_vec_multi_normal(double *a,int n, double *c, int tdc,double eps,double**r)
// {
// 	NagError fail; 
// 	INIT_FAIL(fail); 
// 	g05eac(a,n,c,tdc,eps,r,&fail); 
// 	// if (fail.code==NE_NOERROR) 
// 	// 	ICMLOG("NAG: "<<fail.message<<"["<<fail.errnum<<"]"); 
// }
//	----------------------------------------------------------------------------------------
// g05cbc
static inline void NAG_random_init_repeatable(int seed)
{
	g05cbc(seed) ;
}
//	----------------------------------------------------------------------------------------
// // g05ddc
static inline double NAG_random_normal(double a,double b)
{
	return g05ddc(a,b) ;
}
//	----------------------------------------------------------------------------------------
// g05dac
static inline double NAG_random_continuous_uniform_ab(double a,double b)
{
	return g05dac(a,b); 
}
//	----------------------------------------------------------------------------------------
// g05cac
static inline double NAG_random_continuous_uniform()
{
	return g05cac() ;
}
// 
// g05ezc
// void	NAG_return_multi_normal(double* z,double*r); 

double SplineInterpol(vector<double>& x,vector<double>& y,double value,double smooth,vector<double>& weights);
double HermiteInterpol(vector<double>& x,vector<double>& y,double value);

double DigitalPriceBS(double& yf,double& spot,double& strike,double& vol,double& vol_epsilon,double& epsilon);

double CorridorPriceBS(double& yf,
					   double& spot,
					   double& K1,
					   double& K2,
					   double& volK1,
					   double& vol_epsilonk1,
					   double& volK2,
					   double& vol_epsilonk2,
					   double& epsilon);

// void GaussLegendre_ComputeAbscissasAndWeights(FILE* fOut,const double& x1,const double& x2,const int& n);
// void GaussHermite_Display();
// void GaussHermite_ComputeAbscissasAndWeights(FILE* fOut, int n);
// void GaussLegendre_Display();

// ---------------------------------------------
// computation of MAX(X1+X2+...+Xn - K,0.)
// ---------------------------------------------
class LinearLognProcess{
public :
	bool	itsIsComputed;
	double  itsCallPrice;

	int		its_intstep;
	double	its_yt;
	double	its_strike;
	double	its_rho;
	vector<double> its_coef;
	vector<double> its_spot;			  
	vector<double> its_vol;

	ICM_Integrator	TheIntegrator;

	double  its_IntegrationValue;

public :
	LinearLognProcess(const double& yt,
		const double& strike,
		const double& correlation,
		const vector<double>& coef,			  
		const vector<double>& spot,			  
		const vector<double>& vol,
		int intstep = 20 ):its_yt(yt),its_coef(coef),its_spot(spot),
		its_strike(strike),its_vol(vol),its_intstep(intstep),its_rho(correlation)
	{}

	void Init(){
		its_yt=0.;
		its_strike=0.;
		itsCallPrice=0.;
		itsIsComputed=false;
		its_intstep=20;
		its_rho=0.;
		its_IntegrationValue=0.;
	}

	double Price(bool call=1);
	double CumulativeProba(const double& x);
	double MultivariateDensity(const double& x);
	double BivariateNormalDensity(const double& x,const double& y,
		const double& rho,const double& volx,const double& voly,const double& mx,const double& my);
	double BivariateNormalDensity_(double x, void* params);

	double IntBND(const double& borne_inf_x,
								 const double& borne_sup_x,
								 const double& fixed_value_y);

	double IntBND2(const double& borne_inf_y,
								 const double& borne_sup_y);

	double IntBND2_Single(const double& borne_inf_y,
								 const double& borne_sup_y);

};

# endif