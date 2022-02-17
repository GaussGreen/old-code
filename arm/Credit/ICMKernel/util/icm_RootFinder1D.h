
#ifndef _ICM_RootFinder1D_H_
#define _ICM_RootFinder1D_H_

#include "icm_functors.h"
#include <float.h>
#include <math.h>

#include <vector>

# define _BAD_RET_VALUE_	1e99

/*********************************************************************************/
/*! \class  RootFinder1D_t icm_RootFinder1D.h "icm_RootFinder1D.h"
 *  \author 
 *	\version 1.0
 *	\date   May 2004
 *	\brief  Numerical Recipes p430 - Chapter 10. Minimization or Maximization of Functions
 *	\brief  Numerical Recipes
 *	\brief  9.4 Newton-Raphson Method Using Derivative 362 
 *	\brief  NewtonRaphson 
 *	\brief  Using the Newton-Raphson method, find the root of a function known to lie in the interval
 *	\brief  [x1, x2]. The root rtnewt will be refined until its accuracy is known within ±xacc. funcd
 *	\brief  is a user-supplied routine that returns both the function value and the first derivative of the
 *	\brief  function at the point x.
/***********************************************************************************/

template <class F> class RootFinder1D_t
{
public:
	RootFinder1D_t(F f, F df):m_f(f), m_df(df){}
	RootFinder1D_t(F f):m_f(f), m_df(f){}
	virtual ~RootFinder1D_t() {}

	// Numerical Recipes
	// 9.4 Newton-Raphson Method Using Derivative 362 
	// NewtonRaphson 
	// Using the Newton-Raphson method, find the root of a function known to lie in the interval
	// [x1, x2]. The root rtnewt will be refined until its accuracy is known within ±xacc. funcd
	// is a user-supplied routine that returns both the function value and the first derivative of the
	// function at the point x.
	double NewtonRaphson(double x1,
						 double x2,
						 unsigned int NbIterMax=100,
						 double accuracy=0.00001);

	double SimpleNewtonRaphson	(double			x1,
								double			x2,
								unsigned int	NbIterMax=100,
								double			accuracy=0.00001,
								double			slope=0.01,
								double			guess = CREDIT_DEFAULT_VALUE,
								bool			xception = true);

	double NewtonRaphsonModify(double			x0,
							   double			_inf,
							   double			_sup,
							   unsigned int		IterMax = 100,
							   double			acc_dx = 1.E-3,
							   double			acc_f = 1. E-3);

	// 9.4 Newton-Raphson Method Using Derivative 362 
	// NewtonRaphsonWithBisection
	// Using a combination of Newton-Raphson and bisection, find the root of a function bracketed
	// between x1 and x2. The root, returned as the function value rtsafe, will be refined until
	// its accuracy is known within ±xacc. funcd is a user-supplied routine that returns both the
	// function value and the first derivative of the function.
	double NewtonRaphsonWithBisection(double x1,
									 double x2,
									 unsigned int NbIterMax=100,
									 double delta_shape=0.01,
									 double accuracy=0.00001,
									 double f_accuracy=0.01);


	// NewtonRaphson with several guess & choice of working in relative or not. (thx A. Chaix)
	double	NewtonRaphson ( double						target,
							const std::vector<double>&	init,
							double						bound_inf,				
							double						bound_sup,				
							unsigned int				nbitermax,			
							double						tol_fct,
							double						tol_x,
							bool						is_rel_fct = true,
							bool						is_rel_x   = true);


	//Given a function fx defined on the interval from x1-x2 subdivide the interval into n equally
	//spaced segments, and search f	// Dichotomic
	// Using bisection, find the root of a function func known to lie between x1 and x2. The root,
	// returned as double, will be refined until its accuracy is ±accuracy.
	double Dichotomy(double x1,
					 double x2,					 
					 unsigned int NbIterMax=100,
					 double accuracy= 1.E-5,
					 double accuracy_v = 1.E-1,
					 bool	xception = true);

	double Dichotomy(double x1,
					 double x2,
					 int& nb_iter,
					 unsigned int NbIterMax=100,
					 double accuracy= 1.E-5,
					 double accuracy_v = 1.E-1);

	// Given a function func and an initial guessed range x1 to x2, the routine expands the range
	// geometrically until a root is bracketed by the returned values x1 and x2 (in which case zbrac
	// returns 1) or until the range becomes unacceptably large (in which case zbrac returns 0).
	bool ZeroBracket(double& x1,
					double& x2,
					unsigned int NbIterMax=100,
					double factor=1.6);
	
	bool ZeroBracketDecreasing(double& x1,
					double& x2,
					unsigned int NbIterMax=100,
					double factor=1.6, 
					const double& inf = 0.0,
					const double& sup = 10000.0);

	//Given a function fx defined on the interval from x1-x2 subdivide the interval into n equally
	//spaced segments, and search for zero crossings of the function. nb is input as the maximum number
	//of roots sought, and is reset to the number of bracketing pairs xb1[0..nb-1], xb2[0..nb-1]
	//that are found.
	void ZerosBracket(double x1,
					double x2,
					int NbDivisions,
					double* xb1,
					double* xb2,
					int& nb);

	// Given a function func and an initial guessed range x1 to x2, the routine expands the range
	// geometrically until a root is bracketed by the returned values x1 and x2 (in which case zbrac
	// returns 1) or until the range becomes unacceptably large (in which case zbrac returns 0).
	bool ZeroBracketMax(double& x1,
						double& x2,
						double inf,
						double sup,
						unsigned int NbIterMax=100,
						double factor=1.6);


	static bool test(std::string& errStr);

private :
	F		m_f;
	F		m_df;
};
//----------------------------------------------------------------------------
//	Helper. 
template <class F>
inline RootFinder1D_t<F>
RootFinder1D(F f, F df)
{ 
	return RootFinder1D_t<F>(f, df) ; 
}
//----------------------------------------------------------------------------
//	Helper. 
template <class F>
inline RootFinder1D_t<F>
RootFinder1D(F f)
{ 
	return RootFinder1D_t<F>(f, f) ; 
}
//----------------------------------------------------------------------------
template <class F>
bool
RootFinder1D_t<F>::ZeroBracket(double& x1, double& x2, unsigned int NbIterMax, double factor)
{
//	EP_ASSERT_MSG(x1<x2, "Bad initial range in zbrac")
//	EP_ASSERT_MSG(factor>1., "Expansion factor should be greater than one.")
	
	double f1 = m_f(x1);
	double f2 = m_f(x2);
	for (int j=0;j<NbIterMax;j++) {
		if (f1*f2 < 0.) return true;
		else if (fabs(f1) < fabs(f2)) f1 = m_f(x1 += factor*(x1-x2)); 
		else f2 = m_f(x2 += factor * (x2-x1));		
	}

	return false;
}

//----------------------------------------------------------------------------
//		ZeroBracket with check of decreasing value , one >0 and the other <0
//----------------------------------------------------------------------------
template <class F>
bool
RootFinder1D_t<F>::ZeroBracketDecreasing(double& x1, double& x2, unsigned int NbIterMax, double factor, const double& inf, const double& sup)
{
//	EP_ASSERT_MSG(x1<x2, "Bad initial range in zbrac")
//	EP_ASSERT_MSG(factor>1., "Expansion factor should be greater than one.")
	x1 = MAX(x1,0.0);
	double f1 = m_f(x1);
	double f2 = m_f(x2);
	for (int j=0;j<NbIterMax;j++) {
		if (((f1> 0.) && (f2<0.)) || (f1*f2 == 0.) ) return true;
		if ( f2 > 0) {
			x2 += factor * (x2-x1);
			x2 = MIN(x2,sup);
			f2 = m_f(x2);
		} if (f1 < 0.){
			if (f1 < f2 ){ // avance 
				x1 += factor*(x2-x1);
				x1 = MAX(x1,inf);
				f1 = m_f(x1);
			} else { // recule
				x1 += factor*(x1-x2);
				x1 = MAX(x1,inf);
				f1 = m_f(x1);
			}
		}
	}

	return false;
}
//----------------------------------------------------------------------------
template <class F>
bool
RootFinder1D_t<F>::ZeroBracketMax(double& x1, double& x2,double inf,double sup, unsigned int NbIterMax, double factor)
{
//	EP_ASSERT_MSG(x1<x2, "Bad initial range in zbrac")
//	EP_ASSERT_MSG(factor>1., "Expansion factor should be greater than one.")
	
	double x1_ = MAX(x1,inf);
	double x2_ = MIN(x2,sup);

	double f1 = m_f(x1_);
	double f2 = m_f(x2_);
	for (int j=0;j<NbIterMax;j++) 
	{
		if (f1*f2 < 0.) 
			return true;
		else if ((fabs(f1) < fabs(f2)) && (!CHECK_EQUAL(x1_,inf))) 
		{f1 = m_f(x1_ = MAX(x1_ + factor*(x1_-x2_),inf));}
		else
		{f2 = m_f(x2_ = MIN(x2_ + factor * (x2_-x1_),sup));}

		x1 = x1_;
		x2 = x2_;

		if (CHECK_EQUAL(x1_,inf) && CHECK_EQUAL(x2_,sup) && (f1*f2 >= 0.))
			return false;
	}

	return false;
}
//----------------------------------------------------------------------------
template <class F>
void
RootFinder1D_t<F>::ZerosBracket(double x1,
								double x2,
								int NbDivisions,
								double* xb1,
								double* xb2,
								int& nb)
{
	double x(0.);
	int nbb=0;

	//Determine the spacing appropriate to the mesh.
	double dx = (x2-x1)/((double) NbDivisions); 
	double fp = m_f(x=x1);
	
	//Loop over all intervals
	for (int i=0;i<NbDivisions;i++) { 
		double fc = m_f(x += dx);
		//If a sign change occurs then record values for the bounds.
		if (fc*fp <= 0.) { 
			xb1[nbb] = x-dx;
			xb2[nbb] = x;
			nbb++;
			if(nb == nbb) break;
		}
		fp = fc;
	}
	
	// reset number of roots found
	nb = nbb;
}
//----------------------------------------------------------------------------
template <class F>
double
RootFinder1D_t<F>::Dichotomy(double x1,
							 double x2,
							 unsigned int NbIterMax,
							 double accuracy,
							 double accuracy_v,
							 bool	xception)

{
//	EP_ASSERT(x1<x2)
//	EP_ASSERT(accuracy>0.)

	double f = m_f(x1);
	double fmid = m_f(x2);
	double xmid(x1), dx(0.);
	
	if (f*fmid>0.0)
	{
		if (!xception)
		{
			#ifdef _DEBUG
			ICMLOG("Dichotomy failed for (f*fmid>0.0) "<<" x1="<<x1<<" x2="<<x2<<" f="<<f<<" fmid="<<fmid);
			#endif
			return CREDIT_DEFAULT_VALUE;
		}
		else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Root must be bracketed for dichotomy.");
	}

	//Orient the search so that f>0 lies at x+dx. 
	double ret = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2); 
	for (int nb_iter=0;nb_iter<NbIterMax;nb_iter++) {
		//Bisection loop.
		fmid = m_f(xmid= ret+(dx *= 0.5)); 
		if (fmid <= 0.0) ret = xmid;
		else if (((fabs(dx) < accuracy) && (fabs(fmid)<accuracy_v)) || fmid == 0.0) 
			return ret;
	}
//	EP_THROWEX("Too many iteration in dichotomy");
	return CREDIT_DEFAULT_VALUE; //Never get here.
}
//----------------------------------------------------------------------------//----------------------------------------------------------------------------
template <class F>
double
RootFinder1D_t<F>::Dichotomy(double x1,
							 double x2,
							 int& nb_iter,
							 unsigned int NbIterMax,
							 double accuracy,
							 double accuracy_v)

{
//	EP_ASSERT(x1<x2)
//	EP_ASSERT(accuracy>0.)

	double f = m_f(x1);
	double fmid = m_f(x2);
	double xmid(x1), dx(0.);
	if(f == 0.)
		return x1;
	if(fmid == 0.)
		return x2;

	if (f*fmid>0.0)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Root must be bracketed for dichotomy.");

	//Orient the search so that f>0 lies at x+dx. 
	double ret = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2); 
	for (nb_iter=0;nb_iter<NbIterMax;nb_iter++) {
		//Bisection loop.
		fmid = m_f(xmid=ret+(dx *= 0.5)); 
		if (fmid <= 0.0) ret = xmid;
		else if (((fabs(dx) < accuracy) && (fabs(fmid)<accuracy_v)) || fmid == 0.0) 
			return ret;
	}
//	EP_THROWEX("Too many iteration in dichotomy");
	return _BAD_RET_VALUE_; //Never get here.
}

//----------------------------------------------------------------------------
template <class F>
double
RootFinder1D_t<F>::NewtonRaphson(double			x1,
								 double			x2,
								 unsigned int	NbIterMax,
								 double			accuracy)
{
	if (x1>=x2)
      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"x1>=x2");

	if (accuracy<0.)
      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"accuracy<0.");

	//Initial guess
	double ret = 0.5*(x1+x2);

	for (int j=0;j<NbIterMax;j++) 
	{
		double f	= m_f(ret);
		double df	= m_df(ret);

	if (fabs(df)<=DBL_EPSILON)
      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"fabs(df)<=DBL_EPSILON");

		double dx = f/df;
		ret -= dx;

	if ((x1-ret)*(ret-x2) < 0.)
      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Jumped out of brackets in NewtonRaphson");

	//Convergence
	if ((fabs(dx) < accuracy) && (fabs(f) < 1.)) 
		return ret; 
	}	

	return (ret);
}

template <class F>
double
RootFinder1D_t<F>::SimpleNewtonRaphson	(double			x1,
										double			x2,
										unsigned int	NbIterMax,
										double			accuracy,
										double			slope,
										double			guess,
										bool			xception)
{
	if (x1>=x2)
	{if (!xception) 
	{
		#ifdef _DEBUG
		ICMLOG("SimpleNewtonRaphson failed for x1>=x2 "<<" x1="<<x1<<" x2="<<x2);
		#endif
		return CREDIT_DEFAULT_VALUE;
	}
	  else 
		  throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"x1>=x2");}

	if (accuracy<0.)
	{if (!xception)
	{return CREDIT_DEFAULT_VALUE;}
	  else
		  throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"accuracy<0.");}

	//Initial guess
	double ret = 0.5*(x1+x2);
	if (guess!=CREDIT_DEFAULT_VALUE) ret=guess;

	for (int j=0;j<NbIterMax;j++) 
	{
		double f	= m_f(ret);
		double df	= (m_f(ret+slope) - f)/slope;

	if (fabs(df)<=DBL_EPSILON)
	{ if (!xception)
	{
		#ifdef _DEBUG
		ICMLOG("SimpleNewtonRaphson failed for fabs(df)<=DBL_EPSILON "<<" df="<<df<<" mid="<<ret);
		#endif
		return CREDIT_DEFAULT_VALUE;
	}
	 else
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"fabs(df)<=DBL_EPSILON");}

		double dx = f/df;
		ret -= dx;

	if ((x1-ret)*(ret-x2) < 0.)
	{ if (!xception) 
	{
		#ifdef _DEBUG
		ICMLOG("SimpleNewtonRaphson failed for (x1-ret)*(ret-x2)<0. "<<" x1="<<x1<<" x2="<<x2<<" mid="<<ret);
		#endif
		return CREDIT_DEFAULT_VALUE;
	}
	else	
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Jumped out of brackets in NewtonRaphson");}

	//Convergence
	if ((fabs(dx) < accuracy) && (fabs(f) < 1.)) 
		return ret; 
	}	

	return (ret);
}


template <class F>
double
RootFinder1D_t<F>::NewtonRaphsonModify(double			x0,
									   double			_inf,
									   double			_sup,
									   unsigned int		IterMax,
									   double			acc_dx,
									   double			acc_f)
{
	if (_inf>=_sup)
      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"_inf>=_sup");

	if (acc_dx<0.)
      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"acc_dx<0.");

	if (acc_f<0.)
      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"acc_f<0.");

	//Initial guess
	double ret = x0;

	for (int j=0;j<IterMax;j++) 
	{
		double f	= m_f(ret);
		double df	= m_df(ret);

	if (fabs(df)<=DBL_EPSILON)
      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"fabs(df)<=DBL_EPSILON");

		double dx = f/df;
		ret -= dx;

	if (ret < _inf)
		ret = _inf;
	else if (ret > _sup)
		ret = _sup;

	//Convergence
	if ((fabs(dx) < acc_dx) && (fabs(f) < acc_f)) 
		return ret; 
	}	

	return (ret);
}

//----------------------------------------------------------------------------
template <class F>
double
RootFinder1D_t<F>::NewtonRaphson (  double						target,
									const std::vector<double>&	init,
									double						bound_inf,				
									double						bound_sup,				
									unsigned int				nbitermax,			
									double						tol_fct,
									double						tol_x,
									bool						is_rel_fct,
									bool						is_rel_x)								
{	
	int				p	= 0;
	double			phi	= init[0];
	unsigned int	cpt	= 1;
		
	while (true) {	
		double f	= m_f (phi) - target ;
		double df	= m_df(phi);
		if (is_rel_fct) {
			if (fabs(f/target) <tol_fct) break;
		}
		else if (fabs(f) <tol_fct) break;
			
		phi -= f/df;
		
		if (is_rel_x) {
			if ( fabs(f/df/phi) < tol_x  &&  fabs(df)>DBL_EPSILON) break;  // sortie de la boucle while
		}
		else if ( fabs(f/df) < tol_x &&  fabs(df)>DBL_EPSILON) break;
				
		if( phi<=bound_inf							|| 
			phi>=bound_sup							|| 
			(is_rel_x && fabs(phi)<DBL_EPSILON) || 
			fabs(df)<DBL_EPSILON || p==nbitermax) {	
			if(cpt<init.size())	{
				phi	= init[cpt];
				cpt ++ ;
				p	= -1;
			}
			else EP_THROWEX ("NewtonRaphson::SearchRoot : all initial values failed."); 
		}
		p++;
	}
	
	return phi;
}
//----------------------------------------------------------------------------
template <class F>
double
RootFinder1D_t<F>::NewtonRaphsonWithBisection(double x1,
											  double x2,
											  unsigned int NbIterMax,
											  double delta_shape,
											  double accuracy,
											  double f_accuracy)
{	
	double xl(0.), xh(0.);

	double fl = m_f(x1);
	double fh = m_f(x2);
	if(fl == 0.)
		return x1;
	if(fh ==0.)
		return x2;
	
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Root must be bracketed in rtsafe");
		//EP_THROWEX("Root must be bracketed in rtsafe")
	else if	(fl == 0.0) return x1;
	else if (fh == 0.0) return x2;
	else if (fl < 0.0) { //Orient the search so that f(xl) < 0.
		xl=x1;
		xh=x2;
	} 
	else {
		xh=x1;
		xl=x2;
	}

	//Initialize the guess for root
	double ret		= 0.5*(x1+x2);
	double dxold	= fabs(x2-x1); //the “stepsize before last,”
	double dx		= dxold; // and the last step.

	double f	= m_f(ret);
	double df	= (m_f(ret+delta_shape) - f)/delta_shape;
//	EP_ASSERT(fabs(df)>DBL_EPSILON)

	//Loop over allowed iterations.
	for (int j=0;j<NbIterMax;j++) { 

		if ( ( ((ret-xh)*df-f) * ((ret-xl)*df-f ) > 0.0) || //Bisect if Newton out of range,
			 ( fabs(2.*f) > fabs(dxold*df)) ) { //or not decreasing fast enough.
			dxold	= dx;
			dx		= 0.5*(xh-xl);
			ret		= xl+dx;
		}
		else { //Newton step acceptable. Take it.
			dxold	= dx;
			dx		= f/df;
			ret		-= dx;
		}
		//Convergence criterion.
		if ( (fabs(dx) < accuracy) && fabs(f) < f_accuracy) 
			return ret;
		
		f = m_f(ret);
		df = (m_f(ret+delta_shape) - f)/delta_shape;
		//EP_ASSERT(fabs(df)>DBL_EPSILON)
		
		//The one new function evaluation per iteration.
		if (f < 0.0) xl = ret; //Maintain the bracket on the root.
		else xh = ret;
	}
	
//	EP_THROWEX("Maximum number of iterations exceeded in NewtonRaphsonSafe");

	return ret; //Never get here.
}

//----------------------------------------------------------------------------
// For test purpose
class ATest {
public: 
	ATest(){};
	virtual ~ATest(){};
	double f(double x) {return (x-3.)*(x-10.);}
	double df(double x) {return 2.*x-13.;}
};
//----------------------------------------------------------------------------
template <class F>
bool
RootFinder1D_t<F>::test(std::string& errStr)
{
	bool ret(true);
	
	errStr = "Testing RootFinder1D";
	
	ATest a;

	double val = RootFinder1D(f1::mem_call(&ATest::f, a), f1::mem_call(&ATest::df, a)).NewtonRaphson(0., 9., 100, 0.000001);	
//	SDMTEST((val>2.9999)&&(val<3.0001))
	
	val = RootFinder1D(f1::mem_call(&ATest::f, a), f1::mem_call(&ATest::df, a)).NewtonRaphsonWithBisection(0., 9., 100, 0.000001);	
//	SDMTEST((val>2.9999)&&(val<3.0001))
	
	val = RootFinder1D(f1::mem_call(&ATest::f, a)).Dichotomy(0., 9., 100, 0.000001);	
//	SDMTEST((val>2.9999)&&(val<3.0001))
	
	std::vector<double> init;
	init.push_back(-30.); init.push_back(1.);
	val = RootFinder1D(f1::mem_call(&ATest::f, a), f1::mem_call(&ATest::df, a)).NewtonRaphson ( 0., init, -4., 5., 20, 0.000001, 0.000001, false);	
//	SDMTEST((val>2.9999)&&(val<3.0001))
	
	double x1(1.), x2(2.);
	bool res(false);
	res = RootFinder1D(f1::mem_call(&ATest::f, a)).ZeroBracket(x1, x2, 100, 2);	
//	SDMTEST(res)
//	SDMTEST((x1<3.)&&(x2>3.))
	x1 = -10.;
	x2 = 50;
	int nb = 3;
	double* xb1 = new double [3];
	double* xb2 = new double [3]; 
	RootFinder1D(f1::mem_call(&ATest::f, a)).ZerosBracket(x1, x2, 1000, xb1, xb2, nb);	
//	SDMTEST(nb==2)
//	SDMTEST((xb1[0]<3.)&&(xb2[0]>3.)&&(xb2[0]<10.))
//	SDMTEST((xb1[1]>3.)&&(xb1[1]<10.)&&(xb2[1]>10.))
	delete [] xb1;
	delete [] xb2;

	return ret;
}
#endif	// _F1_RootFinder1D_H_


