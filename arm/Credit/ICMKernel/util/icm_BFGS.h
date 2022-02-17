
#ifndef _ICM_BFGS_H_
#define _ICM_BFGS_H_

#include "icm_functors.h"

/*********************************************************************************/
/*! \class  BFGS_t icm_BFGS.h "icm_BFGS.h"
 *  \author 
 *	\version 1.0
 *	\date   June 2004
 *	\brief  Numerical Recipes p430 - Chapter 10. Minimization or Maximization of Functions
 *	\brief	Given a starting point p[1..n] that is a vector of length n, the Broyden-Fletcher-Goldfarb-
 *	\brief	Shanno variant of Davidon-Fletcher-Powell minimization is performed on a function func, using
 *	\brief  its gradient as calculated by a routine dfunc.
 *	\brief	The convergence requirement on zeroing the gradient is input as gtol.
 *	\brief	Returned quantities are p[1..n] (by ptr, the location of the minimum), and return the minimum value of the function.
 *	\brief	The routine lnsrch is called to perform approximate line minimizations.
/***********************************************************************************/

template <class F, class DF>
class BFGS_t
{
public:
	BFGS_t(F f, DF df): m_dim(0),
						m_f	(f),
						m_df(df),
						m_dg(NULL),
						m_g	(NULL),
						m_hdg(NULL),
						m_pnew(NULL),
						m_xi(NULL),
						m_hessin(NULL) {}

	virtual ~BFGS_t() {clear();}

	bool	Minimize(double* p,
					int		n,
					double	gtol,
					int*	NbIter,
					double*	MinimizeVal,
					int		NbIterMax	= 200,
					double	epsilon		= 1e-8,
					double	precision	= 1e-7,
					double	stepmax		= 0.1);


	bool	line_search(double*	xold,
					   double&	f,
					   double*	g,
					   double*	p,
					   double*	x,
//					   double&	f,
					   double	stepmax, 
					   double	epsilon,
					   double	precision);

protected :
	void	init(int n);
	void	clear();

public:
	static bool test(std::string& errStr);

private :
	F	m_f;	// double m_f(double*)
	DF	m_df;	// void m_df(double* x, double* df)
	int			m_dim;
	double*		m_dg;
	double*		m_g;
	double*		m_hdg;
	double*		m_pnew;
	double*		m_xi;
	double**	m_hessin;
};
//--------------------------------------------------------------------------------------------------------------------------------------
//	Helper Function.
template <class F, class DF> inline
BFGS_t<F, DF>
BFGS(F f, DF df) { return BFGS_t<F, DF>(f, df); }
//--------------------------------------------------------------------------------------------------------------------------------------
template <class F, class DF> 
void
BFGS_t<F, DF>::init(int n)
{
	clear();
	if (n<=0) ICMTHROW(ERR_INVALID_ARGUMENT,"BFGS: n should be >0 "<<n); 
	m_dim = n;

	m_dg	 = new double	[m_dim];
	m_g		 = new double	[m_dim];
	m_hdg	 = new double	[m_dim];
	m_pnew	 = new double	[m_dim];
	m_xi	 = new double	[m_dim];
	m_hessin = new double*	[m_dim];
	for(int j=0;j<n;j++) m_hessin[j] = new double [m_dim];
}
//--------------------------------------------------------------------------------------------------------------------------------------
template <class F, class DF> 
void
BFGS_t<F, DF>::clear()
{
	if(m_dg!=NULL)	 {delete [] m_dg;	m_dg=NULL;}
	if(m_g!=NULL)	 {delete [] m_g;	m_g=NULL;}
	if(m_pnew!=NULL) {delete [] m_pnew; m_pnew=NULL;}
	if(m_hdg!=NULL)	 {delete [] m_hdg;	m_dg=NULL;}
	if(m_xi!=NULL)	 {delete [] m_xi;	m_xi=NULL;}
	for(int k=0;k<m_dim;k++) {
		 if(m_hessin[k]!=NULL) {delete [] m_hessin[k]; m_hessin[k]=NULL;} 
	}
	if(m_hessin!=NULL)	{delete [] m_hessin; m_hessin=NULL;} 
}
//--------------------------------------------------------------------------------------------------------------------------------------
template <class F, class DF> 
bool
BFGS_t<F, DF>::Minimize(double*	p,
						int		n,
						double	gtol,
						int*	NbIter,
						double*	MinimizeVal,
						int		NbIterMax,
						double	epsilon,
						double	precision,
						double	stp
						)
{
	init(n);

	// Calculate starting function value and gradient,
	double fp = m_f(p); 
	m_df(p,m_g);

	// initialize the inverse Hessian to the unit matrix.
	double sum(0.);
	for (int i=0;i<m_dim;i++) {
		for (int j=0;j<m_dim;j++) m_hessin[i][j] = 0.;
		m_hessin[i][i] = 1.;
		//Initial line direction.
		m_xi[i] = -m_g[i]; 
		sum += p[i]*p[i];
	}

	double stpmax = stp * std::_cpp_max(sqrt(sum), (double)m_dim);

	// Main loop over the iterations.
	for (int its=1;its<=NbIterMax;its++) { 

		line_search(p, fp, m_g, m_xi, m_pnew, stpmax, epsilon, precision);
		// The new function evaluation occurs in lnsrch;
		// save the function value in fp for the next line search.
		// It is usually safe to ignore the value of check.
		for (i=0;i<m_dim;i++) {
			m_xi[i]	= m_pnew[i]-p[i]; //Update the line direction,
			p[i]	= m_pnew[i];		//and the current point.
		}
		//Test for convergence on .x.
		double test = 0.; 
		for (i=0;i<m_dim;i++) {
			double temp = fabs(m_xi[i])/std::_cpp_max(fabs(p[i]),1.);
			if (temp > test) test = temp;
		}
		if (test < precision) {
			ICMLOG("BFGS::Minimize : Convergence achieved in " << its << " iterations.")
			
			*NbIter	=	its;
			*MinimizeVal	=	fp;
			return true;
		}

		for (i=0;i<m_dim;i++) m_dg[i] = m_g[i]; //Save the old gradient,
		m_df(p,m_g); //and get the new gradient.

		//Test for convergence on zero gradient.
		test = 0.; 
		double den = std::_cpp_max(fp, 1.);
		for (i=0;i<m_dim;i++) {
			double temp = fabs(m_g[i])*std::_cpp_max(fabs(p[i]),1.0)/den;
			if (temp > test) test=temp;
		}
		if (test < gtol) {
			ICMLOG("BFGS::Minimize : Convergence achieved in " << its << " iterations.") 

			*NbIter	=	its;
			*MinimizeVal	=	fp;
			return true;
		}
		for (i=0;i<m_dim;i++) m_dg[i] = m_g[i]-m_dg[i]; //Compute difference of gradients,
		for (i=0;i<m_dim;i++) { //and difference times current matrix.
			m_hdg[i] = 0.;
			for (int j=0;j<m_dim;j++) m_hdg[i] += m_hessin[i][j]*m_dg[j];
		}
		// Calculate dot products for the denominators.
		double fac(0.), fae(0.), sumdg(0.), sumxi(0.);
		for (i=0;i<m_dim;i++) {
			fac		+= m_dg[i]*m_xi[i];
			fae		+= m_dg[i]*m_hdg[i];
			sumdg	+= m_dg[i]*m_dg[i];
			sumxi	+= m_xi[i]*m_xi[i];
		}
		//Skip update if fac not sufficiently positive.
		if (fac > sqrt(epsilon*sumdg*sumxi)) { 
			fac = 1./fac;
			double fad = 1./fae;
			//The vector that makes BFGS different from DFP:
			for (i=0;i<m_dim;i++) m_dg[i]=fac*m_xi[i]-fad*m_hdg[i];
			//The BFGS updating formula:
			for (i=0;i<m_dim;i++) { 
				for (int j=i;j<m_dim;j++) {
					m_hessin[i][j] += fac*m_xi[i]*m_xi[j] -fad*m_hdg[i]*m_hdg[j]+fae*m_dg[i]*m_dg[j];
					m_hessin[j][i] = m_hessin[i][j];			
				}
			}
		}

		//Now calculate the next direction to go,
		for (i=0;i<m_dim;i++) { 
			m_xi[i] = 0.;
			for (int j=0;j<m_dim;j++) m_xi[i] -= m_hessin[i][j]*m_g[j];	
		}
	}	//and go back for another iteration.

//	ICMTHROW(ERR_INVALID_ARGUMENT,"too many iterations in BFGS::Minimize");

	*NbIter	=	its;
	*MinimizeVal	=	fp;
	return false;

}
//--------------------------------------------------------------------------------------------------------------------------------------
/*Given an n-dimensional point xold[1..n], the value of the function and gradient there, fold
and g[1..n], and a direction p[1..n], finds a new point x[1..n] along the direction p from
xold where the function func has decreased “sufficiently.”
The new function value is returned in f.
Stepmax is an input quantity that limits the length of the steps so that you do not try to
evaluate the function in regions where it is undefined or subject to overflow.
p is usually the Newton direction.
The output quantity check is false (0) on a normal exit. It is true (1) when
x is too close to xold.
In a minimization algorithm, this usually signals convergence and can
be ignored. However, in a zero-finding algorithm the calling program should check whether the
convergence is spurious.*/
template <class F, class DF> 
bool
BFGS_t<F, DF>::line_search(double*	xold,
					  double&	fold,
					  double*	g,
					  double*	p,
					  double*	x, // ret
					  //double&	f,
					  double	stpmax,
					  double	epsilon,
					  double	precision)
{
	double alam2(0.), f2(0.), tmplam(0.);

	double sum(0.);
	for (int i=0;i<m_dim;i++) sum += p[i]*p[i];
	sum = sqrt(sum);

	//Scale if attempted step is too big.
	if (sum > stpmax) for (i=0;i<m_dim;i++) p[i] *= stpmax/sum; 
	
	double slope(0.);
	for(i=0;i<m_dim;i++) slope += g[i]*p[i];
//	if (slope <0) ICMTHROW(ERR_INVALID_ARGUMENT,"Roundoff problem in lnsrch."); 
	// Compute fmin.
	double test = 0.; 
	for (i=0;i<m_dim;i++) {
		double temp = fabs(p[i])/std::_cpp_max(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}

	double alamin = epsilon/test;
	//Always try full Newton step .rst.
	double alam = 1.; 
	for (;;) { //Start of iteration loop.
		for (i=0;i<m_dim;i++) x[i] = xold[i] + alam*p[i];
		double f = m_f(x);
		if (alam < alamin) { //Convergence on Dx. For zero finding, the calling program shouldverify the convergence.
			for (i=0;i<m_dim;i++) x[i] = xold[i];
			fold = f;
			return true;
		}
		else if (f <= fold+precision*alam*slope) {
			fold = f;
			return false; //Sufficient function decrease.
		}
		else { //Backtrack.
			if (alam == 1.0) tmplam = -slope/(2.*(f-fold-slope)); //First time.
			else { //Subsequent backtracks.
				double rhs1 = f-fold-alam*slope;
				double rhs2 = f2-fold-alam2*slope;
				double a	= (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				double b	= (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.) tmplam = -slope/(2.*b);
				else {
					double disc = b*b-3.*a*slope;
					if (disc < 0.)		tmplam = 0.5*alam;
					else if (b <= 0.)	tmplam = (-b+sqrt(disc))/(3.*a);
					else				tmplam = -slope/(b+sqrt(disc));
				}
				if (tmplam > 0.5*alam) tmplam = 0.5*alam; // Lambda<= 0.5 Lambda1.
			}
		}
		alam2	= alam;
		f2		= f;
		alam	= std::_cpp_max(tmplam, 0.1*alam); // Lambda >= 0.1 Lambda1.
	} //Try again.

}
//--------------------------------------------------------------------------------------------------------------------------------------
// For test purpose
class Banana {
public: 
	Banana(){};
	virtual ~Banana(){};
	double	f(double* x) {return 100.*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0])+(1.-x[0])*(1.-x[0]);}
	void	df(double* x, double* ret) {
		ret[0] = -400.* x[0] * (x[1]-x[0]*x[0]) - 2.*(1.-x[0]);
		ret[1] = 200.* (x[1]-x[0]*x[0]);
	}
};
//----------------------------------------------------------------------------
template <class F, class DF>
bool
BFGS_t<F, DF>::test(std::string& errStr)
{
	bool ret(true);
	
	errStr = "Testing BFGS";
	Banana a;
	double* P = new double [2];
	P[0] = -1.2; P[1] = 5; 

//	ff1::mem_call_t<double, ABFGS , double*> f = ff1::mem_call(&Banana::f, a);
//	ff1::mem_call_t<double*, ABFGS , double*> df = ff1::mem_call(&Banana::df, a);

	double val = BFGS(ff1::mem_call(&Banana::f, a), f1::mem_call(&Banana::df, a)).Minimize(P, 2, 0.00001);	
//	SDMTEST((P[0]>0.9999)&&(P[0]<1.0001))
//	SDMTEST((P[1]>0.9999)&&(P[1]<1.0001))
//	SDMTEST((val>-0.000001)&&(val<0.000001))

	delete [] P;

	return ret;
}
#endif	// _BFGS_H_

