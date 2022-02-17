
/********************************************************************************/
/*! \file icm_RootFinder1D.h
 *  \brief  Newton-Raphson Method for Nonlinear Systems of Equations
 *  \author 
 *	\version 1.0
 *	\date   May 2004
 /*********************************************************************************/


#ifndef _ICM_RootFinderND_H_
#define _ICM_RootFinderND_H_

#include "icm_functors.h"
#include <float.h>
#include <math.h>
#include "ICMKernel\glob\icm_enums.h" 
#include "ICMKernel\util\icm_macro.h" 
#include "ICMKernel\optim\uncoptim.h"
#include "ICMKernel\optim\lm.h"
#include "ICMKernel\optim\pswarm.h"
#include "ICMKernel\random\icm_randomnag.h"
#include "ARMKernel\util\rand-gen.h"
#include "gpnumlib\levmarq.h"

#ifdef MPI
#include "mpi.h"
#endif
#ifdef MPE
#include "mpe.h"
#endif
#ifdef AMPL
#include "nlp.h"
#include "asl.h"
#include "getstub.h"
#else
#include <stdio.h>
#endif

#include <signal.h>



#ifdef AMPL
char *pswarm_opt_s(Option_Info *, keyword *, char *);
char *pswarm_opt_i(Option_Info *, keyword *, char *);
char *pswarm_opt_d(Option_Info *, keyword *, char *);
#endif



#include <vector>
#include <nage04.h>	

#ifndef MIN
# define MIN(a,b) (((a)<(b))?(a):(b)) 
#endif
#ifndef MAX
# define MAX(a,b) (((a)>(b))?(a):(b)) 
#endif

/*********************************************************************************/
/*! \class  RootFinderND_t icm_RootFinderND.h "icm_RootFinderND.h"
 *  \author 
 *	\version 1.0
 *	\date   May 2004
 *	\file   icm_RootFinderND.h
 *	\brief  Numerical Recipes
 *	\brief  9.6 Newton-Raphson Method for Nonlinear Systems of Equations
 *	\brief  Given an initial guess x[1..n] for a root in n dimensions, take ntrial Newton-Raphson steps
 *	\brief  to improve the root. Stop if the root converges in either summed absolute variable increments
 *	\brief  tolx or summed absolute function values tolf. */
/***********************************************************************************/

typedef void (* OptFunction)(int*,double*,double*);
typedef void (*LMFunction)(double *p, double *hx, int m, int n, void *adata);
typedef void (NAG_CALL *NAG_E04UNC_FUN)(long M, long N, double X[],double F[],double fjac[],long tdfjac, Nag_Comm* comm);
typedef double (*PSWGENFunction)(int, double *, double *, double *);

void GradientLM(double *p, double *hx, int m, int n, void *adata);

void NAG_CALL FuncConfun(long n, long ncnlin, long needc[], double x[],
                         double conf[], double cjac[],
                         Nag_Comm* comm);
//static void __stdcall confun(long n, long ncnlin, long needc[], double x[], double conf[], double conjac[], Nag_Comm *comm);

namespace OptimTools
{
	class Opt: public UnconstrainedOptim 
	{
	public:

		NAG_E04JBC_FUN functionNag;
		Nag_Comm*	   Nagcomm;

		Opt(int n,int lsmethod):UnconstrainedOptim(n,0,0,lsmethod) 
		{}

		virtual void F() 
		{functionNag(n,x,&f,g,Nagcomm);};

		virtual void G() 
		{
			F();

			double f0 = f;

			for (int i=0; i<n; i++)
			{	
			x[i] += GRAD_EPSILON;
			F();
			g[i] = (f - f0)/GRAD_EPSILON;
			x[i] -= GRAD_EPSILON;
			}

			f = f0;
		}

	};

	static double *scale_;
	static void* GansoLibContext;
	static LMFunction LMfunction;

	double  NagMultiOptimisator(void* context,
							  NAG_E04JBC_FUN f,
							  std::vector<double>& X,
							  std::vector<double>& bound_inf_X,
							  std::vector<double>& bound_sup_X,
							  double tol = 1.e-5,
							  int maxiter = 200,
							  double step_max = 1.e-5,
							  double linesearch=0.5);

	double  NagMultiOptimisator_BFGS(void* context,
							  NAG_E04UNC_FUN f,
							  std::vector<double>& X,
							  std::vector<double>& bound_inf_X,
							  std::vector<double>& bound_sup_X,
							  int nbFuncs=1,
							  double tol = 1.e-5,
							  int maxiter = 200);

	double  GansoMultiOptimisator(void* context,
							  OptFunction f,
							  std::vector<double>& X,
							  std::vector<double>& bound_inf_X,
							  std::vector<double>& bound_sup_X,
							  double tol = 1.e-5,
							  int maxiter = 200,
							  double step_max = 1.e-5,
							  double linesearch=0.5);

	double  UNCMultiOptimisator(void* context,
							  NAG_E04JBC_FUN f,
							  std::vector<double>& X,
							  std::vector<double>& bound_inf_X,
							  std::vector<double>& bound_sup_X,
							  double tol = 1.e-5,
							  int maxiter = 200,
							  double step_max = 1.e-5,
							  double linesearch=0.5);

	double  LMMultiOptimisator(void* context,
							  LMFunction f,
							  std::vector<double>& X,
							  std::vector<double>& bound_inf_X,
							  std::vector<double>& bound_sup_X,
							  int nbFuncs,
							  double tol = 1.e-5,
							  int maxiter = 200,
							  double step_max = 1.e-5,
							  double linesearch=0.5);

	double  PSwarmMultiOptimisator(void* &context,
							  PSWGENFunction f,
							  std::vector<double>& X,
							  std::vector<double>& bound_inf_X,
							  std::vector<double>& bound_sup_X,
							  int nbFuncs,
							  double tol = 1.e-5,
							  int maxiter = 200,
							  double step_max = 1.e-5,
							  double linesearch=0.5);

	// ---------------------------------------------------
	// Generic Optimizer
	// ---------------------------------------------------
	void GenerateX0(ARM_MMTGenerator& RandGen,
					std::vector<double>& X,
					std::vector<double>& bound_inf_X,
					std::vector<double>& bound_sup_X);

	double  MultiOptimisator(qOPTIMIZE_TYPE& type,
							  void* context,
							  void* f,
							  std::vector<double>& X,
							  std::vector<double>& bound_inf_X,
							  std::vector<double>& bound_sup_X,
							  int nbFuncs=1,
							  double tol = 1.e-5,
							  int maxiter = 200,
							  double step_max = 1.e-5,
							  double linesearch=0.5);

	double AleaOptimSearch(qOPTIMIZE_TYPE& type,
							  void* context,
							  void* f,
							  std::vector<double>& X,
							  std::vector<double>& bound_inf_X,
							  std::vector<double>& bound_sup_X,
							  int nbFuncs=1,
							  long nbsimuls = 100,
							  double maxtime = 10000,
							  double tol = 1.e-5,
							  int maxiter = 200,
							  double step_max = 1.e-5,
							  double linesearch=0.5,
							  double TRACE = 0.);

}


template <class F> class RootFinderND_t
{
public:

	RootFinderND_t(F f):m_dim(0),
						m_f(f),
						m_indx(NULL),
						m_p(NULL),
						m_vv(NULL),
						m_fvec(NULL),
						m_fjac(NULL),
						m_g(NULL),
						m_xold(NULL) {}

	virtual ~RootFinderND_t() {clear();}

	// Numerical Recipes
	// 9.6 Newton-Raphson Method for Nonlinear Systems of Equations
	// Given an initial guess x[1..n] for a root in n dimensions, take ntrial Newton-Raphson steps
	// to improve the root. Stop if the root converges in either summed absolute variable increments
	// tolx or summed absolute function values tolf.
	double	NewtonMultiDim(double* x,
							int n,
							int ntrial=300,
							double tolx=1.e-7,
							double tolf=1.e-7);

	double	GloballyConvergentNewton(double* x,
									 int n,
									 bool& check,
									 int nbtrial=300,
									 double x_precision=1.e-7,
									 double f_precision=1.e-7,
									 double step_max=1.,
									 double tol_min=1.e-7);

	static bool test(std::string& errStr);

protected:
	void	init(int n);
	void	clear();
	void	lubksb(double **a, int n, int *indx, double *b);
	void	ludcmp(double **a, int n, int *indx, int& d);
	double	fmin(double* x);
	void	fdjac(double* x, double* fvec, double** fjac, double eps=1.e-7);
	bool	line_search(double*	xold,
					   double&	f,
					   double*	g,
					   double*	p,
					   double*	x,
					   double	stepmax, 
					   double	x_precision,
					   double	f_precision);

private :
	F		m_f;

	int			m_dim;
	int*		m_indx;
	double*		m_p;
	double*		m_g;
	double*		m_xold;
	double*		m_fvec;
	double*		m_vv;
	double**	m_fjac;
};
//----------------------------------------------------------------------------
//	Helper. 
template <class F>
inline RootFinderND_t<F>
RootFinderND(F f)
{ 
	return RootFinderND_t<F>(f) ; 
}
//--------------------------------------------------------------------------------------------------------------------------------------
template <class F> 
void
RootFinderND_t<F>::init(int n)
{
	clear();
//	EP_ASSERT(n>0)
	m_dim = n;

	m_indx	 = new int		[m_dim];
	m_p		 = new double	[m_dim];
	m_vv	 = new double	[m_dim];
	m_fvec	 = new double	[m_dim];
	m_g		 = new double	[m_dim];
	m_xold	 = new double	[m_dim];
	m_fjac	 = new double*	[m_dim];

	for(int j=0;j<m_dim;j++) m_fjac[j] = new double [m_dim];
}
//--------------------------------------------------------------------------------------------------------------------------------------
template <class F> 
void
RootFinderND_t<F>::clear()
{
	if(m_indx!=NULL)	{delete [] m_indx;	m_indx	=NULL;}
	if(m_p!=NULL)		{delete [] m_p;		m_p		=NULL;}
	if(m_vv!=NULL)		{delete [] m_vv;	m_vv	=NULL;}
	if(m_xold!=NULL)	{delete [] m_xold;	m_xold	=NULL;}
	if(m_g!=NULL)		{delete [] m_g;		m_g		=NULL;}
	if(m_fvec!=NULL)	{delete [] m_fvec;	m_fvec	=NULL;}
	for(int k=0;k<m_dim;k++) {
		 if(m_fjac[k]!=NULL) {delete [] m_fjac[k]; m_fjac[k]=NULL;} 
	}
	if(m_fjac!=NULL)	{delete [] m_fjac; m_fjac=NULL;} 
}
//--------------------------------------------------------------------------------------------------------------------------------------
//Return 0.5 F · F at x. 
template <class F> 
double
RootFinderND_t<F>::fmin(double* x)
{
	m_f(x, m_fvec);

	double sum(0.);
	for (int i=0;i<m_dim;i++) sum += m_fvec[i] * m_fvec[i];

	return 0.5*sum;
}
//--------------------------------------------------------------------------------------------------------------------------------------
//Computes forward-difference approximation to Jacobian. On input, x[1..n] is the point at
//which the Jacobian is to be evaluated, fvec[1..n] is the vector of function values at the
//point, and vecfunc(n,x,f) is a user-supplied routine that returns the vector of functions at
//x. On output, df[1..n][1..n] is the Jacobian array.
template <class F> 
void
RootFinderND_t<F>::fdjac(double* x,
							  double* fvec,
							  double** fjac,
							  double eps)
{
	for (int j=0;j<m_dim;j++) {

		double temp = x[j];
		double h = std::_cpp_max (eps*fabs(temp), eps);
		x[j] += h; //Trick to reduce finite precision error.
		m_f(x, m_vv);
		x[j] = temp;
		// Forward difference formula.
		for (int i=0;i<m_dim;i++) fjac[i][j] = (m_vv[i]-m_fvec[i])/h; 
	}
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
template <class F> 
bool
RootFinderND_t<F>::line_search(double*	xold,
								  double&	fold,
								  double*	g,
								  double*	p,
								  double*	x, // ret
								  double	stpmax,
								  double	x_precision,
								  double	f_precision)
{
	double alam2(0.), f2(0.), tmplam(0.);

	double sum(0.);
	for (int i=0;i<m_dim;i++) sum += p[i]*p[i];
	sum = sqrt(sum);

	//Scale if attempted step is too big.
	if (sum > stpmax) for (i=0;i<m_dim;i++) p[i] *= stpmax/sum; 
	
	double slope(0.);
	for(i=0;i<m_dim;i++) slope += g[i]*p[i];
    if (slope < 0.) ICMTHROW(ERR_INVALID_ARGUMENT,"Roundoff problem in lnsrch.");

	// Compute fmin.
	double test = 0.; 
	for (i=0;i<m_dim;i++) {
		double temp = fabs(p[i])/std::_cpp_max(fabs(m_xold[i]),1.0);
		if (temp > test) test=temp;
	}

	double alamin = x_precision/test;
	//Always try full Newton step .rst.
	double alam = 1.; 
	for (;;) { //Start of iteration loop.
		for (i=0;i<m_dim;i++) x[i] = m_xold[i] + alam*p[i];
		double f = fmin(x);
		if (alam < alamin) { //Convergence on Dx. For zero finding, the calling program shouldverify the convergence.
			for (i=0;i<m_dim;i++) x[i] = m_xold[i];
			fold = f;
			return true;
		}
		else if (f <= fold+f_precision*alam*slope) {
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
/*Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
indx[1..n] is an output vector that records the row permutation e.ected by the partial
pivoting; d is output as ±1 depending on whether the number of row interchanges was even
or odd, respectively. This routine is used in combination with lubksb to solve linear equations
or invert a matrix.*/
template <class F> 
void
RootFinderND_t<F>::ludcmp(double **a,
							  int n,
							  int *indx,
							  int& d)
{
	int imax(-1);
	double big(0.), dum, sum(0.), temp;
	
	d = 1; //No row interchanges yet.
	for (int i=0;i<n;i++) { //Loop over rows to get the implicit scaling information.
		big = 0.;
		for (int j=0;j<n;j++) if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.) //EP_THROWEX("Singular matrix in routine ludcmp");
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Singular matrix in routine ludcmp");	
		//No nonzero largest element.
		m_vv[i] = 1./big; //Save the scaling.
	}

	for (int j=0;j<n;j++) { //This is the loop over columns of Crout’s method.
		for (i=0;i<j;i++) { //This is equation (2.3.12) except for i = j.
			sum = a[i][j];
			for (int k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
		}
		big = 0.; //Initialize for the search for largest pivot element.
		for (i=j;i<n;i++) { //This is i = j of equation (2.3.12) and i = j+1. . .N of equation (2.3.13).
			sum = a[i][j];
			for (int k=0;k<j;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=m_vv[i]*fabs(sum)) >= big) { //Is the figure of merit for the pivot better than the best so far?
				big = dum;
				imax =i;
			}
		}
	
		if (j != imax) { //Do we need to interchange rows?
			for (int k=0;k<n;k++) { //Yes, do so...
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			d = -d; //...and change the parity of d.
			m_vv[imax] = m_vv[j]; //Also interchange the scale factor.
		}
		indx[j] = imax;
		if (a[j][j] == 0.) a[j][j] = 1e-18;
		//If the pivot element is zero the matrix is singular (at least to the precision of the
		//algorithm). For some applications on singular matrices, it is desirable to substitute TINY for zero.
		
		if (j != n-1) { //Now, finally, divide by the pivot element.
			dum = 1./(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	} //Go back for the next column in the reduction.

}
//----------------------------------------------------------------------------
/*RootFinderND_t<F>::lubksb(float **a, int n, int *indx, float b[])
Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix
A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
B, and returns with the solution vector X. a, n, and indx are not modi.ed by this routine
and can be left in place for successive calls with di.erent right-hand sides b. This routine takes
into account the possibility that b will begin with many zero elements, so it is e.cient for use
in matrix inversion.*/
template <class F>
void
RootFinderND_t<F>::lubksb(double **a,
							  int n,
							  int* indx,
							  double* b)
{
	int ii=-1;
	double sum(0.);
	for (int i=0;i<n;i++) {
		//When ii is set to a positive value, it will become the index of the first nonvanishing element of b. Wenow
		//do the forward substitution, equation (2.3.6). The only new wrinkle is to unscramble the permutation as we go.
		int ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii!=-1) for (int j = ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i; //A nonzero element was encountered, so from now on we
							// will have to do the sums in the loop above.
		b[i]=sum;
	}

	//Now we do the backsubstitution, equation (2.3.7).
	for (i=n-1;i>=0;i--) { 
		sum = b[i];
		for (int j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i] = sum/a[i][i]; //Store a component of the solution vector X.
	}
}
//----------------------------------------------------------------------------
//Given an initial guess x[1..n] for a root in n dimensions, take ntrial Newton-Raphson steps
//to improve the root. Stop if the root converges in either summed absolute variable increments
//tolx or summed absolute function values tolf.
template <class F>
double
RootFinderND_t<F>::NewtonMultiDim(double* x,
								int n,
								int ntrial,
								double tolx,
								double tolf)
{
	init(n);

	#ifdef _DEBUG
	char TAB[500];
	sprintf(TAB,"c:\\temp\\Optimisation.txt");
	FILE *stream = fopen(TAB, "w+");
	fprintf(stream,"BEGIN ---------------------------------------------\n");
	#endif 

	for (int k=0;k<ntrial;k++) {

		//User function supplies function values at x in fvec and Jacobian matrix in fjac. 
		m_f(x, m_fvec);
		fdjac(x, m_fvec, m_fjac);

		//Check function convergence.
		double errf = 0.;
		for (int i=0;i<m_dim;i++) errf += fabs(m_fvec[i]); 

		#ifdef _DEBUG
		fprintf(stream,"(errf,tolf) = (%f,%f)\n",errf,tolf);
		#endif 

		if (errf <= tolf) {
			//ICMTHROW(ERR_INVALID_ARGUMENT,"RootFinderND::NewtonMultiDim : Convergence achieved in " << k << " iterations.");
			#ifdef _DEBUG
			fprintf(stream,"(errf <= tolf) = (%f <= %f)\n",errf,tolf);
			fprintf(stream,"TERMINATE -----------------------------------------\n");
			fclose(stream);
			#endif 
			return sqrt(2.*fmin(x));
		}

		//Right-hand side of linear equations.
		for (i=0;i<m_dim;i++) m_p[i] = -m_fvec[i]; 
		
		//Solve linear equations using LU decomposition.
		int d;
		ludcmp(m_fjac, m_dim, m_indx, d); 
		lubksb(m_fjac, m_dim, m_indx, m_p);

		//Check root convergence.
		double errx = 0.; 
		for (i=0;i<m_dim;i++) { //Update solution.
			errx += fabs(m_p[i]);
			x[i] += m_p[i];
		}

		if (errx <= tolx) {
			#ifdef _DEBUG
			fprintf(stream,"(errx <= tolx) = (%f <= %f)\n",errx,tolx);
			fprintf(stream,"TERMINATE -----------------------------------------\n");
			fclose(stream);
			#endif 
			//ICMTHROW(ERR_INVALID_ARGUMENT,"RootFinderND::NewtonMultiDim : Convergence achieved in " << k << " iterations.");
			return sqrt(2.*fmin(x));
		}
		
	}

	#ifdef _DEBUG
	fprintf(stream,"Maximum number of iterations exceeded in NewtonMultiDim\n");
	fprintf(stream,"TERMINATE -----------------------------------------\n");
	fclose(stream);
	#endif 
	
	//ICMTHROW(ERR_INVALID_ARGUMENT,"Maximum number of iterations exceeded in NewtonMultiDim");
	return 1.e20;
}
//----------------------------------------------------------------------------
//Given an initial guess x[1..n] for a root in n dimensions, find the root by a globally convergent
//Newton’s method.
//The vector of functions to be zeroed, called m_fvec[1..n] in the routine below, is returned by the user-supplied routine m_f.
//The output quantity check is false (0) on a normal return and true (1) if the routine has converged to a local
//minimum of the function fmin defined below.
//In this case try restarting from a different initial guess.
template <class F>
double
RootFinderND_t<F>::GloballyConvergentNewton(double* x,
												 int n,
												 bool& check,
												 int nbtrial,
												 double x_precision,
												 double f_precision,
												 double step_max,
												 double tol_min)
{
	init(n);

	#ifdef _DEBUG
	char TAB[500];
	sprintf(TAB,"c:\\temp\\Optimisation.txt");
	FILE *stream = fopen(TAB, "w+");
	fprintf(stream,"BEGIN ---------------------------------------------\n");
	#endif 

	double fold(0.), temp(0.);
	
	double f = fmin(x); //fvec is also computed by this call.
	
	//Test for initial guess being a root. Use more stringent test than simply TOLF.
	double test = 0.; 
	for (int i=0;i<m_dim;i++) if (fabs(m_fvec[i]) > test) test = fabs(m_fvec[i]);
	
	if (test < 0.01*f_precision) {
		check = false;
	#ifdef _DEBUG
		fprintf(stream,"(test <= 0.01*f_precision) = (%f <= %f)\n",test,0.01*f_precision);
		fprintf(stream,"Convergence achieved in 0 iterations.\n");
		fprintf(stream,"TERMINATE -----------------------------------------\n");
		fclose(stream);
		#endif 

		return sqrt(2.*fmin(x));
	}
	double sum(0.);
	for (i=0;i<m_dim;i++) sum += x[i]*x[i];
	
	//Calculate stpmax for line searches.
	double stpmax = step_max * std::_cpp_max (sqrt(sum), (double)m_dim);
	//Start of iteration loop.
	for (int its=1;its<=nbtrial;its++) { 

		fdjac(x, m_fvec, m_fjac); //If analytic Jacobian is available, you can replace the routine fdjac below with your own routine.
		for (i=0;i<m_dim;i++) { //Compute gradf for the line search.
			m_g[i] = 0.;
			for (int j=0;j<m_dim;j++) m_g[i] += m_fjac[j][i]*m_fvec[j];
		}

		//Store x and f
		for (i=0;i<m_dim;i++) m_xold[i] = x[i]; 
		fold = f;

		for (i=0;i<m_dim;i++) m_p[i] = -m_fvec[i]; //Right-hand side for linear equations.
		
		//Solve linear equations by LU decomposition.
		int d;
		ludcmp(m_fjac, m_dim, m_indx, d); 
		lubksb(m_fjac, m_dim, m_indx, m_p);

		//lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.
		//lnsrch(m_dim, xold, fold, g, m_p, x, &f, stpmax, check, fmin);
		check = line_search(m_xold, fold, m_g, m_p, x, stpmax, x_precision, f_precision);

		//Test for convergence on function values.
		test = 0.; 
		for (i=0;i<m_dim;i++) if (fabs(m_fvec[i]) > test) test = fabs(m_fvec[i]);
		if (test < f_precision) {
			#ifdef _DEBUG
			fprintf(stream,"(test <= f_precision) = (%f <= %f)\n",test,f_precision);
			fprintf(stream,"Convergence achieved in %i iterations.\n",its);
			fprintf(stream,"TERMINATE -----------------------------------------\n");
			fclose(stream);
			#endif 
			check = false;
	
			return sqrt(2.*fmin(x));
		}
		if (check) { //Check for gradient of f zero, i.e., spurious convergence.
			test = 0.;
			double den = std::_cpp_max (f, 0.5*m_dim);
			for (i=0;i<m_dim;i++) {
				temp = fabs(m_g[i]) * std::_cpp_max (fabs(x[i]), 1.)/den;
				if (temp > test) test=temp;
			}
			check = (test < tol_min ? true : false);

			#ifdef _DEBUG
			fprintf(stream,"(test <= f_precision) = (%f <= %f)\n",test,f_precision);
			fprintf(stream,"Convergence achieved in %i iterations.\n",its);
			fprintf(stream,"TERMINATE -----------------------------------------\n");
			fclose(stream);
			#endif 

			return sqrt(2.*fmin(x));
		}

		//Test for convergence on delta_x.
		test = 0.; 
		for (i=0;i<m_dim;i++) {
			temp = (fabs(x[i]-m_xold[i])) / std::_cpp_max (fabs(x[i]), 1.);
			if (temp > test) test=temp;
		}
		if (test < x_precision) {
			#ifdef _DEBUG
			fprintf(stream,"(test <= f_precision) = (%f <= %f)\n",test,f_precision);
			fprintf(stream,"Convergence achieved in %i iterations.\n",its);
			fprintf(stream,"TERMINATE -----------------------------------------\n");
			fclose(stream);
			#endif 
			return sqrt(2.*fmin(x));
		}
	}

	#ifdef _DEBUG
	fprintf(stream,"Maximum number of iterations exceeded in NewtonMultiDim\n");
	fprintf(stream,"TERMINATE -----------------------------------------\n");
	fclose(stream);
	#endif 

	return 1.e20;
}
//----------------------------------------------------------------------------
template <class F>
bool
RootFinderND_t<F>::test(std::string& errStr)
{
	bool ret(true);
	
	errStr = "Testing RootFinderND";
	
	//ATest a;
	
	/*std::vector<double> init;
	init.push_back(-30.); init.push_back(1.);
	val = RootFinder1D(f1::mem_call(&ATest::f, a), f1::mem_call(&ATest::df, a)).NewtonRaphson ( 0., init, -4., 5., 20, 0.000001, 0.000001, false);	
	SDMTEST((val>2.9999)&&(val<3.0001))
	*/
	return ret;
}


#endif	// _F1_RootFinder1D_H_


